"""
Self-hosted proof-of-work CAPTCHA - a replacement for django-simple-captcha's
distorted-letter image (see CLAUDE.md's "Contact form anti-spam layers": "the
exactly kind of challenge automated solvers handle routinely"). Instead of
asking a human to read an image, the browser has to do real CPU work (brute-
force search a SHA-256 preimage) before the form can be submitted at all.

No external service, no API key, no new Python dependency: the `altcha`
PyPI package needs Python >=3.9 (this project is pinned to 3.8 - see
CLAUDE.md's "biopython==1.70" note), and its current wire protocol (nonce +
KDF-derived key prefixes, canonical JSON) is intricate enough that hand-
matching it without the library wasn't worth the risk. This instead
implements the same idea - and the same simple "hash of salt+number must
start with N zero hex digits" challenge ALTCHA itself used before that -
end to end in code this project owns: `django.core.signing` (already used
for the Workspace permalink and AntiSpamFormMixin's form-age token) issues
and verifies a signed, expiring challenge, plain `hashlib` checks the
proof, and a small vendor-free vanilla-JS snippet (native
`crypto.subtle.digest`, no vendored library) solves it in the browser.

Honest limitation, same one every client-side PoW captcha has (ALTCHA
included): a targeted attacker who re-implements the hash loop natively
solves it far faster than a browser can. What this defends against is
mass, naive scripted spam (the kind this project has actually seen - see
CLAUDE.md's qq.com-registration discussion) by giving every submission a
real, non-zero CPU cost, not a sophisticated, individually-targeted one.

Difficulty was picked by actually measuring solve times in headless
Chrome (native `crypto.subtle.digest`, 16 concurrent in-flight hashes),
not guessed: difficulty 5 (20 bits, ~1e6 expected attempts) solved in
well under a second up to a few seconds across repeated trials - an
acceptable "please wait a moment" cost, not a hang.
"""
import hashlib
import secrets

from django.conf import settings
from django.core import signing
from django.core.exceptions import ValidationError
from django.forms import CharField, HiddenInput
from django.utils.html import format_html

CHALLENGE_SALT = 'utils.powcaptcha.challenge'
# An unsolved challenge sitting open is only ever meant to be solved right
# after the page that issued it loads - deliberately much shorter than
# antispam.MAX_FORM_AGE_SECONDS, which bounds a filled-in form, not an
# unsolved PoW challenge.
CHALLENGE_MAX_AGE_SECONDS = 10 * 60

TOKEN_FIELD = 'pow_token'
NONCE_FIELD = 'pow_nonce'

# How many concurrent in-flight crypto.subtle.digest() calls the browser
# keeps outstanding while searching - see the module docstring for how
# this was calibrated alongside the difficulty default below.
SOLVE_CONCURRENCY = 16


def _difficulty():
    return getattr(settings, 'NGPHYLO_POW_CAPTCHA_DIFFICULTY', 5)


def _digest_hex(salt, number):
    return hashlib.sha256(('%s%s' % (salt, number)).encode('utf-8')).hexdigest()


def make_challenge(difficulty=None):
    """Returns (salt, difficulty, signed_token) for a fresh challenge."""
    if difficulty is None:
        difficulty = _difficulty()
    salt = secrets.token_hex(16)
    token = signing.dumps({'salt': salt, 'difficulty': difficulty},
                           salt=CHALLENGE_SALT)
    return salt, difficulty, token


def solve_challenge(salt, difficulty):
    """
    A plain-Python reference solver - what the browser's inline script
    does with crypto.subtle.digest(), done here with hashlib instead.
    Not used by the real form flow (the browser solves it), only by
    tests that need a genuinely valid nonce to submit, and as a sanity
    check that the server and client sides agree on the same hash.
    """
    number = 0
    while not _digest_hex(salt, number).startswith('0' * difficulty):
        number += 1
    return number


def pow_proof_for_form(form):
    """
    Test helper: given a freshly-instantiated (unbound) form that mixes
    in PowCaptchaFormMixin, solves its challenge and returns the
    {pow_token: ..., pow_nonce: ...} dict to merge into POST data for a
    submission that should pass validation.
    """
    token = form.initial[TOKEN_FIELD]
    attrs = form.fields[TOKEN_FIELD].widget.attrs
    nonce = solve_challenge(attrs['data-pow-salt'], attrs['data-pow-difficulty'])
    return {TOKEN_FIELD: token, NONCE_FIELD: str(nonce)}


class PowCaptchaFormMixin:
    """
    Adds the hidden token/nonce fields and their validation to a form -
    same shape as utils.antispam.AntiSpamFormMixin (which this is meant
    to be combined with, not replace): must come before the form's own
    base class in the MRO so its __init__ runs and its clean_* hooks are
    found.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields[TOKEN_FIELD] = CharField(widget=HiddenInput)
        self.fields[NONCE_FIELD] = CharField(
            required=False, max_length=24, widget=HiddenInput)
        if not self.is_bound:
            salt, difficulty, token = make_challenge()
            self.initial[TOKEN_FIELD] = token
            # Not secret - only the signed token above is what verification
            # actually trusts. Exposed here purely so the inline script
            # below knows what to solve.
            self.fields[TOKEN_FIELD].widget.attrs.update({
                'data-pow-salt': salt,
                'data-pow-difficulty': difficulty,
            })

    def pow_layout_html(self):
        """
        A crispy-forms HTML() snippet: the visible "verifying" indicator
        plus the solver script, parameterized with this field's own
        auto id (Django's default 'id_<name>', unchanged here) so it
        finds the right inputs without assuming only one form per page.
        """
        token_id = self[TOKEN_FIELD].auto_id
        nonce_id = self[NONCE_FIELD].auto_id
        return format_html('''
<div class="pow-captcha-status" id="{token_id}-status">Verifying your browser&hellip;</div>
<script>
(function() {{
    var tokenInput = document.getElementById("{token_id}");
    var nonceInput = document.getElementById("{nonce_id}");
    var statusEl = document.getElementById("{token_id}-status");
    if (!tokenInput || !nonceInput || !statusEl) {{ return; }}
    // Already solved (e.g. the form was redisplayed after an unrelated
    // validation error) - nothing to do.
    if (nonceInput.value) {{
        statusEl.textContent = "Verified.";
        statusEl.className = "pow-captcha-status pow-captcha-status-ok";
        return;
    }}
    var form = tokenInput.closest("form");
    var submitBtn = form ? form.querySelector('[type=submit]') : null;
    if (submitBtn) {{ submitBtn.disabled = true; }}
    if (!window.crypto || !window.crypto.subtle) {{
        statusEl.textContent =
            "Your browser does not support the verification this form needs - please use a modern browser.";
        statusEl.className = "pow-captcha-status pow-captcha-status-error";
        return;
    }}
    var salt = tokenInput.getAttribute("data-pow-salt");
    var prefix = new Array(parseInt(tokenInput.getAttribute("data-pow-difficulty"), 10) + 1).join("0");
    var concurrency = {concurrency};
    // A hard ceiling only meant to catch something fundamentally broken
    // (not the normal case - see the server-side module docstring for
    // how the default difficulty was measured) rather than ever hang
    // the tab indefinitely.
    var maxAttempts = 50000000;
    function sha256Hex(s) {{
        return crypto.subtle.digest("SHA-256", new TextEncoder().encode(s)).then(function(buf) {{
            return Array.prototype.map.call(new Uint8Array(buf), function(b) {{
                return b.toString(16).padStart(2, "0");
            }}).join("");
        }});
    }}
    function solveFrom(nonce) {{
        var batch = [];
        for (var i = 0; i < concurrency; i++) {{ batch.push(nonce + i); }}
        return Promise.all(batch.map(function(n) {{ return sha256Hex(salt + n); }})).then(function(hashes) {{
            for (var i = 0; i < hashes.length; i++) {{
                if (hashes[i].indexOf(prefix) === 0) {{ return batch[i]; }}
            }}
            if (nonce > maxAttempts) {{ return null; }}
            return solveFrom(nonce + concurrency);
        }});
    }}
    solveFrom(0).then(function(solution) {{
        if (solution === null) {{
            statusEl.textContent = "Could not verify your browser - please reload the page.";
            statusEl.className = "pow-captcha-status pow-captcha-status-error";
            return;
        }}
        nonceInput.value = String(solution);
        statusEl.textContent = "Verified.";
        statusEl.className = "pow-captcha-status pow-captcha-status-ok";
        if (submitBtn) {{ submitBtn.disabled = false; }}
    }});
}})();
</script>
''', token_id=token_id, nonce_id=nonce_id, concurrency=SOLVE_CONCURRENCY)

    def clean_pow_token(self):
        token = self.cleaned_data.get(TOKEN_FIELD, '')
        try:
            self._pow_challenge = signing.loads(
                token, salt=CHALLENGE_SALT, max_age=CHALLENGE_MAX_AGE_SECONDS)
        except signing.BadSignature:
            raise ValidationError(
                "Your browser verification expired - please reload the "
                "page and try again.")
        return token

    def clean_pow_nonce(self):
        nonce = self.cleaned_data.get(NONCE_FIELD, '')
        challenge = getattr(self, '_pow_challenge', None)
        if challenge is None:
            # clean_pow_token() already raised for this - no need for a
            # second, redundant error.
            return nonce
        if not nonce:
            raise ValidationError(
                "Please wait for the browser verification to finish "
                "(or enable JavaScript) before submitting.")
        if not _digest_hex(challenge['salt'], nonce) \
                .startswith('0' * challenge['difficulty']):
            raise ValidationError(
                "Browser verification failed - please reload the page "
                "and try again.")
        return nonce
