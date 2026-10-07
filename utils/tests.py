import time
from unittest.mock import patch

from django import forms
from django.core import signing
from django.test import TestCase

from utils.powcaptcha import (
    CHALLENGE_SALT, NONCE_FIELD, SOLVE_CONCURRENCY, TOKEN_FIELD,
    PowCaptchaFormMixin, make_challenge, pow_proof_for_form, solve_challenge)


class _Form(PowCaptchaFormMixin, forms.Form):
    """A minimal form with nothing but the mixin's own fields, for
    testing the mixin in isolation from either real form it's mixed
    into (surveys.forms.FeedbackForm, account.forms.AccountCreationForm)."""


class MakeChallengeTest(TestCase):

    def test_returns_a_salt_matching_the_signed_token(self):
        salt, difficulty, token = make_challenge(difficulty=4)
        decoded = signing.loads(token, salt=CHALLENGE_SALT)
        self.assertEqual(decoded, {'salt': salt, 'difficulty': 4})

    def test_defaults_to_the_configured_difficulty(self):
        with self.settings(NGPHYLO_POW_CAPTCHA_DIFFICULTY=7):
            _, difficulty, _ = make_challenge()
        self.assertEqual(difficulty, 7)

    def test_two_challenges_get_different_salts(self):
        salt1, _, _ = make_challenge()
        salt2, _, _ = make_challenge()
        self.assertNotEqual(salt1, salt2)


class SolveChallengeTest(TestCase):

    def test_solution_actually_satisfies_the_difficulty(self):
        salt, difficulty, _ = make_challenge(difficulty=4)
        nonce = solve_challenge(salt, difficulty)
        import hashlib
        digest = hashlib.sha256(('%s%s' % (salt, nonce)).encode()).hexdigest()
        self.assertTrue(digest.startswith('0' * difficulty))


class PowCaptchaFormMixinTest(TestCase):
    """
    utils.powcaptcha.PowCaptchaFormMixin, exercised directly through a
    minimal form rather than through surveys.forms.FeedbackForm/
    account.forms.AccountCreationForm (both covered by their own
    app's tests) - this is the one place the mixin's own behavior is
    tested in isolation.
    """

    def test_unbound_form_exposes_a_solvable_challenge(self):
        form = _Form()
        token = form.initial[TOKEN_FIELD]
        self.assertTrue(token)
        attrs = form.fields[TOKEN_FIELD].widget.attrs
        self.assertIn('data-pow-salt', attrs)
        self.assertIn('data-pow-difficulty', attrs)
        decoded = signing.loads(token, salt=CHALLENGE_SALT)
        self.assertEqual(decoded['salt'], attrs['data-pow-salt'])
        self.assertEqual(decoded['difficulty'], attrs['data-pow-difficulty'])

    def test_a_correct_solution_validates(self):
        unbound = _Form()
        proof = pow_proof_for_form(unbound)
        form = _Form(data=proof)
        self.assertTrue(form.is_valid(), form.errors)

    def test_a_wrong_nonce_fails_validation(self):
        unbound = _Form()
        proof = pow_proof_for_form(unbound)
        proof[NONCE_FIELD] = str(int(proof[NONCE_FIELD]) + 1)
        form = _Form(data=proof)
        self.assertFalse(form.is_valid())
        self.assertIn(NONCE_FIELD, form.errors)

    def test_a_missing_nonce_fails_validation_with_a_distinct_message(self):
        unbound = _Form()
        token = unbound.initial[TOKEN_FIELD]
        form = _Form(data={TOKEN_FIELD: token, NONCE_FIELD: ''})
        self.assertFalse(form.is_valid())
        self.assertIn('wait', form.errors[NONCE_FIELD][0])

    def test_a_tampered_token_fails_with_a_reload_message(self):
        unbound = _Form()
        token = unbound.initial[TOKEN_FIELD] + 'x'
        form = _Form(data={TOKEN_FIELD: token, NONCE_FIELD: '0'})
        self.assertFalse(form.is_valid())
        self.assertIn('reload', form.errors[TOKEN_FIELD][0])
        # The nonce field shouldn't also raise its own, redundant error
        # once the token itself is already invalid.
        self.assertNotIn(NONCE_FIELD, form.errors)

    def test_an_expired_token_fails_validation(self):
        # signing.dumps() stamps its own creation time, so the clock
        # has to be backdated at creation (not just the payload), same
        # technique as surveys.tests.FeedbackAntiSpamTest's own
        # expired-token test.
        with patch('django.core.signing.time.time',
                   return_value=time.time() - 60 * 60):
            salt, difficulty, token = make_challenge()
        nonce = solve_challenge(salt, difficulty)
        form = _Form(data={TOKEN_FIELD: token, NONCE_FIELD: str(nonce)})
        self.assertFalse(form.is_valid())
        self.assertIn('reload', form.errors[TOKEN_FIELD][0])

    def test_a_bound_redisplay_keeps_an_already_solved_nonce_valid(self):
        # Simulates the form being redisplayed after an unrelated
        # validation error elsewhere - the already-submitted token/nonce
        # pair must still validate on its own, without needing a fresh
        # challenge or another solve.
        unbound = _Form()
        proof = pow_proof_for_form(unbound)
        first = _Form(data=proof)
        self.assertTrue(first.is_valid())
        second = _Form(data=proof)
        self.assertTrue(second.is_valid())

    def test_layout_html_embeds_the_configured_concurrency(self):
        form = _Form()
        html = form.pow_layout_html()
        self.assertIn('var concurrency = %d;' % SOLVE_CONCURRENCY, html)
        self.assertIn('id_pow_token', html)
        self.assertIn('id_pow_nonce', html)
