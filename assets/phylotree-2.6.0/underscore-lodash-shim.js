// phylotree@2.6.0's UMD bundle (dist/phylotree.js) declares BOTH
// underscore and lodash as external browser globals - but its own build
// config never gave them distinct global names, so both defaulted to
// "_"; rollup only avoided a *local parameter name* collision inside its
// own factory function by suffixing the second one to "_$1" - that's
// purely a local variable name, not what it actually reads off the
// global object. Verified directly in the unminified bundle: the real
// UMD header is
//   factory(global.phylotree = ..., global._, global._$1)
// i.e. it genuinely expects a global literally named "_$1" for lodash -
// nothing ever defines that on its own. This is a real upstream
// packaging bug (both externals collapsed onto the conventional "_"
// name), not something fixable by loading things in a different order
// alone.
//
// Load underscore, then lodash (both target window._ when loaded via a
// plain <script> tag, so lodash's load overwrites underscore's), then
// this shim - swaps them into the two distinct globals phylotree.js
// actually reads.
window._$1 = window._; // currently lodash (loaded last)
window._ = window.__phylotreeUnderscore; // restore underscore
