function simdSupported() {
  try {
    return WebAssembly.validate(new Uint8Array([
      0, 97, 115, 109, 1, 0, 0, 0, 1, 5, 1, 96, 0, 1, 123, 3,
      2, 1, 0, 10, 10, 1, 8, 0, 65, 0, 253, 15, 253, 98, 11]));
  } catch (e) {
    return false;
  }
}

var modulePromise = null;

/** Instantiate the wasm module */
function load() {
  if (modulePromise === null) {
    var factory = simdSupported()
      ? require("./tornadofold-simd.js")
      : require("./tornadofold.js");
    modulePromise = factory();
  }
  return modulePromise;
}

/**
 * Fold an RNA sequence. Whitespace is removed, T is read as U and other
 * characters are treated as N, so positions line up with the input.
 *
 * @param {string} sequence
 * @returns {Promise<{ structure: string, mfe: number }>} dot-bracket structure
 *          and its free energy in kcal/mol
 */
function fold(sequence) {
  return load().then(function (mod) {
    var result = mod.foldSeq(sequence);
    var tab = result.lastIndexOf("\t");
    return {
      structure: result.slice(0, tab),
      mfe: parseInt(result.slice(tab + 1), 10) / 100,
    };
  });
}

module.exports = { fold: fold, load: load, simdSupported: simdSupported };
