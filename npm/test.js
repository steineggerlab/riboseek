const fs = require("fs");
const path = require("path");

const GOLDEN = path.join(__dirname, "..", "tests", "data", "archiveii_smoke.tsv");

function foldWith(mod, sequence) {
  const result = mod.foldSeq(sequence);
  const tab = result.lastIndexOf("\t");
  return {
    structure: result.slice(0, tab),
    mfe: parseInt(result.slice(tab + 1), 10) / 100,
  };
}

async function checkBuild(label, factory, rows) {
  const mod = await factory();
  let fails = 0;
  for (const [name, seq, , mfe, db] of rows) {
    const got = foldWith(mod, seq);
    const gotMfe = got.mfe.toFixed(2);
    if (gotMfe !== mfe || got.structure !== db) {
      fails++;
      console.log(`FAIL ${label} ${name} (len ${seq.length}): MFE ${gotMfe} vs ${mfe}; ` +
        `structure ${got.structure === db ? "ok" : "DIFFERS"}`);
    }
  }
  console.log(`${label}: ${rows.length - fails}/${rows.length} sequences reproduce smoke MFE and structure`);

  const dna = foldWith(mod, "GGGGTTTTCCCC");
  if (dna.structure !== "((((....))))") {
    fails++;
    console.log(`FAIL ${label} dna input: ${dna.structure}`);
  }
  const empty = foldWith(mod, "");
  if (empty.structure !== "") {
    fails++;
    console.log(`FAIL ${label} empty input: ${JSON.stringify(empty.structure)}`);
  }
  return fails;
}

async function main() {
  const rows = fs.readFileSync(GOLDEN, "utf8").trim().split("\n").slice(1)
    .map((l) => l.split("\t"));

  // test both build
  let fails = await checkBuild("scalar", require("./tornadofold.js"), rows);
  fails += await checkBuild("simd", require("./tornadofold-simd.js"), rows);

  // and the package entry point
  const { fold, simdSupported } = require("./index.js");
  const trna = await fold("GGGCUAUUAGCUCAGUUGGUUAGAGCGCACCCCUGAUAAGGGUGAGGUCGCUGAUUCGAAUUCAGCAUAGCCCA");
  if (trna.structure.length !== 74 || !(trna.mfe < 0)) {
    fails++;
    console.log(`FAIL api: ${JSON.stringify(trna)}`);
  }
  console.log(`api: tRNA-Phe folds to ${trna.mfe} kcal/mol (simd: ${simdSupported()})`);

  if (fails) {
    process.exit(1);
  }
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
