// Batch worker: loads the wasm module once, folds one sequence per message.
// The main thread passes the module to load (?m=tornadofold[-simd].js) so the
// worker matches the page's SIMD-support decision.
const module = new URLSearchParams(self.location.search).get("m") || "tornadofold.js";
importScripts(module);

const ready = tornadofold();

self.onmessage = (e) => {
  const { id, seq } = e.data;
  ready.then((mod) => {
    const res = mod.foldSeq(seq);
    const tab = res.lastIndexOf("\t");
    self.postMessage({ id, db: res.slice(0, tab), mfe: parseInt(res.slice(tab + 1), 10) });
  });
};
