ckt = { lr: 30e-6, lm: 1000e-6, cr: 47e-9, chb: 500e-12, vbus: 410, vout: 200 };
let states = findSteadyState(2e-6, ckt, 500e-9, 100e3);
data = sample(states, 20000);

window.onload = test;
window.app = [];
window.app.resultIndices = [];

function test() {
  document.querySelector("#test-output").innerText = "I'm here.";
  // document.querySelector("#test-output").innerText = `dv = ${dv}`;
  drawSingleOP(data);
}
