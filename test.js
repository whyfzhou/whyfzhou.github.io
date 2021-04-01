let ckt = { lr: 30e-6, lm: 1000e-6, cr: 47e-9, chb: 500e-12, vbus: 410, vout: 200 };

// let states = simulate(201, ckt, [500e-9, 500e-9]).slice(-1)[0];

// let states = findSteadyState(500e-9, ckt, 500e-9, 100e3);
// let perf = evaluateSwitchingPeriod(states);

let tf, states, evaluated, pmax;
[tf, states, evaluated, pmax] = evaluateOperatingPoint(40, ckt, 500e-9, 100e3);

let data = sample(states, 20000);
let perfsampled = evaluateSampled(states, data);

window.onload = test;
window.app = [];
window.app.resultIndices = [];

function test() {
  document.querySelector("#test-output").innerText = "I'm here.";
  // document.querySelector("#test-output").innerText = `dv = ${dv}`;
  drawSingleOP(data);
}
