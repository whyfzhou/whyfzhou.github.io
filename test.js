ckt = { lr: 30e-6, lm: 1000e-6, cr: 47e-9, chb: 500e-12, vbus: 410, vout: 200 };
let dv;
let states;
[dv, states] = simulate(200, ckt, [2e-6, 500e-9]);
data = sample(states);

window.onload = test;

function test() {
  document.querySelector("#test-output").innerText = "I'm here.";
  document.querySelector("#test-output").innerText = `dv = ${dv}`;
}
