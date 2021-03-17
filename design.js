window.onload = init;

// ----------------------------------------------------------------------------

function init() {
  document.querySelector("#action-analyze").onclick = analyze;
}

// ----------------------------------------------------------------------------

function analyze() {
  let vbus = document.querySelector("#vbus").valueAsNumber;
  let vload = document.querySelector("#vload").valueAsNumber;
  let pload = document.querySelector("#pload").valueAsNumber;

  let t12 = document.querySelector("#t12").valueAsNumber * 1e-9;
  let nps = document.querySelector("#nps").valueAsNumber;
  let cr = document.querySelector("#cr").valueAsNumber * 1e-9;

  let lrmin = document.querySelector("#lrmin").valueAsNumber * 1e-6;
  let lrmax = document.querySelector("#lrmax").valueAsNumber * 1e-6;
  let lpri1 = document.querySelector("#lpri1").valueAsNumber * 1e-6;
  let lpri2 = document.querySelector("#lpri2").valueAsNumber * 1e-6;
  let lpri3 = document.querySelector("#lpri3").valueAsNumber * 1e-6;

  if ([vbus, vload, pload, t12, nps, cr, lrmin, lrmax].some((e) => isNaN(e))) {
    alert("Input Error: missing input.");
    return;
  }
  if ([lpri1, lpri2, lpri3].every((e) => isNaN(e))) {
    alert("Input Error: should at least give one Lpri.");
    return;
  }

  let vout = vload * nps;

  const n = 25;
  let lrs = [...Array(n).keys()].map((e) => lrmin + ((lrmax - lrmin) * e) / (n - 1));
  let lpris = [lpri1, lpri2, lpri3].filter((e) => !isNaN(e));

  let y00 = [];
  let y01 = [];
  let y10 = [];
  let y11 = [];
  for (let lpri of lpris) {
    function evaluateLr(lr) {
      let lm = lpri - lr;
      let iRMSPri;
      let psiPri;
      let qLowHigh;
      let fsw;
      try {
        let steadyState = steady_state_pout(pload, lr, lm, cr, vbus, vout, t12);
        let sampledData = sample(steadyState, 2000);
        let measurement = evaluateSampled(steadyState, sampledData);
        iRMSPri = measurement.iRMSPri;
        psiPri = measurement.psiPri;
        qLowHigh = measurement.qLowHigh;
        fsw = measurement.fsw;
      } catch (e) {
        iRMSPri = psiPri = qLowHigh = fsw = NaN;
      } finally {
        return [iRMSPri, psiPri, qLowHigh, fsw];
      }
    }
    let evaluated = lrs.map(evaluateLr);
    y00.push(evaluated.map((e) => e[0]));
    y01.push(evaluated.map((e) => e[1]));
    y10.push(evaluated.map((e) => e[2]));
    y11.push(evaluated.map((e) => e[3]));
  }
  drawDesign(lrs, [y00, y01, y10, y11], lpris);
}
