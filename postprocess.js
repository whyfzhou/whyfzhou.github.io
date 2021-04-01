function sample(arcs, nsample) {
  let tsw = arcs.map((e) => e.dt).reduce((acc, e) => acc + e, 0);
  let t = [...Array(nsample).keys()].map((e) => (e / nsample) * tsw);
  let ir = Array(nsample);
  let vr = Array(nsample);
  let im = Array(nsample);
  let vhb = Array(nsample);
  let ta = 0;
  let tb = 0;
  for (let i = 0; i < arcs.length; ++i) {
    tb = ta + arcs[i].dt;
    let piecewiseIndices = [...t.keys()].filter((e, k) => ta <= t[k] && t[k] <= tb);
    let start = piecewiseIndices[0];
    let count = piecewiseIndices.length;
    if (Array.isArray(piecewiseIndices) && count) {
      let i0 = arcs[i].i0;
      // let v0 = arcs[i].v0;
      let r = arcs[i].r;
      let phi = arcs[i].phi;
      let w = arcs[i].w;
      let z = arcs[i].z;
      let vg = arcs[i].vg;
      let km = arcs[i].km;
      let hb = arcs[i].hb;
      ir.splice(start, count, ...t.slice(start, start + count).map((tau) => (r * Math.cos(w * (tau - ta) + phi)) / z));
      vr.splice(start, count, ...t.slice(start, start + count).map((tau) => vg + r * Math.sin(w * (tau - ta) + phi)));
      if (Math.abs(km) > 1e-9) {
        im.splice(start, count, ...t.slice(start, start + count).map((tau) => i0 + km * (tau - ta)));
      } else {
        im.splice(start, count, ...ir.slice(start, start + count));
      }
      vhb.splice(start, count, ...Array(count).fill(hb));
    }
    ta = tb;
  }
  return { t: t, ir: ir, vr: vr, im: im, vhb: vhb };
}

function trapz(x, t) {
  let integral = x
    .slice(1)
    .map((e, k) => ((e + x[k]) * (t[k + 1] - t[k])) / 2)
    .reduce((acc, e) => acc + e, 0);
  return integral;
}

function evaluateSampled(arcs, sampled) {
  let tsw = arcs.map((e) => e.dt).reduce((acc, e) => acc + e, 0);
  let dutyHighSide = arcs[2].dt / tsw;
  let dutyDiode = arcs[0].dt / (arcs[0].dt + arcs[1].dt + arcs[3].dt);

  let iDCOut =
    trapz(
      sampled.im.map((e, k) => e - sampled.ir[k]),
      sampled.t
    ) / tsw;

  let iRMSPri =
    (trapz(
      sampled.ir.map((e) => e ** 2),
      sampled.t
    ) /
      tsw) **
    0.5;
  let iPPeakPri = Math.max(...sampled.ir);
  let iNPeakPri = Math.min(...sampled.ir);

  let vMeanRes = trapz(sampled.vr, sampled.t) / tsw;
  let vACRMSRes =
    (trapz(
      sampled.vr.map((e) => (e - vMeanRes) ** 2),
      sampled.t
    ) /
      tsw) **
    0.5;
  let vRMSPri =
    (trapz(
      sampled.vhb.map((e, k) => (e - sampled.vr[k]) ** 2),
      sampled.t
    ) /
      tsw) **
    0.5;

  let iRMSOut =
    (trapz(
      sampled.im.map((e, k) => (e - sampled.ir[k]) ** 2),
      sampled.t
    ) /
      tsw) **
    0.5;
  let iPPeakOut = Math.max(...sampled.im.map((e, k) => e - sampled.ir[k]));

  let psiPri =
    trapz(
      sampled.vhb.map((e, k) => Math.abs(e - sampled.vr[k])),
      sampled.t
    ) / 4;

  let qLowHigh = trapz(
    sampled.ir.map((e, k) => (e < 0 && sampled.vhb[k] > 0 ? e : 0)),
    sampled.t
  );

  return {
    fsw: 1 / tsw,
    dutyHighSide: dutyHighSide,
    dutyDiode: dutyDiode,
    iDCOut: iDCOut,
    iRMSPri: iRMSPri,
    iPPeakPri: iPPeakPri,
    iNPeakPri: iNPeakPri,
    vMeanRes: vMeanRes,
    vACRMSRes: vACRMSRes,
    vRMSPri: vRMSPri,
    iRMSOut: iRMSOut,
    iPPeakOut: iPPeakOut,
    psiPri: psiPri,
    qLowHigh: qLowHigh,
  };
}
