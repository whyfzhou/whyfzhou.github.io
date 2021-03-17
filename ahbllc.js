function mod(a, b) {
  let c = a % b;
  if ((c < 0 && 0 < b) || (c > 0 && 0 > b)) {
    c += b;
  }
  return c;
}

function ridder(eq, x0, y0, x1, y1, xtol, ytol, maxIter) {
  let ans = Number.NaN;

  for (let numTries = 0; numTries < maxIter; ++numTries) {
    let x2 = (x0 + x1) / 2;
    let y2 = eq(x2);
    let s = (y2 * y2 - y0 * y1) ** 0.5;
    if (s === 0) {
      return ans;
    }

    let x3 = x2 + ((x2 - x0) * (y0 - y1 >= 0 ? 1 : -1) * y2) / s;
    if (Math.abs(x3 - ans) <= xtol) {
      return ans;
    }
    let y3 = eq(x3);
    ans = x3;
    if (Math.abs(y3) <= ytol) {
      return ans;
    }

    if (Math.sign(y2) !== Math.sign(y3)) {
      [x0, y0, x1, y1] = [x2, y2, x3, y3];
    } else if (Math.sign(y0) !== Math.sign(y3)) {
      [x1, y1] = [x3, y3];
    } else if (Math.sign(y1) !== Math.sign(y3)) {
      [x0, y0] = [x3, y3];
    } else {
      return Number.NaN;
    }
  }
  return Number.NaN;
}

function brent(eq, x0, y0, x1, y1, xtol, ytol, maxIter) {
  const eps = 7 / 3 - 4 / 3 - 1;

  let lastBisection = true;
  if (Math.abs(y0) < Math.abs(y1)) {
    [x0, x1, y0, y1] = [x1, x0, y1, y0];
  }
  let x2 = x0;
  let y2 = y0;
  let xs = x2;
  let ys = y2;
  let xd = x0;

  let numTries;
  for (numTries = 0; numTries < maxIter; ++numTries) {
    if (Math.abs(ys) <= ytol) {
      break;
    }
    if (Math.abs(y0 - y2) > 2 * eps && Math.abs(y1 - y2) > 2 * eps) {
      xs =
        (x0 * y1 * y2) / (y0 - y1) / (y0 - y2) +
        (x1 * y0 * y2) / (y1 - y0) / (y1 - y2) +
        (x2 * y0 * y1) / (y2 - y0) / (y2 - y1);
    } else {
      xs = x1 - (y1 * (x1 - x0)) / (y1 - y0);
    }

    let delta = Math.abs(2 * eps * x0);
    if (
      !(((3 * x0 + x1) / 4 < xs && xs < x1) || ((3 * x0 + x1) / 4 > xs && xs > x1)) ||
      (lastBisection && Math.abs(xs - x1) >= Math.abs(x1 - x2) / 2) ||
      (!lastBisection && Math.abs(xs - x1) >= Math.abs(x2 - xd) / 2) ||
      (lastBisection && Math.abs(x1 - x2) < delta) ||
      (!lastBisection && Math.abs(x2 - xd) < delta)
    ) {
      xs = (x0 + x1) / 2;
      lastBisection = true;
    } else {
      lastBisection = false;
    }

    ys = eq(xs);
    if (Math.abs(ys) <= ytol) {
      return xs;
    }

    [xd, x2, y2] = [x2, x1, y1];
    if ((y0 < 0 && 0 < ys) || (y0 > 0 && 0 > ys)) {
      [x1, y1] = [xs, ys];
    } else {
      [x0, y0] = [xs, ys];
    }
    if (Math.abs(y0) < Math.abs(y1)) {
      [x0, x1, y0, y1] = [x1, x0, y1, y0];
    }
  }

  if (numTries >= maxIter) {
    return Number.NaN;
  }

  return Math.abs(y1) < Math.abs(ys) ? x1 : xs;
}

function nsolve(eq, x0, x1, solver) {
  const maxIter = 20;
  const xtol = 1e-9;
  const ytol = 1e-9;

  let y0 = eq(x0);
  let y1 = eq(x1);
  if (Math.abs(y0) <= ytol) {
    return x0;
  }
  if (Math.abs(y1) <= ytol) {
    return x1;
  }
  if (Math.sign(y0) === Math.sign(y1)) {
    return Number.NaN;
  }

  return solver(eq, x0, y0, x1, y1, xtol, ytol, maxIter);
}

function calc_v0(lr, lm, vout) {
  let v0 = (vout / lm) * (lr + lm);
  return v0;
}

function calc_t1(i0, lr, lm, cr, vout) {
  let w = (lr * cr) ** -0.5;
  let z = (lr / cr) ** 0.5;

  let v0 = calc_v0(lr, lm, vout);
  let r = Math.hypot(v0 - vout, i0 * z);
  let phi = Math.atan2(v0 - vout, i0 * z);

  let tb = mod(2 * Math.PI - phi, 2 * Math.PI) / w;
  let ta = tb / 2;

  while ((r * Math.cos(w * ta + phi)) / z > i0 - (vout * ta) / lm && ta > 1e-9) {
    ta /= 2;
  }
  let t1 = ta;
  if (ta > 1e-9) {
    t1 = nsolve((t) => i0 - (vout * t) / lm - (r * Math.cos(w * t + phi)) / z, ta, tb, brent);
  }

  let i1 = (r * Math.cos(w * t1 + phi)) / z;
  let v1 = r * Math.sin(w * t1 + phi) + vout;

  return [t1, i1, v1];
}

function steady_state_i0(i0, lr, lm, cr, vbus, vout, t12) {
  // "01": t0 == 0 <= t <= t1: low-side device on, diode on
  let v0 = calc_v0(lr, lm, vout);
  let t01, i1, v1;
  [t01, i1, v1] = calc_t1(i0, lr, lm, cr, vout);

  // when diode is off
  const w1 = ((lr + lm) * cr) ** -0.5;
  const z1 = ((lr + lm) / cr) ** 0.5;

  // "12": t1 <= t <= t2: low-side device on, diode off
  let r12 = Math.hypot(v1, i1 * z1);
  let phi12 = Math.atan2(v1, i1 * z1);
  let i2 = (r12 * Math.cos(w1 * t12 + phi12)) / z1;
  let v2 = r12 * Math.sin(w1 * t12 + phi12);

  // "23": t2 <= t <= t3, high-side device on, diode off
  let r23 = Math.hypot(v2 - vbus, i2 * z1);

  // "30": t3 <= t <= Tsw, low-side device on, diode off
  let r30 = Math.hypot(v0, i0 * z1);
  // find the intersection of two ellipses that are similar
  let v3 = (r30 ** 2 - r23 ** 2 + vbus ** 2) / (2 * vbus);
  let i3;
  if (v3 < r30) {
    i3 = (r30 ** 2 - v3 ** 2) ** 0.5 / z1;
  } else {
    i3 = 0;
  }
  let t30 = (Math.atan2(v0, i0 * z1) - Math.atan2(v3, i3 * z1)) / w1;

  // when diode is on
  const w0 = (lr * cr) ** -0.5;
  const z0 = (lr / cr) ** 0.5;
  let r01 = Math.hypot(v0 - vout, i0 * z0);
  let phi01 = Math.atan2(v0 - vout, i0 * z0);

  let phi23 = Math.atan2(v2 - vbus, i2 * z1);
  let t23 = (Math.atan2(v3 - vbus, i3 * z1) - phi23) / w1;

  let phi30 = Math.atan2(v3, i3 * z1);

  return [
    { dt: t01, i0: i0, v0: v0, r: r01, phi: phi01, w: w0, z: z0, vg: vout, km: -vout / lm, hb: 0 },
    { dt: t12, i0: i1, v0: v1, r: r12, phi: phi12, w: w1, z: z1, vg: 0, km: 0, hb: 0 },
    { dt: t23, i0: i2, v0: v2, r: r23, phi: phi23, w: w1, z: z1, vg: vbus, km: 0, hb: vbus },
    { dt: t30, i0: i3, v0: v3, r: r30, phi: phi30, w: w1, z: z1, vg: 0, km: 0, hb: 0 },
  ];
}

function calc_pout_ss(ss, vout) {
  let qout = 0;
  let tsw = 0;
  for (let state of ss) {
    tsw += state.dt;
    if (state.km !== 0) {
      let q1 = ((state.i0 + state.i0 + state.km * state.dt) * state.dt) / 2;
      let q2 = (state.r / state.z / state.w) * (Math.sin(state.phi + state.w * state.dt) - Math.sin(state.phi));
      qout += q1 - q2;
    }
  }
  let iout = qout / tsw;
  let pout = iout * vout;
  return pout;
}

function calc_pout_i0(i0, lr, lm, cr, vbus, vout, t12) {
  let ss = steady_state_i0(i0, lr, lm, cr, vbus, vout, t12);
  return calc_pout_ss(ss, vout);
}

function steady_state_pout(pout, lr, lm, cr, vbus, vout, t12) {
  const eq_max_i0 = (i) => steady_state_i0(i, lr, lm, cr, vbus, vout, t12).slice(-1)[0].dt;
  let i0_max = vbus / (lr / cr) ** 0.5;
  if (((lr + lm) * cr) ** -0.5 * t12 >= Math.PI) {
    i0_max *= 100;
  } else {
    while (eq_max_i0(i0_max) > 0) {
      i0_max *= 2;
    }
    i0_max = nsolve(eq_max_i0, 1e-6, i0_max, brent);
  }

  const eq = (i) => calc_pout_i0(i, lr, lm, cr, vbus, vout, t12) - pout;
  let i0 = nsolve(eq, 1e-6, i0_max, brent);
  if (Number.isNaN(i0)) {
    let pmax = calc_pout_i0(i0_max, lr, lm, cr, vbus, vout, t12);
    throw { name: "ValueError", message: "Failed to find a steady-state.", value: pmax };
  }

  return steady_state_i0(i0, lr, lm, cr, vbus, vout, t12);
}

function multiple_diode_conduction(ss) {
  let r12 = ss[1].r;
  let v0 = ss[0].v0;
  if (r12 > v0) {
    let a = ss[1].phi;
    let iz = (r12 ** 2 - v0 ** 2) ** 0.5;
    let b = Math.atan2(v0, iz);
    let w1 = ss[1].w;
    let t = mod(b - a, 2 * Math.PI) / w1;
    let t12 = ss[1].dt;
    if (t12 >= t) {
      return true;
    }
  }
  return false;
}

/*
(() => {
  const fmt = (x) => {
    if (x === 0 || x === -0.0 || x === 0.0) {
      return "0";
    }
    let m = Math.floor(Math.log10(Math.abs(x)) / 3);
    let b = x / 10 ** (3 * m);
    let metricPrefix;
    switch (m) {
      case 8:
        metricPrefix = "Y";
        break;
      case 7:
        metricPrefix = "Z";
        break;
      case 6:
        metricPrefix = "E";
        break;
      case 5:
        metricPrefix = "P";
        break;
      case 4:
        metricPrefix = "T";
        break;
      case 3:
        metricPrefix = "G";
        break;
      case 2:
        metricPrefix = "M";
        break;
      case 1:
        metricPrefix = "k";
        break;
      case 0:
        metricPrefix = "";
        break;
      case -1:
        metricPrefix = "m";
        break;
      case -2:
        metricPrefix = "μ";
        break;
      case -3:
        metricPrefix = "n";
        break;
      case -4:
        metricPrefix = "p";
        break;
      case -5:
        metricPrefix = "f";
        break;
      case -6:
        metricPrefix = "a";
        break;
      case -7:
        metricPrefix = "z";
        break;
      case -8:
        metricPrefix = "y";
        break;
      default:
        metricPrefix = `e${(3 * m).toLocaleString(undefined, { maximumFractionDigits: 0 })}`;
    }
    let s = b.toLocaleString(undefined, { maximumFractionDigits: 3 }) + metricPrefix;
    return s;
  };
  let lr = 30e-6;
  let lm = 1e-3;
  let cr = 68e-9;
  let t12 = 800e-9;
  const test1 = () => {
    let vbus = 410;
    let vout = 150;
    for (let i0 of [0.001, 0.01, 0.1, 0.2, 0.5, 1, 2, 5]) {
      let ss = steady_state_i0(i0, lr, lm, cr, vbus, vout, t12);
      console.log(
        `i0=${i0}`,
        ss.reduce((acc, el) => acc + " " + fmt(el.dt), "")
      );
    }
    let ss = steady_state_pout(2, lr, lm, cr, vbus, vout, t12);
    console.log(calc_pout_ss(ss, vout));
  };
  const test2 = () => {
    let vbus = 410;
    // let vout = 77 * 5;
    // let p = 40;
    let t12 = 1500e-9;
    let vout = 80 * 3.5;
    let p = 80;
    try {
      let ss = steady_state_pout(p, lr, lm, cr, vbus, vout, t12);
      console.log(calc_pout_ss(ss, vout));
    } catch (e) {
      console.error(e.name);
      console.error(e.message);
      console.error(`Maximum possible output power: ${e.value}.`);
    }
  };
  // test1();
  test2();
})();
*/
