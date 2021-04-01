const CONTROL = {
  MINIMUM_TIME: 500e-9,
  MINIMUM_VOLTAGE: 1e-6,
  MINIMUM_CURRENT: 1e-6,
  MAXIMUM_COMMUTATION_TIME: 600e-9,
  MINIMUM_HIGH_SIDE_ON_TIME: 100e-9,
  MINIMUM_FORWARD_TIME: 500e-9,
  MAXIMUM_T12: 20e-6,
  MINIMUM_DETECTABLE_DIODE_ON_TIME: 500e-9,
};

function stateL1(i0, v0, vhb0, im0, ckt, con) {
  const chb = Infinity;
  const cr = ckt.cr;
  const ctot = 1 / (1 / cr + 1 / chb);
  const ltot = ckt.lr;
  const w = (ltot * ctot) ** -0.5;
  const z = (ltot / ctot) ** 0.5;
  const km = ckt.vout / ckt.lm;

  vhb0 = 0;
  const vout = ckt.vout;
  const vcen = vhb0 + vout;
  const r = Math.hypot(v0 - vcen, i0 * z);
  const phi = Math.atan2(v0 - vcen, i0 * z);

  let tb = parg(-phi) / w;
  let ta = tb;
  let dt;
  while ((r * Math.cos(w * ta + phi)) / z > im0 - km * ta && ta > CONTROL.MINIMUM_TIME) {
    ta /= 2;
  }
  if (ta > CONTROL.MINIMUM_TIME) {
    dt = nsolve((t) => im0 - km * t - (r * Math.cos(w * t + phi)) / z, ta, tb, brent);
  } else {
    dt = ta;
  }

  const nextState = stateL0;
  const i1 = (r * Math.cos(w * dt + phi)) / z;
  const v1 = r * Math.sin(w * dt + phi) + vcen;
  const vhb1 = vhb0;
  const im1 = i1;
  return [
    nextState,
    {
      state: "l1",
      dt: dt,
      i0: i0,
      v0: v0,
      vhb0: vhb0,
      im0: im0,
      i1: i1,
      v1: v1,
      vhb1: vhb1,
      im1: im1,
      r: r,
      phi: phi,
      w: w,
      z: z,
      cr: cr,
      chb: chb,
      vcen: vcen,
      vout: vout,
      km: -km,
    },
  ];
}

function stateL0(i0, v0, vhb0, im0, ckt, con) {
  const t12 = con;

  const chb = Infinity;
  const cr = ckt.cr;
  const ctot = 1 / (1 / cr + 1 / chb);
  const ltot = ckt.lr + ckt.lm;
  const w = (ltot * ctot) ** -0.5;
  const z = (ltot / ctot) ** 0.5;
  const km = 0;

  vhb0 = 0;
  const vout = 0;
  const vcen = vhb0 + vout;
  const r = Math.hypot(v0 - vcen, i0 * z);
  const phi = Math.atan2(v0 - vcen, i0 * z);

  const vdon = (ckt.vout / ckt.lm) * (ckt.lr + ckt.lm);
  let dtdon = Infinity;
  if (-r <= vdon - vcen && vdon - vcen <= r) {
    let a = Math.asin((vdon - vcen) / r);
    dtdon = Math.min(parg(a - phi), parg(Math.PI - a - phi)) / w;
  }

  const it12 = (r * Math.cos(w * t12 + phi)) / z;
  const dtizc = parg(Math.PI / 2 - phi) / w;

  let dt;
  let nextState;
  if (dtdon < t12 || (it12 > 0 && dtdon < dtizc)) {
    dt = dtdon;
    nextState = stateL1;
  } else if (it12 <= 0) {
    dt = t12;
    nextState = stateC0;
  } else {
    dt = dtizc;
    nextState = stateC0;
  }

  const i1 = (r * Math.cos(w * dt + phi)) / z;
  const v1 = r * Math.sin(w * dt + phi) + vcen;
  const vhb1 = vhb0;
  const im1 = i1;
  return [
    nextState,
    {
      state: "l0",
      dt: dt,
      i0: i0,
      v0: v0,
      vhb0: vhb0,
      im0: im0,
      i1: i1,
      v1: v1,
      vhb1: vhb1,
      im1: im1,
      r: r,
      phi: phi,
      w: w,
      z: z,
      cr: cr,
      chb: chb,
      vcen: vcen,
      vout: vout,
      km: -km,
    },
  ];
}

function stateH0(i0, v0, vhb0, im0, ckt, con) {
  let dt = con;

  const chb = Infinity;
  const cr = ckt.cr;
  const ctot = 1 / (1 / cr + 1 / chb);
  const ltot = ckt.lr + ckt.lm;
  const w = (ltot * ctot) ** -0.5;
  const z = (ltot / ctot) ** 0.5;
  const km = 0;

  vhb0 = ckt.vbus;
  const vout = 0;
  const vcen = vhb0 + vout;
  const r = Math.hypot(v0 - vcen, i0 * z);
  const phi = Math.atan2(v0 - vcen, i0 * z);

  if (isNaN(dt)) {
    dt = parg(-Math.PI / 2 - phi) / w;
  }

  nextState = stateC0;
  const i1 = (r * Math.cos(w * dt + phi)) / z;
  const v1 = r * Math.sin(w * dt + phi) + vcen;
  const vhb1 = vhb0;
  const im1 = i1;
  return [
    nextState,
    {
      state: "h0",
      dt: dt,
      i0: i0,
      v0: v0,
      vhb0: vhb0,
      im0: im0,
      i1: i1,
      v1: v1,
      vhb1: vhb1,
      im1: im1,
      r: r,
      phi: phi,
      w: w,
      z: z,
      cr: cr,
      chb: chb,
      vcen: vcen,
      vout: vout,
      km: -km,
    },
  ];
}

function stateC0(i0, v0, vhb0, im0, ckt, con) {
  const chb = ckt.chb;
  const cr = ckt.cr;
  const ctot = 1 / (1 / cr + 1 / chb);
  const ltot = ckt.lr + ckt.lm;
  const w = (ltot * ctot) ** -0.5;
  const z = (ltot / ctot) ** 0.5;
  const km = 0;

  const vout = 0;
  const vcen = vout;
  const vcap0 = v0 - vhb0;
  const r = Math.hypot(vcap0 - vcen, i0 * z);
  const phi = Math.atan2(vcap0 - vcen, i0 * z);

  let dt;
  let nextState;
  if (Math.abs(vhb0) < CONTROL.MINIMUM_VOLTAGE) {
    let v = ckt.vbus;
    v -= -(vout - v0) / (1 + chb / cr) + vhb0 / (1 + cr / chb);
    v *= 1 + chb / cr;
    if (-r <= v && v <= r) {
      let a = Math.asin(v / -r);
      dt = Math.min(parg(a - phi), parg(Math.PI - a - phi)) / w;
    } else {
      dt = ((-Math.PI / 2 - phi) % (2 * Math.PI)) / w;
    }
    nextState = stateH0;
  } else {
    let dtcom;
    let v = 0;
    v -= -(vout - v0) / (1 + chb / cr) + vhb0 / (1 + cr / chb);
    v *= 1 + chb / cr;
    if (-r <= v && v <= r) {
      let a = Math.asin(v / -r);
      dtcom = Math.min(parg(a - phi), parg(Math.PI - a - phi)) / w;
    } else {
      dtcom = parg(Math.PI / 2 - phi) / w;
    }

    let dtdon;
    v = (ckt.vout / ckt.lm) * (ckt.lr + ckt.lm);
    v -= vcen;
    if (-r <= v && v <= r) {
      let a = Math.asin(v / r);
      dtdon = Math.min(parg(a - phi), parg(Math.PI - a - phi)) / w;
    } else {
      dtdon = Infinity;
    }

    if (dtdon < dtcom) {
      dt = dtdon;
      nextState = stateC1;
    } else {
      dt = dtcom;
      nextState = stateL0;
    }
  }

  const i1 = (r * Math.cos(w * dt + phi)) / z;
  const v1 = (r * Math.sin(w * dt + phi)) / (1 + cr / chb) + v0 / (1 + chb / cr) + (vhb0 + vout) / (1 + cr / chb);
  let vhb1 = (-r * Math.sin(w * dt + phi)) / (1 + chb / cr) - (vout - v0) / (1 + chb / cr) + vhb0 / (1 + cr / chb);
  if (Math.abs(vhb1 - ckt.vbus) < CONTROL.MINIMUM_VOLTAGE) {
    vhb1 = ckt.vbus;
  } else if (Math.abs(vhb1) < CONTROL.MINIMUM_VOLTAGE) {
    vhb1 = 0;
  }
  const im1 = im0;
  return [
    nextState,
    {
      state: "c0",
      dt: dt,
      i0: i0,
      v0: v0,
      vhb0: vhb0,
      im0: im0,
      i1: i1,
      v1: v1,
      vhb1: vhb1,
      im1: im1,
      r: r,
      phi: phi,
      w: w,
      z: z,
      cr: cr,
      chb: chb,
      vcen: vcen,
      vout: vout,
      km: -km,
    },
  ];
}

function stateC1(i0, v0, vhb0, im0, ckt, con) {
  const chb = ckt.chb;
  const cr = ckt.cr;
  const ctot = 1 / (1 / cr + 1 / chb);
  const ltot = ckt.lr;
  const w = (ltot * ctot) ** -0.5;
  const z = (ltot / ctot) ** 0.5;
  const km = ckt.vout / ckt.lm;

  const vout = 0;
  const vcen = vout;
  const vcap0 = v0 - vhb0;
  const r = Math.hypot(vcap0 - vcen, i0 * z);
  const phi = Math.atan2(vcap0 - vcen, i0 * z);

  let v = 0;
  v -= -(vout - v0) / (1 + chb / cr) + vhb0 / (1 + cr / chb);
  v *= 1 + chb / cr;
  let dtcom;
  if (-r <= v && v <= r) {
    let a = Math.asin(v / -r);
    dtcom = Math.min(parg(a - phi), parg(Math.PI - a - phi)) / w;
  } else {
    dtcom = parg(Math.PI / 2 - phi) / w;
  }

  let tb = parg(-phi) / w;
  let ta = tb;
  while ((r * Math.cos(w * ta + phi)) / z > im0 - km * ta && ta > CONTROL.MINIMUM_TIME) {
    ta /= 2;
  }
  let dtdof;
  if (ta > CONTROL.MINIMUM_TIME) {
    dtdof = nsolve((t) => im0 - km * t - (r * Math.cos(w * t + phi)) / z, ta, tb, brent);
  } else {
    dtdof = ta;
  }
  let nextState;
  if (dtdof < dtcom) {
    dt = dtdof;
    nextState = stateC0;
  } else {
    dt = dtcom;
    nextState = stateL1;
  }

  const i1 = (r * Math.cos(w * dt + phi)) / z;
  const v1 = (r * Math.sin(w * dt + phi)) / (1 + cr / chb) + v0 / (1 + chb / cr) + (vhb0 + vout) / (1 + cr / chb);
  let vhb1 = (-r * Math.sin(w * dt + phi)) / (1 + chb / cr) - (vout - v0) / (1 + chb / cr) + vhb0 / (1 + cr / chb);
  if (vhb1 - ckt.vbus < CONTROL.MINIMUM_VOLTAGE) {
    vhb1 = ckt.vbus;
  } else if (Math.abs(vhb1) < CONTROL.MINIMUM_VOLTAGE) {
    vhb1 = 0;
  }
  let im1 = im0 - km * dt;
  if (Math.abs(im1 - i1) < 1e-6) {
    im1 = i1;
  }
  return [
    nextState,
    {
      state: "c1",
      dt: dt,
      i0: i0,
      v0: v0,
      vhb0: vhb0,
      im0: im0,
      i1: i1,
      v1: v1,
      vhb1: vhb1,
      im1: im1,
      r: r,
      phi: phi,
      w: w,
      z: z,
      cr: cr,
      chb: chb,
      vcen: vcen,
      vout: vout,
      km: -km,
    },
  ];
}

function simulatePhase(isf, cond, i0, v0, vhb0, im0, ckt, con) {
  const isLowSideOnPhase = (isf === stateL0 || isf === stateL1) && !isNaN(con);

  if (isLowSideOnPhase) {
    const t12given = con;
    let t12history = [];

    let nsf;
    let s;
    [nsf, s] = isf(i0, v0, vhb0, im0, ckt, con);

    let states = [s];
    while (cond(nsf)) {
      let csf = nsf;
      [nsf, s] = csf(s.i1, s.v1, s.vhb1, s.im1, ckt, con);

      if (csf === stateL0) {
        t12history.push(s.dt);
      } else {
        if (s.dt < CONTROL.MINIMUM_DETECTABLE_DIODE_ON_TIME) {
          con = t12given - t12history.reduce((x, y) => x + y)
        } else {
          t12history = [];
          con = t12given;
        }
      }

      states.push(s)
    }

    return [nsf, states]
  } else {
    let nsf;
    let s;
    [nsf, s] = isf(i0, v0, vhb0, im0, ckt, con);
    let states = [s];
    while (cond(nsf)) {
      [nsf, s] = nsf(s.i1, s.v1, s.vhb1, s.im1, ckt, con);
      states.push(s)
    }
    return [nsf, states]
  }
}

function simulate(v0, ckt, con) {
  const isCommutating = (s) => s === stateC0 || s === stateC1;
  const isActive = (s) => !isCommutating(s);

  const maximumDv = (i0, v0) => {
    const l = ckt.lr + ckt.lm;
    const c = 1 / (1 / ckt.cr + 1 / ckt.chb);
    const z = (l / c) ** 0.5;
    const r = Math.hypot(v0, i0 * z);
    const phi = Math.atan2(v0, i0 * z);
    return r * (1 + Math.sin(phi));
  };

  const equationZVon = (t) => {
    let hsOffLast = hsOff.slice(-1)[0];
    let lsOn = simulatePhase(

      lsOnIsf,

           isActive,

      hsOffLast.i1,

      hsOffLast.v1,

          hsOffLast.vhb1,

         hsOffLast.im1,

        ckt,

       t

    ).slice(-1)[0];
    let lsOnLast = lsOn.slice(-1)[0];
    return maximumDv(lsOnLast.i1, lsOnLast.v1) - (1 + ckt.chb / ckt.cr) * ckt.vbus;
  };

  let conv;
  let t12min;
  [conv, t12min] = con;
  const i0 = 0;
  const vhb0 = ckt.vbus;

  let nsf;
  let hsOn0;
  [nsf, hsOn0] = simulatePhase(stateH0, isActive, i0, v0, vhb0, i0, ckt, conv);
  let hsOn0Last = hsOn0.slice(-1)[0];

  let lsOnIsf;
  let hsOff;
  [lsOnIsf, hsOff] = simulatePhase(
    nsf,
    isCommutating,
    hsOn0Last.i1,
    hsOn0Last.v1,
    hsOn0Last.vhb1,
    hsOn0Last.im1,
    ckt,
    NaN
  );
  let hsOffLast = hsOff.slice(-1)[0];

  let lsOn;
  lsOn = simulatePhase(lsOnIsf, isActive, hsOffLast.i1, hsOffLast.v1, hsOffLast.vhb1, hsOffLast.im1, ckt, t12min).slice(

    -1

  )[0];
  let lsOnLast = lsOn.slice(-1)[0];
  if (maximumDv(lsOnLast.i1, lsOnLast.v1) < (1 + ckt.chb / ckt.cr) * ckt.vbus) {
    const t12max = (Math.PI - lsOnLast.phi) / lsOnLast.w;
    let t12ZVon;
    if (t12min < t12max && t12max <= CONTROL.MAXIMUM_T12) {
      if (0 < equationZVon(t12max)) {
        t12ZVon = nsolve(equationZVon, t12min, t12max, brent);
      } else {
        t12ZVon = t12max;
      }
    }
    lsOn = simulatePhase(
      lsOnIsf,
      isActive,
      hsOffLast.i1,
      hsOffLast.v1,
      hsOffLast.vhb1,
      hsOffLast.im1,
      ckt,
      t12ZVon
    ).slice(-1)[0];
  }
  lsOnLast = lsOn.slice(-1)[0];

  let lsOff = simulatePhase(
    stateC0,
    isCommutating,
    lsOnLast.i1,
    lsOnLast.v1,
    lsOnLast.vhb1,
    lsOnLast.im1,
    ckt,
    NaN
  ).slice(-1)[0];
  let lsOffLast = lsOff.slice(-1)[0];

  let hsOn1 = simulatePhase(
    stateH0,
    isActive,
    lsOffLast.i1,
    lsOffLast.v1,
    lsOffLast.vhb1,
    lsOffLast.im1,
    ckt,
    NaN
  ).slice(-1)[0];

  return [hsOn1.slice(-1)[0].v1 - v0, hsOn0.concat(hsOff).concat(lsOn).concat(lsOff).concat(hsOn1)];
}

function evaluateSwitchingPeriod(states) {
  let operation;
  if (
    Math.abs(states[0].i0 - states.slice(-1)[0].i1) < CONTROL.MINIMUM_CURRENT &&
    Math.abs(states[0].v0 - states.slice(-1)[0].v1) < CONTROL.MINIMUM_VOLTAGE
  ) {
    operation = "steady-state";
  } else {
    operation = "dynamics";
  }

  let tsw = 0;
  let qout = 0;
  let i2acc = 0;
  let i2accOut = 0;
  let imax = 0;
  let imaxOut = 0;
  let vacc = 0;
  let v2acc = 0;

  let tDiodeOn = 0;
  let tHsOff = 0;
  let tLsOn = 0;
  let tLsOff = 0;
  let tHsOn = 0;

  const sin = Math.sin;
  const cos = Math.cos;

  for (state of states) {
    tsw += state.dt;
    j = state.r / state.z;
    i2acc += (j ** 2 * state.dt) / 2;
    i2acc += (j ** 2 / (4 * state.w)) * (sin(2 * state.w * state.dt + 2 * state.phi) - sin(2 * state.phi));

    m = state.r / (1 + state.cr / state.chb);
    b = state.v0 / (1 + state.chb / state.cr);
    b += (state.vhb0 + state.vout) / (1 + state.cr / state.chb);
    vacc += (m / state.w) * (cos(state.phi) - cos(state.w * state.dt + state.phi));
    vacc += b * state.dt;
    v2acc += (m ** 2 / 2 + b ** 2) * state.dt;
    v2acc -= (m ** 2 / (4 * state.w)) * (sin(2 * state.w * state.dt + 2 * state.phi) - sin(2 * state.phi));
    v2acc += ((2 * m * b) / state.w) * (cos(state.phi) - cos(state.w * state.dt + state.phi));

    if (state.state.slice(-1)[0] === "1") {
    }
    if (state.state[0] === "c") {
    } else if (state.state[0] === "l") {
    } else if (state.state[0] === "h") {
    }
    if (state.phi <= 0 && 0 <= state.w * state.dt + state.phi) {
    }
  }

  const irms = (i2acc / tsw) ** 0.5;
  const irmsOut = (i2accOut / tsw) ** 0.5;
  const vavg = vacc / tsw;
  const vrms = (v2acc / tsw) ** 0.5;
  const iout = qout / tsw;

  t12 = states.filter((s) => s.state === "l0").slice(-1)[0].dt;

  return {
    operation: operation,
    tsw: tsw,
    irms: irms,
    irmsOut: irmsOut,
    imax: imax,
    imaxOut: imaxOut,
    vavg: vavg,
    vrms: vrms,
    iout: iout,
    tDiodeOn: tDiodeOn,
    tHsOff: tHsOff,
    tLsOn: tLsOn,
    tLsOff: tLsOff,
    tHsOn: tHsOn,
    t12: t12,
  };
}
