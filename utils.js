function fmt(x) {
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
}


function mod(a, b) {
  let c = a % b;
  if ((c < 0 && 0 < b) || (c > 0 && 0 > b)) {
    c += b;
  }
  return c;
}

function parg(angle) {
  return mod(angle, 2 * Math.PI);
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

