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
