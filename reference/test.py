import ahbllc
import ahbllcapp
import plothelper
from plothelper import (fmt, axis_formatter, plot)

import math
import numpy as np
from scipy.optimize import brentq as nsolve
import matplotlib.pyplot as plt


def test():
    ckt = ahbllc.AHBLLCTrafo(25e-6, 1225e-6, 39e-9, 4.3, 410, 40, chb=500e-12)
    print(ckt)
    ss, ev = ahbllc.evaluate_operating_point(.4, ckt)
    ev.set_circuit(ckt)
    print(ev)
    plothelper.plot(ss, show=True)


if __name__ == '__main__':
    test()
