import ahbllc
import ahbllcapp
import plothelper
from plothelper import (fmt, axis_formatter, plot)

import math
import numpy as np
from scipy.optimize import brentq as nsolve
import matplotlib.pyplot as plt


def test1():
    ckt = ahbllc.AHBLLCTrafo(lr=25e-6, lm=1225e-6, cr=39e-9, nps=4.3, vbus=410, vload=80, chb=500e-12)
    print(ckt)
    voff, ss, ev = ahbllc.evaluate_operating_point(1, ckt)
    ev.set_circuit(ckt)
    print(ev)
    plothelper.plot(ss, show=True)
    ss = ahbllc.find_steady_state(voff+1, ckt)
    plothelper.plot(ss, show=True)


def test2():
    ckt = ahbllc.AHBLLCTrafo(lr=25e-6, lm=1225e-6, cr=39e-9, nps=4.3, vbus=410, vload=80, chb=500e-12)
    ahbllcapp.evaluate_operating_point_with_tolerance(20, ckt)


if __name__ == '__main__':
    test1()
