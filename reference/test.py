import ahbllc
import ahbllcapp
import plothelper
from plothelper import (fmt, axis_formatter, plot)

import math
import numpy as np
from scipy.optimize import brentq as nsolve
import matplotlib.pyplot as plt


def test1():
    ckt = ahbllc.AHBLLCTrafo(lr=25e-6, lm=1225e-6, cr=39e-9, nps=4.3, vbus=410, vload=40, chb=500e-12)
    print(ckt)
    voff, ss, ev = ahbllc.evaluate_operating_point(1, ckt)
    ev.set_circuit(ckt)
    print(ev)
    print(ss[0].i0, '-->', ahbllc.evaluate_switching_period(ss).iout * ckt.vout, fmt(1/ahbllc.evaluate_switching_period(ss).tsw))
    plothelper.plot(ss, show=True)
    # for dv in [.001, .01, .1, .2, .5, 1, 2, 5, 10]:
    #     ss = ahbllc.find_steady_state(voff+dv, ckt)
    #     print(ss[0].i0, '-->', ahbllc.evaluate_switching_period(ss).iout * ckt.vout, fmt(1/ahbllc.evaluate_switching_period(ss).tsw))
    #     plothelper.plot(ss, show=True)


def test2():
    ckt = ahbllc.AHBLLCTrafo(lr=25e-6, lm=1225e-6, cr=39e-9, nps=4.3, vbus=410, vload=80, chb=500e-12)
    ahbllcapp.evaluate_operating_point_with_tolerance(20, ckt)


def test3():
    ckt = ahbllc.AHBLLCTrafo(lr=25e-6, lm=1225e-6, cr=39e-9, nps=4.3, vbus=410, vload=40, chb=500e-12)
    dvoff = 10
    # dv, ss = ahbllc.sim(330, ckt, (dvoff, 500e-9))
    # print(f'v0 = {fmt(330)}V ==> dv = {fmt(dv, 5)}V')
    # plothelper.plot(ss, show=True)
    for i, v0 in enumerate(range(-100, 400, 10)):
        dv, ss = ahbllc.sim(v0, ckt, (dvoff, 500e-9))
        print(f'v0 = {fmt(v0)}V ==> dv = {fmt(dv, 5)}V')
        plothelper.plot(ss, show=False, filename=f'{i:03d}')


if __name__ == '__main__':
    test3()
