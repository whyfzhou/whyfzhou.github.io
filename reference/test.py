import ahbllc
import ahbllcapp
from plothelper import (fmt, plot)

import math
import numpy as np
from scipy.optimize import (brentq as nsolve, minimize)
import matplotlib.pyplot as plt


def test_case_01():
    ckt = ahbllc.AHBLLCTrafo(lr=25e-6, lm=1225e-6, cr=39e-9, nps=4.3, vbus=410, vload=20, chb=500e-12)
    con = 200e-9
    v0 = nsolve(lambda v: ahbllc.sim(v, ckt, (con, 500e-9))[0], -ckt.vbus, ckt.vout / ckt.lm * (ckt.lr + ckt.lm))
    _, ss = ahbllc.sim(v0, ckt, (con, 500e-9))
    fig, *_ = plot(ss, show=True)
    ev = ahbllc.evaluate_switching_period(ss)
    print(ev.set_circuit(ckt))
    plt.close(fig)
    # for con in np.arange(1e-6, 4e-6, 1e-6):
    #     v0 = nsolve(lambda v: ahbllc.sim(v, ckt, (con, 500e-9))[0], -ckt.vbus, ckt.vout / ckt.lm * (ckt.lr + ckt.lm))
    #     _, ss = ahbllc.sim(v0, ckt, (con, 500e-9))
    #     fig, *_ = plot(ss, show=True)
    #     plt.close(fig)


def test_case_02():
    ckt = ahbllc.AHBLLCTrafo(lr=25e-6, lm=1225e-6, cr=39e-9, nps=4.3, vbus=410, vload=20, chb=500e-12)
    ss = ahbllc.find_steady_state(2500e-9, ckt)
    ev = ahbllc.evaluate_switching_period(ss)
    print(ev.set_circuit(ckt))
    fig, *_ = plot(ss, show=True)
    plt.close(fig)


def test_case_03():
    ckt = ahbllc.AHBLLCTrafo(lr=31.5e-6, lm=945e-6, cr=49.4e-9, nps=1, vbus=390, vload=330, chb=630e-12)
    tf, ss, ev, pmax = ahbllc.evaluate_operating_point(133, ckt)
    print(f'Control variable: high-side forward time: {fmt(tf, 4)}s')
    ev = ahbllc.evaluate_switching_period(ss)
    print(ev.set_circuit(ckt))
    fig, *_ = plot(ss, show=True)
    plt.close(fig)


def test_case_04():
    ckt = ahbllc.AHBLLC(lr=57e-6, lm=855e-6, cr=44.6e-9, vbus=390, vout=330, chb=570e-12)
    ss = ahbllc.find_steady_state(500e-9, ckt)
    tf, ss, ev, pmax = ahbllc.evaluate_operating_point(pout=150, ckt=ckt)
    print(f'Control variable: high-side forward time: {fmt(tf, 4)}s')
    ev = ahbllc.evaluate_switching_period(ss)
    print(ev.set_circuit(ckt))
    fig, *_ = plot(ss, show=True)
    plt.close(fig)


def test_case_05():
    ckt = ahbllc.AHBLLCTrafo(lr=30e-6, lm=900e-6, cr=47e-9, nps=2, vbus=410, vload=150, chb=600e-12)
    ahbllcapp.evaluate_operating_point_with_tolerance(pout=150, ckt=ckt)


def test_case_06():
    ckt = ahbllc.AHBLLC(lr=30e-6, lm=1000e-6, cr=47e-9, chb=200e-12, vbus=410, vout=80)
    tf, ss, ev, pmax = ahbllc.evaluate_operating_point(pout=1, ckt=ckt)
    print(ev)


if __name__ == '__main__':
    test_case_06()

