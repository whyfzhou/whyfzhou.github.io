import ahbllc
import ahbllcapp
import plothelper
from plothelper import (fmt, axis_formatter, plot)

import math
import numpy as np
from scipy.optimize import brentq as nsolve
import matplotlib.pyplot as plt


def test():
    ckt = ahbllc.AsymmetricalHalfBridgeLLC(25e-6, 1225e-6, 39e-9, 410, 40 * 4.3, chb=500e-12)
    fswmax = 150e3
    i0max = ckt.vbus / ((ckt.lr + ckt.lm) / ckt.cr)**.5

    t12min = .5e-6
    voff = 166
    # voff = 116  # voff 有可能会太小，以至于如果 t12 不增加的话，周期终值电流大于任何初始假设（包括 1mA），会发散，所以存在一个最低 Voff，需要提前找到
    def eq(i):
        di, _, ss = ahbllc.sim(i, ckt, (voff, t12min))
        fsw = 1 / sum(s['dt'] for s in ss)
        if fsw > fswmax:
            def eqf(t):
                _, _, ss = ahbllc.sim(i, ckt, (voff, t))
                fsw = 1 / sum(s['dt'] for s in ss)
                return fsw - fswmax
            t = nsolve(eqf, .5e-6, 1 / fswmax)
            di, _, ss = ahbllc.sim(i, ckt, (voff, t))
        else:
            t = t12min
        dv = ss[-1]['v1'] - voff
        return di, dv, t, ss

    # for v in [170, 160, 150, 140, 130, 120, 110]:
    #     _, _, ss = ahbllc.sim(.4, ckt, (v, t12))
    #     print(f"resonant voltage intial value : {fmt(v)}V, final value: {fmt(ss[-1]['v1'])}V")

    voff = 158
    # for i in np.linspace(1e-3, i0max, 10):
    #     di, dv, ss = eq(i)
    #     print(f'i0={fmt(i)}A, dv={fmt(dv)}V, di={fmt(di)}A')
    #     plothelper.plot(ss, show=True)
    i0 = nsolve(lambda i: eq(i)[0], 1e-3, i0max)
    _, _, t12min, ss = eq(i0)
    plothelper.plot(ss, show=True)

    # _, _, ss = ahbllc.sim(0.466, ckt, (voff, t12))
    # plothelper.plot(ss, show=True)

    # _, _, ss = ahbllc.sim(0.159, ckt, (voff, t12))
    # plothelper.plot(ss, show=True)

    # for i0it in [.001, .01, .1] + list(range(10 + 1 + 1)):
    #     i0 = .159 + i0it / 10 * (.466 - .159)
    #     _, _, ss = ahbllc.sim(i0, ckt, (voff, t12))
    #     plothelper.plot(ss, show=True)

    # i0 = nsolve(eq, .001, i0max)  # TODO: initial bracket
    # # i0 = nsolve(eq, i0max/20, i0max)  # TODO: initial bracket
    # print(f'solved i0 = {fmt(i0)}A')
    # _, t12, ss = ahbllc.sim(i0, ckt, (voff, t12))
    # print(f'corrected t12 = {fmt(t12)}s')
    # plothelper.plot(ss, show=True)

    # sampled_data = ahbllc.sample(ss)
    # perf = ahbllc.evaluate_sampled(*sampled_data, ss)
    # plothelper.plot_sampled(*sampled_data, ss, perf, t12, ckt, show=True)


def test1():
    t12 = .5e-6
    # voff = 166
    voff = 116  # voff 有可能会太小，以至于如果 t12 不增加的话，周期终值电流大于任何初始假设（包括 1mA），会发散，所以存在一个最低 Voff，需要提前找到
    ckt = ahbllc.AsymmetricalHalfBridgeLLC(25e-6, 1225e-6, 39e-9, 410, 40 * 4.3, chb=500e-12)

    fmax = 100e3
    def eq(i):
        def eqf(t):
            _, _, ss = ahbllc.sim(i, ckt, (voff, t))
            fsw = 1 / sum(s['dt'] for s in ss)
            return fsw - fmax
        _, _, ss = ahbllc.sim(i, ckt, (voff, t12))
        if 1 / sum(s['dt'] for s in ss) > fmax:
            t = nsolve(eqf, t12, 1 / fmax)
        else:
            t = t12
        return ahbllc.sim(i, ckt, (voff, t))[0]

    eq(.001)
    eq(.01)
    eq(.1)
    eq(.5)
    i0max = ckt.vbus / ((ckt.lr + ckt.lm) / ckt.cr)**.5
    i0 = nsolve(eq, 1e-3, i0max)
    _, _, ss = ahbllc.sim(i0, ckt, (voff, t12))
    plothelper.plot(ss, show=True)


def test2():
    ckt = ahbllc.AsymmetricalHalfBridgeLLC(25e-6, 1225e-6, 39e-9, 410, 40 * 4.3, chb=500e-12)
    voffs = np.linspace(150, 250, 50)
    pouts = []
    for voff in voffs:
        t12, ss = ahbllc.find_steady_state(voff, ckt)
        print(f'corrected t12 = {fmt(t12)}s', end=', ')
        iout = ahbllc.calculate_iout(ss)
        pout = iout * ckt.vout
        pouts.append(pout)
        print(f'iout = {fmt(iout)}A, pout = {fmt(pout)}W')
        print('  ', ahbllc.evaluate_steady_state(ss))
        # fig = plothelper.plot(ss, show=True)
        # plt.close(fig)
    plt.plot(voffs, pouts)
    plt.show()


def test3():
    ckt = ahbllc.AHBLLCTrafo(25e-6, 1225e-6, 39e-9, 4.3, 410, 40, chb=500e-12)
    print(ckt)
    t12, ss, ev = ahbllc.evaluate_operating_point(40, ckt)
    print(t12)
    ev.set_circuit(ckt)
    print(ev)
    plothelper.plot(ss, show=True)


if __name__ == '__main__':
    test3()
