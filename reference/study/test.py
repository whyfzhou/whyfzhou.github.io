import ahbllc
import ahbllcapp
import plothelper
from plothelper import (fmt, axis_formatter, plot)

import math
import numpy as np
from scipy.optimize import (brentq as nsolve, minimize)
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
    ckt = ahbllc.AHBLLCTrafo(lr=25e-6, lm=1225e-6, cr=39e-9, nps=4.3, vbus=410, vload=80, chb=500e-12)
    vinc = 12
    # dv, ss = ahbllc.sim(161, ckt, (dvoff, 500e-9), constr=ahbllc.state_h0_vinc)
    # print(f'v0 = {fmt(161)}V ==> dv = {fmt(dv, 5)}V')
    # plothelper.plot(ss, show=True)
    for vv in range(333490, 333500, 1):
        v0 = vv / 1000
        dv, ss = ahbllc.sim(v0, ckt, (vinc, 500e-9), constr=ahbllc.state_h0_vinc)
        print(f'v0 = {fmt(v0, 6)}V ==> dv = {fmt(dv, 5)}V')
        fig, *_ = plothelper.plot(ss, show=False, filename='{:0>10}'.format(f'{int(v0 * 1000):_}'),
                                  sphlines=[ckt.vout / ckt.lm * (ckt.lr + ckt.lm)])
        plt.close(fig)


def test4():
    # vinc 控制的不稳定性
    ckt = ahbllc.AHBLLCTrafo(lr=25e-6, lm=1225e-6, cr=39e-9, nps=4.3, vbus=410, vload=80, chb=500e-12)
    vinc = 12
    v0 = 333.496
    sss = []
    for _ in range(100):
        _, ss = ahbllc.sim(v0, ckt, (vinc, 500e-9), constr=ahbllc.state_h0_vinc)
        v0 = ss[-1].v1
        sss += ss
    fig, *_ = plothelper.plot(sss, show=True, sphlines=[ckt.vout / ckt.lm * (ckt.lr + ckt.lm)],
                            spradii=False, splegends=False)
    plt.close(fig)


def test5():
    # vinc 控制的不稳定性
    ckt = ahbllc.AHBLLCTrafo(lr=25e-6, lm=1225e-6, cr=39e-9, nps=4.3, vbus=410, vload=80, chb=500e-12)
    vinc = 49
    v0 = 287
    sss = []
    for _ in range(100):
        _, ss = ahbllc.sim(v0, ckt, (vinc, 500e-9), constr=ahbllc.state_h0_vinc)
        v0 = ss[-1].v1
        print(f'v0 = {fmt(ss[0].v0)}V, actual dvoff = {fmt(ss[0].v1 - ss[0].v0, 4)}V')
        sss += ss
    fig, *_ = plothelper.plot(sss, show=True, sphlines=[ckt.vout / ckt.lm * (ckt.lr + ckt.lm)],
                            spradii=False, splegends=False)
    plt.close(fig)


def test6():
    # 用 dvoff 控制，从随意给定的初始值开始仿真，取稳定后的几个状态做图
    ckt = ahbllc.AHBLLCTrafo(lr=25e-6, lm=1225e-6, cr=39e-9, nps=4.3, vbus=410, vload=80, chb=500e-12)
    dvoff = 1
    i0 = -.5
    v0 = 272
    sss = []
    for _ in range(500):
        _, ss = ahbllc.sim_dvoff(i0, v0, ckt, (dvoff, 500e-9))
        i0, v0 = ss[-1].i1, ss[-1].v1
        sss += ss
    fig, *_ = plothelper.plot(sss[-30:], show=True, sphlines=[ckt.vout / ckt.lm * (ckt.lr + ckt.lm)],
                            spradii=False, splegends=False)
    plt.close(fig)
    print(ss[-1].i1, ss[-1].v1)


def test7():
    # 用 dvoff 控制，尝试用 scipy.optimize.minimize 直接解起始点位置
    ckt = ahbllc.AHBLLCTrafo(lr=25e-6, lm=1225e-6, cr=39e-9, nps=4.3, vbus=410, vload=80, chb=500e-12)
    dvoff = 1
    i0 = -.499
    v0 = 336.75
    eq = lambda x: math.hypot(*ahbllc.sim_dvoff(x[0], x[1], ckt, (dvoff, 500e-9))[0])
    res = minimize(eq, (i0, v0), method='L-BFGS-B')
    print(res)
    print(eq((-0.49902207045580427, 336.75073812079916)))
    i0, v0 = res.x
    _, ss = ahbllc.sim_dvoff(i0, v0, ckt, (dvoff, 500e-9))
    fig, *_ = plothelper.plot(ss, show=True,
                              sphlines=[ckt.vout / ckt.lm * (ckt.lr + ckt.lm)])
    plt.close(fig)


def test8():
    ckt = ahbllc.AHBLLCTrafo(lr=25e-6, lm=1225e-6, cr=39e-9, nps=4.3, vbus=410, vload=80, chb=500e-12)
    con = .4
    v0 = 270
    sss = []
    for _ in range(100):
        _, ss = ahbllc.sim(v0, ckt, (con, 500e-9), constr=ahbllc.state_h0_cpm)
        v0 = ss[-1].v1
        print(f'v0 = {fmt(ss[0].v0)}V, actual dvoff = {fmt(ss[0].v1 - ss[0].v0, 4)}V')
        sss += ss
    fig, *_ = plothelper.plot(sss, show=True, sphlines=[ckt.vout / ckt.lm * (ckt.lr + ckt.lm)],
                            spradii=True, splegends=False)
    plt.close(fig)


def test9():
    ckt = ahbllc.AHBLLCTrafo(lr=25e-6, lm=1225e-6, cr=39e-9, nps=4.3, vbus=410, vload=20, chb=500e-12)
    con = 200e-9
    v0 = nsolve(lambda v: ahbllc.sim(v, ckt, (con, 500e-9))[0], -ckt.vbus, ckt.vout / ckt.lm * (ckt.lr + ckt.lm))
    _, ss = ahbllc.sim(v0, ckt, (con, 500e-9))
    fig, *_ = plothelper.plot(ss, show=True)
    ev = ahbllc.evaluate_switching_period(ss)
    print(ev.set_circuit(ckt))
    plt.close(fig)
    # for con in np.arange(1e-6, 4e-6, 1e-6):
    #     v0 = nsolve(lambda v: ahbllc.sim(v, ckt, (con, 500e-9))[0], -ckt.vbus, ckt.vout / ckt.lm * (ckt.lr + ckt.lm))
    #     _, ss = ahbllc.sim(v0, ckt, (con, 500e-9))
    #     fig, *_ = plothelper.plot(ss, show=True)
    #     plt.close(fig)


if __name__ == '__main__':
    test9()
