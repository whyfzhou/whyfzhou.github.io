import math
import itertools
from scipy.optimize import ridder as nsolve
import ahbllc
import plothelper
import matplotlib.pyplot as plt


def evaluate_operating_point_with_tolerance(pout, ckt, t12min=500e-9, fswmax=100e3):
    assert isinstance(ckt, ahbllc.AHBLLCTrafo)
    minmax = lambda x, rtol: (x * (1 - rtol), x * (1 + rtol))
    lrs = minmax(ckt.lr, ckt.lr_rtol)
    lms = minmax(ckt.lm, ckt.lm_rtol)
    crs = minmax(ckt.cr, ckt.cr_rtol)
    chbs = minmax(ckt.chb, ckt.chb_rtol)
    npss = minmax(ckt.nps, ckt.nps_rtol)
    vbuss = minmax(ckt.vbus, ckt.vbus_rtol)
    vload = ckt.vload

    first = True
    tsw = 0
    sss = []
    evs = []
    for lr, lm, cr, chb, nps, vbus in itertools.product(lrs, lms, crs, chbs, npss, vbuss):
        vout = nps * vload
        cc = ahbllc.AHBLLC(lr=lr, lm=lm, cr=cr, vbus=vbus, vout=vout, chb=chb)
        _, ss, ev = ahbllc.evaluate_operating_point(pout, cc, t12min, fswmax)
        if ev.tsw > tsw:
            tsw = ev.tsw
        if first:
            fig, axes = plothelper.plot(ss, show=False, spradii=False, splegends=False)
            first = False
        else:
            plothelper.plot(ss, fig=fig, axes=axes, show=False, spradii=False, splegends=False)
        ev.set_circuit(cc)
        sss.append(ss)
        evs.append(ev)
        print(cc)
        print(ev)
        print('*' * 120)
        print()
    axes[0].set_xlim((0, tsw))
    plt.show()