import math
import itertools
from scipy.optimize import ridder as nsolve
import ahbllc
from plothelper import (fmt, plot)
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
    failed = []
    failed_pmax = []
    for iter, (vbus, lr, lm, cr, chb, nps) in enumerate(itertools.product(vbuss, lrs, lms, crs, chbs, npss)):
        vout = nps * vload
        cc = ahbllc.AHBLLC(lr=lr, lm=lm, cr=cr, vbus=vbus, vout=vout, chb=chb)
        print(cc)
        tf, ss, ev, pmax = ahbllc.evaluate_operating_point(pout, cc, t12min, fswmax)
        if tf > 0:
            print('{:^60}\n'.format(f'{" Solved: #" + str(iter) + " ":_^30}'))
            if ev.tsw > tsw:
                tsw = ev.tsw
            if first:
                fig, axes = plot(ss, show=False, spradii=False, splegends=False)
                first = False
            else:
                plot(ss, fig=fig, axes=axes, show=False, spradii=False, splegends=False)
            ev.set_circuit(cc)
            sss.append(ss)
            evs.append(ev)
            print(ev)
            print('*' * 120)
            print()
        else:
            print('{:^60}\n'.format(f'{" Failed: #" + str(iter) + " ":_^30}'))
            print('{:^60}\n'.format(f'{" Maximum attainable output power: #" + fmt(pmax) + "W ":_^40}'))
            failed.append(cc)
            failed_pmax.append(pmax)

    axes[0].set_xlim((0, tsw))
    plt.show()

    if failed:
        print(f'Totally {len(failed)} failed condition(s):')
        for f, p in zip(failed, failed_pmax):
            print(f)
            print(f'Maximum attainable output power: {fmt(p)}W.')
    else:
        print("All boundary conditions passed!")