import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FuncFormatter


plt.rcParams.update({'font.family':'Fira Code'})
plt.rcParams.update({'pdf.fonttype':42})
plt.rcParams.update({'ps.fonttype':42})


def fmt(x, n=3):
    if x == 0 or x == -0.0 or x == 0.0:
        return '0'

    m = math.log10(abs(x)) // 3
    b = x / 10**(3*m)
    metric_prefix = {
        8: 'Y',
        7: 'Z',
        6: 'E',
        5: 'P',
        4: 'T',
        3: 'G',
        2: 'M',
        1: 'k',
        0: '',
        -1: 'm',
        -2: 'μ',
        -3: 'n',
        -4: 'p',
        -5: 'f',
        -6: 'a',
        -7: 'z',
        -8: 'y',
    }.get(m, f'e{3*m:.0f}')
    f = f'{{:.{max(0, n - math.floor(math.log10(abs(b))) - 1)}f}}'

    s = f.format(b) + metric_prefix
    return s


@FuncFormatter
def axis_formatter(x, pos):
    return fmt(x)


def plot(ss, iter=None, fig=None, axes=None, show=True, filename=None, highlight=False, sphlines=None, spradii=True, splegends=True):
    n = 500
    if axes is None:
        fig = plt.figure(figsize=(8, 8))
        # fig = plt.figure(figsize=(6, 6))
        gs = GridSpec(2, 2, figure=fig)
        ax_i = fig.add_subplot(gs[0, 0])
        ax_v = fig.add_subplot(gs[1, 0], sharex=ax_i)
        ax_sp = fig.add_subplot(gs[:, 1])
    else:
        ax_i, ax_v, ax_sp = axes

    t = 0
    last_state = None
    lss = ['solid', 'dashed', 'dashdot', 'dotted']
    for k, state in enumerate(ss):
        c = plt.rcParams['axes.prop_cycle'].by_key()['color'][k % len(plt.rcParams['axes.prop_cycle'].by_key()['color'])]
        s1 = state.state[0].upper()
        s2 = state.state[1]
        if iter is None:
            sn = f'$\\rm{s1}_{{{s2}}}$'
        else:
            sn = f'{iter}: $\\rm{s1}_{{{s2}}}$'

        ts = np.linspace(0, state.dt, n)

        iseq = state.r * np.cos(state.w * ts + state.phi) / state.z
        vseq = (state.r * np.sin(state.w * ts + state.phi) / (1 + state.cr / state.chb) +
                state.v0 / (1 + state.chb / state.cr) +
                (state.vhb0 + state['vout']) / (1 + state.cr / state.chb))
        vhbseq = (-state.r * np.sin(state.w * ts + state.phi) / (1 + state.chb / state.cr) -
                  (state['vout'] - state.v0) /(1 + state.chb / state.cr) +
                  state.vhb0 / (1 + state.cr / state.chb))
        ax_i.plot(t + ts, iseq, color=c, label=sn)
        ax_v.plot(t + ts, vseq, color=c, label=sn)
        ax_v.plot(t + ts, vhbseq, color=c)
        if last_state is not None and (last_state[0].upper() == 'L'
                                       and state.state[0].upper() == 'H'
                                       or last_state[0].upper() == 'H'
                                       and state.state[0].upper() == 'L'):
            ax_v.plot([t, t], [last_state, state['hb']], color=c)
        last_state = state.state

        ax_sp.plot(iseq, vseq, label=sn)
        if spradii:
            ax_sp.plot([0, iseq[0]], [state.vcen, vseq[0]], color = '#aaa', linestyle=lss[k % 4], zorder=-10)
            ax_sp.plot([0, iseq[-1]], [state.vcen, vseq[-1]], color = '#aaa', linestyle=lss[k % 4], zorder=-10)

        if state.state[1] == '1':
            imseq = state.im0 + state.km * ts
            ax_i.plot(t + ts, imseq, linestyle=':', color=c)

        t += state.dt

    if highlight:
        ax_i.axhline(ss[0].i0, alpha=.5, c='#f00', linestyle='dashed')
        ax_i.axhline(ss[-1].i1, alpha=.5, c='#0f0', linestyle='dashdot')
    if sphlines is not None:
        for sphline in sphlines:
            ax_sp.axhline(sphline, alpha=.5, c='#06a', linestyle='dashed', zorder=-10, linewidth=2)

    ax_i.set_xlim((0, t))
    ax_i.grid(True)
    ax_i.set_ylabel('currents')
    ax_i.xaxis.set_major_formatter(axis_formatter)
    ax_i.xaxis.set_minor_formatter(axis_formatter)
    ax_v.set_xlim((0, t))
    ax_v.grid(True)
    ax_v.set_ylabel('voltages')
    ax_v.set_xlabel('time')
    if splegends:
        ax_sp.legend()
    ax_sp.grid(True)
    ax_sp.set_xlabel('current')
    ax_sp.set_ylabel('voltage')

    fig.tight_layout()
    if filename is not None:
        plt.savefig(filename)
    if show:
        plt.show()

    return fig, (ax_i, ax_v, ax_sp)

