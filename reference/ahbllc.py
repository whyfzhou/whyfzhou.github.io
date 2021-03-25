import math
from scipy.optimize import brentq as nsolve
from plothelper import fmt


class AHBLLC:
    def __init__(self, lr, lm, cr, vbus, vout, chb=1e-12):
        self.lr = lr
        self.lm = lm
        self.cr = cr
        self.vbus = vbus
        self.vout = vout
        self.chb = chb

    def __str__(self):
        s = ''
        s += f'{" An asymmetrical half-bridge LLC converter ":-^40}\n'
        s += f'{"Lr":>20} = {fmt(self.lr)}H\n'
        s += f'{"Lm":>20} = {fmt(self.lm)}H\n'
        s += f'{"Cr":>20} = {fmt(self.cr)}F\n'
        s += f'{"Chb":>20} = {fmt(self.chb)}F\n'
        s += f'{"DC bus voltage":>20} = {fmt(self.vbus)}V\n'
        s += f'{"Output voltage":>20} = {fmt(self.vout)}V\n'
        return s


class AHBLLCTrafo(AHBLLC):
    def __init__(self, lr, lm, cr, nps, vbus, vload, chb=1e-12):
        self.lr = lr
        self.lm = lm
        self.cr = cr
        self.nps = nps
        self.vbus = vbus
        self.vload = vload
        self.vout = self.nps * self.vload
        self.chb = chb

    def __str__(self):
        s = super().__str__()
        s += f'{"nps":>20} = {fmt(self.nps)}\n'
        s += f'{"Load voltage":>20} = {fmt(self.vload)}V\n'
        return s


class Evaluation:
    def __init__(self, **kwargs):
        self.tsw = kwargs.get('tsw', 0)
        self.irms = kwargs.get('irms', 0)
        self.imax = kwargs.get('imax', 0)
        self.vavg = kwargs.get('vavg', 0)
        self.vrms = kwargs.get('vrms', 0)
        self.iout = kwargs.get('iout', 0)
        self.t_diode_on = kwargs.get('t_diode_on', 0)
        self.t_hs_off = kwargs.get('t_hs_off', 0)
        self.t_ls_on = kwargs.get('t_ls_on', 0)
        self.t_ls_off = kwargs.get('t_ls_off', 0)
        self.t_hs_on = kwargs.get('t_hs_on', 0)
        self.ckt = None

    def set_circuit(self, ckt):
        self.ckt = ckt

    def __str__(self):
        s = ''
        s += f'{"Tsw":>30} = {fmt(self.tsw, 4)}s\n'
        s += f'{"fsw":>30} = {fmt(1 / self.tsw, 4)}Hz\n'
        s += f'{"δhs = T_hs_on / Tsw":>30} ≈ {self.t_hs_on / self.tsw:.1%}\n'
        s += f'{"δd = T_diode_on / T_ls_on":>30} ≈ {self.t_diode_on / self.tsw:.1%}\n'
        s += f'{"Irms,pri":>30} = {fmt(self.irms, 3)}A\n'
        s += f'{"Imax,pri":>30} = {fmt(self.imax, 3)}A\n'
        s += f'{"Vdc,r":>30} = {fmt(self.vavg, 3)}V\n'
        s += f'{"Vacrms,r":>30} = {fmt((self.vrms**2 - self.vavg**2)**.5, 3)}V\n'
        if self.ckt is None:
            s += f'{"Idc,out":>30} = {fmt(self.iout, 3)}A\n'
        else:
            pout = self.iout * self.ckt.vout
            s += f'{"Vdc,bus":>30} = {fmt(self.ckt.vbus, 3)}V\n'
            s += f'{"Pavg,out":>30} = {fmt(pout, 3)}W\n'
            if isinstance(self.ckt, AHBLLCTrafo):
                iload = pout / self.ckt.vload
                s += f'{"Idc,load":>30} = {fmt(iload, 3)}A\n'
                s += f'{"Vdc,load":>30} = {fmt(self.ckt.vload, 3)}V\n'
            else:
                s += f'{"Idc,out":>30} = {fmt(self.iout, 3)}V\n'
                s += f'{"Vdc,out":>30} = {fmt(self.ckt.vout, 3)}V\n'
        return s


class State:
    def __init__(self, **kwargs):
        self.state = kwargs.get('state', 0)
        self.dt = kwargs.get('dt', 0)

        self.i0 = kwargs.get('i0', 0)
        self.v0 = kwargs.get('v0', 0)
        self.vhb0 = kwargs.get('vhb0', 0)
        self.im0 = kwargs.get('im0', 0)

        self.i1 = kwargs.get('i1', 0)
        self.v1 = kwargs.get('v1', 0)
        self.vhb1 = kwargs.get('vhb1', 0)
        self.im1 = kwargs.get('im1', 0)

        self.r = kwargs.get('r', 0)
        self.phi = kwargs.get('phi', 0)
        self.w = kwargs.get('w', 0)
        self.z = kwargs.get('z', 0)

        self.cr = kwargs.get('cr', 0)
        self.chb = kwargs.get('chb', 0)

        self.vcen = kwargs.get('vcen', 0)
        self.vout = kwargs.get('vout', 0)
        self.km = kwargs.get('km', 0)

    def __getitem__(self, key: str):
        return self.__dict__.get(key, 0)


def state_l1(i0, v0, vhb0, im0, ckt, con):
    """Sub-state of:
    l - low-side device on,
    1 - output rectifier on.

    Exit criteria: resonant current == magnetizing current (solve the Kepler's Equation)
    """
    del vhb0, con

    chb = math.inf
    cr = ckt.cr
    ctot = 1 / (1 / cr + 1 / chb)
    ltot = ckt.lr
    w = (ltot * ctot)**-.5
    z = (ltot / ctot)**.5
    km = ckt.vout / ckt.lm

    vhb0 = 0
    vout = ckt.vout
    vcen = vhb0 + vout
    r = math.hypot(v0 - vcen, i0 * z)
    phi = math.atan2(v0 - vcen, i0 * z)
    # i = r * cos(w * t + phi) / z
    # v = r * sin(w * t + phi) + vcen
    # vhb = 0
    # im = im0 - km * t

    ta = tb = (-phi % (2 * math.pi)) / w
    # ta = tb / 2
    # ta = tb = (im0 + r / z) / km
    while r * math.cos(w * ta + phi) / z > im0 - km * ta and ta > 1e-9:
        ta /= 2
    if ta > 1e-9:
        dt = nsolve(lambda t: im0 - km * t - r * math.cos(w * t + phi) / z, ta, tb)
    else:
        dt = ta
    i1 = r * math.cos(w * dt + phi) / z
    v1 = r * math.sin(w * dt + phi) + vout
    next_state = state_l0
    return next_state, State(state='l1', dt=dt,
                             i0=i0, v0=v0, vhb0=vhb0, im0=im0,
                             i1=i1, v1=v1, vhb1=vhb0, im1=i1,
                             r=r, phi=phi, w=w, z=z,
                             cr=cr, chb=chb,
                             vcen=vcen, vout=vout, km=-km) # yapf: disable


def state_l0(i0, v0, vhb0, im0, ckt, con):
    """Sub-state of:
    l - low-side device on,
    0 - output rectifier off.

    Exit criteria:
      (1) resonant voltage high enough to turn on the diode --> l1
      (2) delay timeout --> c0
    """
    del vhb0, im0
    t12 = con

    chb = math.inf
    cr = ckt.cr
    ctot = 1 / (1 / cr + 1 / chb)
    ltot = ckt.lr + ckt.lm
    w = (ltot * ctot)**-.5
    z = (ltot / ctot)**.5

    vhb0 = 0
    vout = 0
    vcen = vhb0 + vout
    r = math.hypot(v0 - vcen, i0 * z)
    phi = math.atan2(v0 - vcen, i0 * z)
    # i = r * cos(w * t + phi) / z
    # v = r * sin(w * t + phi) + vcen
    # vhb = 0
    # im = i

    v_diode_on = ckt.vout / ckt.lm * (ckt.lr + ckt.lm)
    # solve dt in v_diode_on == r * math.sin(w * dt + phi) + vg, or
    # solve i in (i * z)**2 + (v_diode_on - vg)**2 == r**2, and let dt_diode_on = (math.atan2(v_diode_on - vg, i * z) - phi) / w
    if -r <= v_diode_on - vcen <= r:
        dt_diode_on = min((math.asin((v_diode_on - vcen) / r) - phi) % (2 * math.pi),
                          (math.pi - math.asin((v_diode_on - vcen) / r) - phi) % (2 * math.pi)) / w
    else:
        dt_diode_on = math.inf

    i_t12 = r * math.cos(w * t12 + phi) / z
    t_i_to_0 = ((math.pi / 2 - phi) % (2 * math.pi)) / w
    if dt_diode_on < t12 or i_t12 > 0 and dt_diode_on < t_i_to_0:
        dt = dt_diode_on
        next_state = state_l1
    elif i_t12 <= 0:  # i <= 0 at t12
        dt = t12
        next_state = state_c0
    else:
        dt = t_i_to_0
        next_state = state_c0

    # if dt_diode_on < t12:  # diode turns on before t12 timeout
    #     dt = dt_diode_on
    #     next_state = state_l1
    # elif r * math.cos(w * t12 + phi) / z <= 0:  # i <= 0 at t12
    #     dt = t12
    #     next_state = state_c0
    # else:  # i > 0 at t12, postpone t12 to avoid capacitive mode
    #     dt = ((math.pi / 2 - phi) % (2 * math.pi)) / w
    #     next_state = state_c0

    i1 = r * math.cos(w * dt + phi) / z
    v1 = r * math.sin(w * dt + phi) + vcen
    return next_state, State(state='l0', dt=dt,
                             i0=i0, v0=v0, vhb0=vhb0, im0=i0,
                             i1=i1, v1=v1, vhb1=vhb0, im1=i1,
                             r=r, phi=phi, w=w, z=z,
                             cr=cr, chb=chb,
                             vcen=vcen, vout=vout, km=0) # yapf: disable


def state_h0(i0, v0, vhb0, im0, ckt, con):
    """Sub-state of:
    h - high-side device on,
    0 - output rectifier off.

    Exit criteria: high-side turned off when v reaches voff
    """
    del vhb0, im0
    voff = con

    chb = math.inf
    cr = ckt.cr
    ctot = 1 / (1 / cr + 1 / chb)
    ltot = ckt.lr + ckt.lm
    w = (ltot * ctot)**-.5
    z = (ltot / ctot)**.5

    vhb0 = ckt.vbus
    vout = 0
    vcen = vhb0 + vout
    r = math.hypot(v0 - vcen, i0 * z)
    phi = math.atan2(v0 - vcen, i0 * z)
    # i = r * cos(w * t + phi) / z
    # v = r * sin(w * t + phi) + vcen
    # vhb = vbus
    # im = i

    # now solve the following system
    #   (i * z)**2 + (v - vcen)**2 == r**2
    #   v == voff
    if -r <= voff - vcen <= r:
        i1 = (r**2 - (voff - vcen)**2)**.5 / z
        v1 = voff
    else:
        i1 = 0
        v1 = vcen + r
    dt = (math.atan2(v1 - vcen, i1 * z) - phi) / w
    v_diode_on = ckt.vout / ckt.lm * (ckt.lr + ckt.lm)
    if v1 > v_diode_on:
        next_state = state_c1
    else:
        next_state = state_c0
    return next_state, State(state='h0', dt=dt,
                             i0=i0, v0=v0, vhb0=vhb0, im0=i0,
                             i1=i1, v1=v1, vhb1=vhb0, im1=i1,
                             r=r, phi=phi, w=w, z=z,
                             cr=cr, chb=chb,
                             vcen=vcen, vout=vout, km=0) # yapf: disable


def state_c0(i0, v0, vhb0, im0, ckt, con):
    """Sub-state of:
    c - both devices are off, half-bridge behaving like a capacitor,
    0 - output rectifier off,
    This state is only needed for considering the commutation transient.

    Exit criteria:
      (1) half-bridge voltage reaches DC bus --> h0
      (2) half-bridge voltage reaches 0 --> l0
      (3) magnetizing inductor voltage high enough to turn on the diode --> c1
          only possible in high-side to low-side commutation
    """
    del im0, con

    chb = ckt.chb
    cr = ckt.cr
    ctot = 1 / (1 / cr + 1 / chb)
    ltot = ckt.lr + ckt.lm
    w = (ltot * ctot)**-.5
    z = (ltot / ctot)**.5

    vcen = vout = 0
    vcap0 = v0 - vhb0
    r = math.hypot(vcap0 - vcen, i0 * z)
    phi = math.atan2(vcap0 - vcen, i0 * z)
    # i(t) == r * cos(w * t + phi) / z, resonant current
    #   i(0) == i0 <= 0, guaranteed by previous state
    # vcap(t) = r * sin(w * t + phi) + vcen, voltage across two capacitors
    # v(t) == r * sin(w * t + phi) / (1 + cr / chb) +
    #         v0 / (1 + chb / cr) +
    #         (vhb0 + vout) / (1 + cr / chb)
    # vhb(t) == -r * sin(w * t + phi) / (1 + chb / cr) -
    #           (vout - v0) / (1 + chb / cr) +
    #           vhb0 / (1 + cr / chb)

    if -1e-6 < vhb0 < 1e-6:
        v = ckt.vbus
        v -= -(vout - v0) / (1 + chb / cr) + vhb0 / (1 + cr / chb)
        v *= (1 + chb / cr)
        if -r <= v <= r:
            dt = min((math.asin(v / -r) - phi) % (2 * math.pi),
                     (math.pi - math.asin(v / -r) - phi) % (2 * math.pi)) / w
        else:
            dt = (-math.pi / 2 - phi) % (2 * math.pi) / w
        next_state = state_h0
    else:
        v = 0
        v -= -(vout - v0) / (1 + chb / cr) + vhb0 / (1 + cr / chb)
        v *= (1 + chb / cr)
        if -r <= v <= r:
            dt1 = min((math.asin(v / -r) - phi) % (2 * math.pi),
                      (math.pi - math.asin(v / -r) - phi) % (2 * math.pi)) / w
        else:
            dt1 = (math.pi / 2 - phi) % (2 * math.pi) / w

        v = ckt.vout / ckt.lm * (ckt.lr + ckt.lm)
        v -= vcen
        if -r <= v <= r:
            dt2 = min((math.asin(v / r) - phi) % (2 * math.pi),
                      (math.pi - math.asin(v / r) - phi) % (2 * math.pi)) / w
        else:
            dt2 = math.inf

        if dt2 < dt1:
            dt = dt2
            next_state = state_c1
        else:
            dt = dt1
            next_state = state_l0

    i1 = r * math.cos(w * dt + phi) / z
    v1 = (r * math.sin(w * dt + phi) / (1 + cr / chb) +
          v0 / (1 + chb / cr) +
          (vhb0 + vout) / (1 + cr / chb))
    vhb1 = (-r * math.sin(w * dt + phi) / (1 + chb / cr) -
            (vout - v0) / (1 + chb / cr) +
            vhb0 / (1 + cr / chb))
    if -1e-6 < vhb1 - ckt.vbus < 1e-6:
        vhb1 = ckt.vbus
    elif -1e-6 < vhb1 < 1e-6:
        vhb1 = 0

    return next_state, State(state='c0', dt=dt,
                             i0=i0, v0=v0, vhb0=vhb0, im0=i0,
                             i1=i1, v1=v1, vhb1=vhb1, im1=i1,
                             r=r, phi=phi, w=w, z=z,
                             cr=cr, chb=chb,
                             vcen=vcen, vout=vout, km=0) # yapf: disable


def state_c1(i0, v0, vhb0, im0, ckt, con):
    """Sub-state of:
    c - both devices are off, half-bridge behaving like a capacitor,
    1 - output rectifier on.
    This state occurs when
      1) we consider the commutation, and
      2) in heavy load where the output diode turns on immediately after the turning-off of the high-side MOS.

    Exit criteria:
      (1) half-bridge voltage reaches DC bus or 0 --> h1 (prevented in control strategy) or l1
      (2) Exit criteria: resonant current == magnetizing current (solve the Kepler's Equation) --> c0
          very unlikely to happen (but it's *not* impossible, consider the condition where Lm is extremely small)
    """
    del con

    chb = ckt.chb
    cr = ckt.cr
    ctot = 1 / (1 / cr + 1 / chb)
    ltot = ckt.lr
    w = (ltot * ctot)**-.5
    z = (ltot / ctot)**.5

    vcen = vout = ckt.vout
    vcap0 = v0 - vhb0
    r = math.hypot(vcap0 - vcen, i0 * z)
    phi = math.atan2(vcap0 - vcen, i0 * z)
    km = ckt.vout / ckt.lm
    # i(t) == r * cos(w * t + phi) / z, resonant current
    #   i(0) == i0 <= 0, guaranteed by previous state
    # vcap(t) = r * sin(w * t + phi) + vcen, voltage across two capacitors
    # v(t) == r * sin(w * t + phi) / (1 + cr / chb) +
    #         v0 / (1 + chb / cr) +
    #         (vhb0 + vout) / (1 + cr / chb)
    # vhb(t) == -r * sin(w * t + phi) / (1 + chb / cr) -
    #           (vout - v0) / (1 + chb / cr) +
    #           vhb0 / (1 + cr / chb)

    v = 0
    v -= -(vout - v0) / (1 + chb / cr) + vhb0 / (1 + cr / chb)
    v *= (1 + chb / cr)
    if -r <= v <= r:
        dt1 = min((math.asin(v / -r) - phi) % (2 * math.pi),
                  (math.pi - math.asin(v / -r) - phi) % (2 * math.pi)) / w
    else:
        dt1 = (math.pi / 2 - phi) % (2 * math.pi) / w  # hard-switching, babe~

    ta = tb = (-phi % (2 * math.pi)) / w
    # ta = tb / 2
    # ta = tb = (im0 + r / z) / km
    while r * math.cos(w * ta + phi) / z > im0 - km * ta and ta > 1e-9:
        ta /= 2
    if ta > 1e-9:
        dt2 = nsolve(lambda t: im0 - km * t - r * math.cos(w * t + phi) / z, ta, tb)
    else:
        dt2 = ta

    if dt2 < dt1:  # diode turns off before vhb reaches 0
        dt = dt2
        next_state = state_c0
    else:
        dt = dt1
        next_state = state_l1

    i1 = r * math.cos(w * dt + phi) / z
    v1 = (r * math.sin(w * dt + phi) / (1 + cr / chb) +
          v0 / (1 + chb / cr) +
          (vhb0 + vout) / (1 + cr / chb))
    vhb1 = (-r * math.sin(w * dt + phi) / (1 + chb / cr) -
            (vout - v0) / (1 + chb / cr) +
            vhb0 / (1 + cr / chb))
    if -1e-6 < vhb1 - ckt.vbus < 1e-6:
        vhb1 = ckt.vbus
    elif -1e-6 < vhb1 < 1e-6:
        vhb1 = 0
    im1 = im0 - km * dt
    if -1e-6 < im1 - i1 < 1e-6:
        im1 = i1

    return next_state, State(state='c1', dt=dt,
                             i0=i0, v0=v0, vhb0=vhb0, im0=im0,
                             i1=i1, v1=v1, vhb1=vhb1, im1=im1,
                             r=r, phi=phi, w=w, z=z,
                             cr=cr, chb=chb,
                             vcen=vcen, vout=vout, km=-km) # yapf: disable


def _sim1(isf, cond, i0, v0, vhb0, im0, ckt, t12):
    nsf, s = isf(i0, v0, vhb0, im0, ckt, t12)
    states = [s]
    while cond(nsf):
        nsf, s = nsf(s.i1, s.v1, s.vhb1, s.im1, ckt, t12)
        states.append(s)
    return nsf, states


def sim(i0, ckt, con):
    commutating = lambda s: s in {state_c0, state_c1}
    active = lambda s: s not in {state_c0, state_c1}

    def eq_zvon(t, hs_off_fins):
        _, ss = _sim1(isf_ls_on, active, *hs_off_fins, ckt, t)
        ls_on_fins = (ss[-1].i1, ss[-1].v1, ss[-1].vhb1, ss[-1].im1)
        _, ss = _sim1(state_c0, commutating, *ls_on_fins, ckt, t)
        c0 = ss[-1]
        r, chb, cr = c0.r, c0.chb, c0.cr
        vout, v0, vhb0 = c0['vout'], c0.v0, c0.vhb0
        vhbmax = r / (1 + chb/cr) - (vout-v0) / (1 + chb/cr) + vhb0 / (1 + cr/chb)
        return vhbmax - ckt.vbus

    # def eq_vcont(t, hs_off_fins):
    #     """实际上并不是这么做的，而是用了限制最高频率的方式。"""
    #     _, ss = sim1(isf_ls_on, active, *hs_off_fins, ckt, t)
    #     ls_on_fins = (ss[-1].i1, ss[-1].v1, ss[-1].vhb1, ss[-1].im1)
    #     _, ss = sim1(state_c0, commutating, *ls_on_fins, ckt, t)
    #     # ls_off_fin_v = ss[-1].v1
    #     # return ls_off_fin_v - voff*.99  # FIXME: 实际上并不是这么做的
    #     v0 = ss[-1].v1
    #     i0 = ss[-1].i1
    #     z = ((ckt.lr + ckt.lm) / ckt.cr)**.5
    #     vcen = ckt.vbus
    #     r = math.hypot(v0 - vcen, i0 * z)
    #     # phi = math.atan2(v0 - vcen, i0 * z)
    #     # v = r * sin(w * t + phi) + vcen
    #     vmin = -r + vcen
    #     return vmin - voff*.99  # FIXME: incorrect, leading to multiple root in finding i0, 会改变 sim 方程的单调性

    voff, t12min = con
    # starting at the turning-off of the high-side device
    v0 = voff
    vhb0 = ckt.vbus
    v_diode_on = ckt.vout / ckt.lm * (ckt.lr + ckt.lm)
    isf = state_c1 if v0 - vhb0 >= v_diode_on else state_c0

    nsf, hs_off = _sim1(isf, commutating, i0, v0, vhb0, i0, ckt, t12min)
    hs_off_fins = (hs_off[-1].i1, hs_off[-1].v1, hs_off[-1].vhb1, hs_off[-1].im1)

    isf_ls_on = nsf
    nsf, ls_on = _sim1(isf_ls_on, active, *hs_off_fins, ckt, t12min)
    ls_on_fins = (ls_on[-1].i1, ls_on[-1].v1, ls_on[-1].vhb1, ls_on[-1].im1)
    t12capm = ls_on[-1].dt

    nsf, ls_off = _sim1(state_c0, commutating, *ls_on_fins, ckt, t12min)
    ls_off_fins = (ls_off[-1].i1, ls_off[-1].v1, ls_off[-1].vhb1, ls_off[-1].im1)

    t12zvon = 0
    # t12zvon = t12vcont = 0
    t12max = (math.pi - ls_on[-1].phi) / ls_on[-1].w
    if ls_off[-1].vhb1 < ckt.vbus - 1e-6:
        if eq_zvon(t12max, hs_off_fins) >= 0:
            t12zvon = nsolve(lambda t: eq_zvon(t, hs_off_fins), t12min / 2, t12max)
        else:
            t12zvon = t12max
    # if ls_off[-1].v1 > voff:
    #     t12vcont = nsolve(lambda t: eq_vcont(t, hs_off_fins), t12min / 2, t12max)
    # t12 = max(t12min, t12capm, t12zvon, t12vcont)
    t12 = max(t12min, t12capm, t12zvon)

    if t12 > t12min:
        nsf, ls_on = _sim1(isf_ls_on, active, *hs_off_fins, ckt, t12)
        ls_on_fins = (ls_on[-1].i1, ls_on[-1].v1, ls_on[-1].vhb1, ls_on[-1].im1)
        nsf, ls_off = _sim1(state_c0, commutating, *ls_on_fins, ckt, t12)
        ls_off_fins = (ls_off[-1].i1, ls_off[-1].v1, ls_off[-1].vhb1, ls_off[-1].im1)

    if ls_off[-1].v1 > voff:
        hs_on_v0 = ls_off[-1].v1
        hs_on_i0 = ls_off[-1].i1
        hs_on_z = ((ckt.lr + ckt.lm) / ckt.cr)**.5
        hs_on_w = ((ckt.lr + ckt.lm) * ckt.cr)**-.5
        hs_on_vcen = ckt.vbus
        hs_on_r = math.hypot(hs_on_v0 - hs_on_vcen, hs_on_i0 * hs_on_z)
        hs_on_phi = math.atan2(hs_on_v0 - hs_on_vcen, hs_on_i0 * hs_on_z)
        dtmax = 600e-9
        dthson = dtmax - sum(s.dt for s in ls_off)
        if dthson <= 100e-9:
            dthson = 100e-9
        hs_on_v1 = hs_on_r * math.sin(hs_on_w * dthson + hs_on_phi) + hs_on_vcen
        _, hs_on = _sim1(state_h0, active, *ls_off_fins, ckt, hs_on_v1)
    else:
        _, hs_on = _sim1(state_h0, active, *ls_off_fins, ckt, voff)
    # if hs_on[-1].v1 > voff + 1e-6:
    #     res = hs_on[-1].v1 - voff
    # else:
    #     res = hs_on[-1].i1 - i0
    # _, hs_on = sim1(state_h0, active, *ls_off_fins, ckt, voff)
    res = hs_on[-1].i1 - i0
    states = hs_off + ls_on + ls_off + hs_on
    # return hs_on[-1].i1 - i0, t12, states
    return res, t12, states


def find_steady_state(voff, ckt, t12min=500e-9, fswmax=100e3):
    def eq(i):
        di, _, ss = sim(i, ckt, (voff, t12min))
        fsw = 1 / sum(s.dt for s in ss)
        if fsw > fswmax:
            def eqf(t):
                _, _, ss = sim(i, ckt, (voff, t))
                fsw = 1 / sum(s.dt for s in ss)
                return fsw - fswmax
            t12 = nsolve(eqf, .5e-6, 1 / fswmax)
            di, _, ss = sim(i, ckt, (voff, t12))
        else:
            t12 = max(t12min, [s for s in ss if s.state == 'l0'][-1].dt)
        dv = ss[-1].v1 - voff
        return di, dv, t12, ss

    i0max = ckt.vbus / ((ckt.lr + ckt.lm) / ckt.cr)**.5
    i0 = nsolve(lambda i: eq(i)[0], 1e-3, i0max * 2)
    _, _, t12, ss = eq(i0)
    return t12, ss


def evaluate_steady_state(states):
    tsw = 0
    qout = i2acc = imax = vacc = v2acc = 0
    t_diode_on = t_hs_off = t_ls_on = t_ls_off = t_hs_on = 0
    for state in states:
        tsw += state.dt
        j = state.r / state.z
        i2acc += j**2 * state.dt / 2
        i2acc += j**2 / (4 * state.w) * (math.sin(2 * state.w * state.dt + 2 * state.phi) - math.sin(2 * state.phi))
        m = state.r / (1 + state.cr / state.chb)
        b = state.v0 / (1 + state.chb / state.cr)
        b += (state.vhb0 + state['vout']) / (1 + state.cr / state.chb)
        vacc += m / state.w * (math.cos(state.phi) - math.cos(state.w * state.dt + state.phi))
        vacc += b * state.dt
        v2acc += (m**2 / 2 + b**2) * state.dt
        v2acc -= m**2 / (4 * state.w) * (math.sin(2 * state.w * state.dt + 2 * state.phi) - math.sin(2 * state.phi))
        v2acc += 2 * m * b / state.w * (math.cos(state.phi) - math.cos(state.w * state.dt + state.phi))
        if state.state.endswith('1'):
            t_diode_on += state.dt
            tr = (state.im0 + state.im1) * state.dt / 2
            sc = state.cr * (state.v0 - state.v1)
            qout += tr + sc
        if state.state.startswith('c'):
            if state.vhb0 > state.vhb1:
                t_hs_off += state.dt
            else:
                t_ls_off += state.dt
        if state.state.startswith('l'):
            t_ls_on += state.dt
        if state.state.startswith('h'):
            t_hs_on += state.dt
        if state.i1 > imax:
            imax = state.i1

    irms = (i2acc / tsw)**.5
    vavg = vacc / tsw
    vrms = (v2acc / tsw)**.5
    iout = qout / tsw
    return Evaluation(tsw=tsw,
                      irms=irms,
                      imax=imax,
                      vavg=vavg,
                      vrms=vrms,
                      iout=iout,
                      t_diode_on=t_diode_on,
                      t_hs_off=t_hs_off,
                      t_ls_on=t_ls_on,
                      t_ls_off=t_ls_off,
                      t_hs_on=t_hs_on)


def evaluate_operating_point(pout, ckt, t12min=500e-9, fswmax=100e3):
    pmax = evaluate_steady_state(find_steady_state(ckt.vbus, ckt, t12min, fswmax)[-1]).iout * ckt.vout
    if 0 < pout <= pmax:
        voff = nsolve(lambda v: evaluate_steady_state(find_steady_state(v, ckt, t12min, fswmax)[-1]).iout * ckt.vout - pout,
                      1e-3, ckt.vbus)
        t12, ss = find_steady_state(voff, ckt, t12min, fswmax)
        return t12, ss, evaluate_steady_state(ss)
    else:
        return 0, [State()], Evaluation()
