import math
from numpy import result_type
from scipy.optimize import brentq as nsolve
from plothelper import (fmt, plot)


MINIMUM_TIME = 1e-9
MINIMUM_VOLTAGE = 1e-5
MINIMUM_CURRENT = 1e-7
MAXIMUM_COMMUTATION_TIME = 600e-9
MINIMUM_HIGH_SIDE_ON_TIME = 100e-9
MINIMUM_FORWARD_TIME = .5e-6
MAXIMUM_T12 = 20e-6
MINIMUM_DETECTABLE_DIODE_ON_TIME = 500e-9


class AHBLLC:
    def __init__(self, **kwargs):
        self.lr = kwargs.get('lr', 30e-6)
        self.lm = kwargs.get('lm', 1000e-6)
        self.cr = kwargs.get('cr', 47e-9)
        self.vbus = kwargs.get('vbus', 410)
        self.vout = kwargs.get('vout', 200)
        self.chb = kwargs.get('chb', 100e-12)

        self.lr_rtol = kwargs.get('lr_rtol', .05)
        self.lm_rtol = kwargs.get('lm_rtol', .05)
        self.cr_rtol = kwargs.get('cr_rtol', .05)
        self.vbus_rtol = kwargs.get('vbus_rtol', .05)
        self.chb_rtol = kwargs.get('chb_rtol', .05)

    def __str__(self):
        s = ''
        s += f'{" An asymmetrical half-bridge LLC converter ":-^60}\n'
        s += f'{"Lr":>28} = {fmt(self.lr)}H\n'
        s += f'{"Lm":>28} = {fmt(self.lm)}H\n'
        s += f'{"Cr":>28} = {fmt(self.cr)}F\n'
        s += f'{"Chb":>28} = {fmt(self.chb)}F\n'
        s += f'{"DC bus voltage":>28} = {fmt(self.vbus)}V\n'
        s += f'{"Output voltage":>28} = {fmt(self.vout)}V\n'
        return s


class AHBLLCTrafo(AHBLLC):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.vload = kwargs.get('vload', 40)
        self.nps = kwargs.get('nps', 3.8)
        self.vout = self.nps * self.vload

        self.nps_rtol = kwargs.get('nps_rtol', .10)

    def __str__(self):
        s = super().__str__()
        s += f'{"nps":>28} = {fmt(self.nps)}\n'
        s += f'{"Load voltage":>28} = {fmt(self.vload)}V\n'
        return s


class Evaluation:
    def __init__(self, **kwargs):
        self.operation = kwargs.get('operation', 'unknown')
        self.tsw = kwargs.get('tsw', 0)
        self.irms = kwargs.get('irms', 0)
        self.irms_out = kwargs.get('irms_out', 0)
        self.imax = kwargs.get('imax', 0)
        self.imax_out = kwargs.get('imax_out', 0)
        self.vavg = kwargs.get('vavg', 0)
        self.vrms = kwargs.get('vrms', 0)
        self.iout = kwargs.get('iout', 0)
        self.t_diode_on = kwargs.get('t_diode_on', 0)
        self.t_hs_off = kwargs.get('t_hs_off', 0)
        self.t_ls_on = kwargs.get('t_ls_on', 0)
        self.t_ls_off = kwargs.get('t_ls_off', 0)
        self.t_hs_on = kwargs.get('t_hs_on', 0)
        self.t12 = kwargs.get('t12', 0)
        self.ckt = None

    def set_circuit(self, ckt):
        self.ckt = ckt
        return self

    def __str__(self):
        s = ''
        if self.operation == 'unknown':
            s += f'{" Failed to evaluate, operation unknown ":-^60}\n'
        elif self.operation == 'steady-state':
            s += f'{" Converter operates in steady-state ":-^60}\n'
        else:
            s += f'{" Converter stabilizing... ":-^60}\n'
        s += f'{"Tsw":>28} = {fmt(self.tsw, 4)}s\n'
        s += f'{"fsw":>28} = {fmt(1 / self.tsw, 4)}Hz\n'
        s += f'{"(Corrected) T12":>28} = {fmt(self.t12, 4)}s\n'
        s += f'{"δhs = T_hs_on / Tsw":>28} ≈ {self.t_hs_on / self.tsw:.1%}\n'
        s += f'{"δd = T_diode_on / T_ls_on":>28} ≈ {self.t_diode_on / self.tsw:.1%}\n'
        s += f'{"Irms,pri":>28} = {fmt(self.irms, 3)}A\n'
        s += f'{"Imax,pri":>28} = {fmt(self.imax, 3)}A\n'
        s += f'{"Irms,out":>28} = {fmt(self.irms_out, 3)}A\n'
        s += f'{"Imax,out":>28} = {fmt(self.imax_out, 3)}A\n'
        s += f'{"Vdc,r":>28} = {fmt(self.vavg, 3)}V\n'
        s += f'{"Vacrms,r":>28} = {fmt((self.vrms**2 - self.vavg**2)**.5, 3)}V\n'
        if self.ckt is None:
            s += f'{"Idc,out":>28} = {fmt(self.iout, 3)}A\n'
        else:
            s += f'{"Ψ,pri":>28} = {fmt((self.ckt.lr + self.ckt.lm) * self.imax, 3)}Wb\n'
            pout = self.iout * self.ckt.vout
            s += f'{"Vdc,bus":>28} = {fmt(self.ckt.vbus, 3)}V\n'
            s += f'{"Pavg,out":>28} = {fmt(pout, 3)}W\n'
            if isinstance(self.ckt, AHBLLCTrafo):
                iload = pout / self.ckt.vload
                s += f'{"Idc,load":>28} = {fmt(iload, 3)}A\n'
                s += f'{"Vdc,load":>28} = {fmt(self.ckt.vload, 3)}V\n'
                s += f'{"Irms,diode":>28} = {fmt(self.irms_out * self.ckt.nps, 3)}A\n'
                s += f'{"Imax,diode":>28} = {fmt(self.imax_out * self.ckt.nps, 3)}A\n'
            else:
                s += f'{"Idc,out":>28} = {fmt(self.iout, 3)}A\n'
                s += f'{"Vdc,out":>28} = {fmt(self.ckt.vout, 3)}V\n'
        return s


class State:
    def __init__(self, **kwargs):
        self.state = kwargs.get('state', 'unknown')
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
        # for backward compatibility
        return self.__dict__.get(key, 0)


def state_l1(i0, v0, vhb0, im0, ckt, con):
    """Sub-state of:
    l - low-side device on,
    1 - output rectifier on.

    Exit criteria: resonant current == magnetizing current (solve the Kepler's Equation)
    """
    del vhb0, con  # vhb == 0

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
    while r * math.cos(w * ta + phi) / z > im0 - km * ta and ta > MINIMUM_TIME:
        ta /= 2
    if ta > MINIMUM_TIME:
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
    del vhb0, im0  # vhb == 0, im == i
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
    if -r <= v_diode_on - vcen <= r:
        dt_diode_on = min((math.asin((v_diode_on - vcen) / r) - phi) % (2 * math.pi),
                          (math.pi - math.asin((v_diode_on - vcen) / r) - phi) % (2 * math.pi)) / w
    else:
        dt_diode_on = math.inf

    i_t12 = r * math.cos(w * t12 + phi) / z
    dt_izc = ((math.pi / 2 - phi) % (2 * math.pi)) / w
    if dt_diode_on < t12 or i_t12 > 0 and dt_diode_on < dt_izc:
        dt = dt_diode_on
        next_state = state_l1
    elif i_t12 <= 0:
        dt = t12
        next_state = state_c0
    else:
        dt = dt_izc
        next_state = state_c0

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

    tfwd - fixed forward on-time control
    """
    del vhb0, im0
    dt = con

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

    if dt is None:
        # 0 = r * cos(w * t + phi) / z, w * t + phi = -pi / 2
        dt = ((-math.pi / 2 - phi) % (2 * math.pi)) / w
    i1 = r * math.cos(w * dt + phi) / z
    v1 = r * math.sin(w * dt + phi) + vcen

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
    del im0, con  # im == i

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
    # i(t) = r * cos(w * t + phi) / z, resonant current
    #   i(0) = i0 <= 0, guaranteed by previous state
    # vcap(t) = r * sin(w * t + phi) + vcen, voltage across two capacitors
    # v(t) = r * sin(w * t + phi) / (1 + cr / chb) +
    #        v0 / (1 + chb / cr) +
    #        (vhb0 + vout) / (1 + cr / chb)
    # vhb(t) = -r * sin(w * t + phi) / (1 + chb / cr) -
    #          (vout - v0) / (1 + chb / cr) +
    #          vhb0 / (1 + cr / chb)

    if -MINIMUM_VOLTAGE < vhb0 < MINIMUM_VOLTAGE:  # low-side turning off
        # output cannot turn on when low-side is turning off
        # this is guaranteed by:
        #   1) low-side turning off only happens in l0, and
        #   2) when low-side is turned off, the resonant voltage v reduces
        #      while the half-bridge voltage vhb increases, so that the voltage across
        #      the magnetizing inductor cannot rise to v_diode_on
        v = ckt.vbus  # target of resonant voltage v
        v -= -(vout - v0) / (1 + chb / cr) + vhb0 / (1 + cr / chb)
        v *= (1 + chb / cr)
        if -r <= v <= r:
            dt = min((math.asin(v / -r) - phi) % (2 * math.pi),
                     (math.pi - math.asin(v / -r) - phi) % (2 * math.pi)) / w
        else:  # hard-switching, which should be avoided but is still allowed
            dt = (-math.pi / 2 - phi) % (2 * math.pi) / w
        next_state = state_h0
    else:  # high-side turning off
        v = 0  # target of resonant voltage v
        v -= -(vout - v0) / (1 + chb / cr) + vhb0 / (1 + cr / chb)
        v *= (1 + chb / cr)
        if -r <= v <= r:
            dtcomm = min((math.asin(v / -r) - phi) % (2 * math.pi),
                         (math.pi - math.asin(v / -r) - phi) % (2 * math.pi)) / w
        else:  # hard-switching, which is allowed
            dtcomm = (math.pi / 2 - phi) % (2 * math.pi) / w

        v = ckt.vout / ckt.lm * (ckt.lr + ckt.lm)
        v -= vcen
        if -r <= v <= r:
            dtdion = min((math.asin(v / r) - phi) % (2 * math.pi),
                         (math.pi - math.asin(v / r) - phi) % (2 * math.pi)) / w
        else:
            dtdion = math.inf

        if dtdion < dtcomm:  # diode turns on before commutation finishes
            dt = dtdion
            next_state = state_c1
        else:
            dt = dtcomm
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

    v = 0  # target of the resonant voltage v
    v -= -(vout - v0) / (1 + chb / cr) + vhb0 / (1 + cr / chb)
    v *= (1 + chb / cr)
    if -r <= v <= r:
        dtcomm = min((math.asin(v / -r) - phi) % (2 * math.pi),
                     (math.pi - math.asin(v / -r) - phi) % (2 * math.pi)) / w
    else:  # hard-switching, which is allowed
        dtcomm = (math.pi / 2 - phi) % (2 * math.pi) / w

    ta = tb = (-phi % (2 * math.pi)) / w
    while r * math.cos(w * ta + phi) / z > im0 - km * ta and ta > MINIMUM_TIME:
        ta /= 2
    if ta > MINIMUM_TIME:
        dtdiof = nsolve(lambda t: im0 - km * t - r * math.cos(w * t + phi) / z, ta, tb)
    else:
        dtdiof = ta

    if dtdiof < dtcomm:  # diode turns off before commutation finishes (vhb reduces to 0), almost impossible
        dt = dtdiof
        next_state = state_c0
    else:
        dt = dtcomm
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


def _sim_phase_simple(isf, cond, i0, v0, vhb0, im0, ckt, con):
    nsf, s = isf(i0, v0, vhb0, im0, ckt, con)
    states = [s]
    while cond(nsf):
        nsf, s = nsf(s.i1, s.v1, s.vhb1, s.im1, ckt, con)
        states.append(s)
    return nsf, states


def _sim_phase(isf, cond, i0, v0, vhb0, im0, ckt, con):
    is_ls_on_phase = isf in {state_l0, state_l1} and con is not None
    if is_ls_on_phase:
        t12given = con
        t12hist = []

    nsf, s = isf(i0, v0, vhb0, im0, ckt, con)

    states = [s]
    while cond(nsf):
        csf = nsf
        nsf, s = csf(s.i1, s.v1, s.vhb1, s.im1, ckt, con)

        if is_ls_on_phase:
            if csf is state_l0:
                t12hist.append(s.dt)
            else:
                if s.dt < MINIMUM_DETECTABLE_DIODE_ON_TIME:  # l1 duration too short
                    con = t12given - sum(t12hist)
                else:
                    t12hist = []
                    con = t12given
        states.append(s)

    return nsf, states


def sim(v0, ckt, con):
    commutating = lambda s: s in {state_c0, state_c1}
    active = lambda s: s not in {state_c0, state_c1}

    def max_dv(i0, v0):
        l = ckt.lr + ckt.lm
        c = 1 / (1 / ckt.cr + 1 / ckt.chb)
        z = (l / c)**.5
        r = math.hypot(v0, i0 * z)
        phi = math.atan2(v0, i0 * z)
        return r * (1 + math.sin(phi))

    def eq_zvon(t):
        _, ls_on = _sim_phase(ls_on_inf, active, hs_off[-1].i1, hs_off[-1].v1, hs_off[-1].vhb1, hs_off[-1].im1, ckt, t)
        return max_dv(ls_on[-1].i1, ls_on[-1].v1) - (1 + ckt.chb / ckt.cr) * ckt.vbus

    conv, t12min = con
    i0 = 0
    vhb0 = ckt.vbus

    # high-side on phase
    nsf, hs_on = _sim_phase(state_h0, active, i0, v0, vhb0, i0, ckt, conv)

    # high-side off phase
    ls_on_inf, hs_off = _sim_phase(nsf, commutating, hs_on[-1].i1, hs_on[-1].v1, hs_on[-1].vhb1, hs_on[-1].im1, ckt, None)

    # low-side on phase
    _, ls_on = _sim_phase(ls_on_inf, active, hs_off[-1].i1, hs_off[-1].v1, hs_off[-1].vhb1, hs_off[-1].im1, ckt, t12min)
    if max_dv(ls_on[-1].i1, ls_on[-1].v1) < (1 + ckt.chb / ckt.cr) * ckt.vbus:
        t12max = (math.pi - ls_on[-1].phi) / ls_on[-1].w
        if t12min < t12max <= MAXIMUM_T12:
            if 0 < eq_zvon(t12max):
                t12_zvon = nsolve(eq_zvon, t12min, t12max)
            else:
                t12_zvon = t12max
            _, ls_on = _sim_phase(ls_on_inf, active, hs_off[-1].i1, hs_off[-1].v1, hs_off[-1].vhb1, hs_off[-1].im1, ckt, t12_zvon)

    # low-side off phase
    _, ls_off = _sim_phase(state_c0, commutating, ls_on[-1].i1, ls_on[-1].v1, ls_on[-1].vhb1, ls_on[-1].im1, ckt, None)

    # high-side on phase
    _, hs_on2 = _sim_phase(state_h0, active, ls_off[-1].i1, ls_off[-1].v1, ls_off[-1].vhb1, ls_off[-1].im1, ckt, None)

    return hs_on2[-1].v1 - v0, hs_on + hs_off + ls_on + ls_off + hs_on2


def evaluate_switching_period(states):
    if (abs(states[0].i0 - states[-1].i1) < MINIMUM_CURRENT and
        abs(states[0].v0 - states[-1].v1) < MINIMUM_VOLTAGE):
        operation = 'steady-state'
    else:
        operation = 'dynamics'

    tsw = 0
    qout = i2acc = i2acc_out = imax = imax_out = vacc = v2acc = 0
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
            km = state.km
            im0 = state.im0
            phi1 = state.w * state.dt + state.phi
            i2acc_out += j**2 * state.dt / 2
            i2acc_out += j**2 / (4 * state.w) * (math.sin(2 * phi1) - math.sin(2 * state.phi))
            i2acc_out -= 2 * j * im0 / state.w * (math.sin(phi1) - math.sin(state.phi))
            i2acc_out -= 2 * j * km / state.w**2 * (math.cos(phi1) - math.cos(state.phi) + state.w * state.dt * math.sin(phi1))
            i2acc_out += im0**2 * state.dt + km * state.dt**2 + km**2 * state.dt**3 / 3
            th0 = math.pi - math.asin(-km / (j * state.w))  # iout max's out at math.pi - theta, where theta is a small angle
            if state.phi <= th0 <= phi1:
                if im0 + km * (th0 - state.phi) / state.w - j * math.cos(th0) > imax_out:
                    imax_out = im0 + km * (th0 - state.phi) / state.w - j * math.cos(th0)
            else:
                if im0 + km * (phi1 - state.phi) / state.w - j * math.cos(phi1) > imax_out:
                    imax_out = im0 + km * (phi1 - state.phi) / state.w - j * math.cos(phi1)
        if state.state.startswith('c'):
            if state.vhb0 > state.vhb1:
                t_hs_off += state.dt
            else:
                t_ls_off += state.dt
        if state.state.startswith('l'):
            t_ls_on += state.dt
        if state.state.startswith('h'):
            t_hs_on += state.dt
        if state.phi <= 0 <= state.w * state.dt + state.phi:
            imax = j

    irms = (i2acc / tsw)**.5
    irms_out = (i2acc_out / tsw)**.5
    vavg = vacc / tsw
    vrms = (v2acc / tsw)**.5
    iout = qout / tsw

    t12 = [s for s in states if s.state == 'l0'][-1].dt
    return Evaluation(operation=operation,
                      tsw=tsw,
                      irms=irms,
                      irms_out=irms_out,
                      imax=imax,
                      imax_out=imax_out,
                      vavg=vavg,
                      vrms=vrms,
                      iout=iout,
                      t_diode_on=t_diode_on,
                      t_hs_off=t_hs_off,
                      t_ls_on=t_ls_on,
                      t_ls_off=t_ls_off,
                      t_hs_on=t_hs_on,
                      t12=t12)


def find_steady_state(tfwd, ckt, t12min=500e-9, fswmax=100e3):
    # print(f'finding steady-state for tfwd={fmt(tfwd)}s...')
    def eq(v0):
        # print(f'  trying from v0={fmt(v0, 4)}V')
        dv, ss = sim(v0, ckt, (tfwd, t12min))
        fsw = 1 / sum(s.dt for s in ss)
        if fsw > fswmax:
            def eqf(t12):
                _, ss = sim(v0, ckt, (tfwd, t12))
                fsw = 1 / sum(s.dt for s in ss)
                return fsw - fswmax
            t12 = nsolve(eqf, t12min, 1 / fswmax)
            dv, ss = sim(v0, ckt, (tfwd, t12))
        return dv, ss

    v0max = ckt.vbus - 1
    v0min = -ckt.vout / ckt.lm * (ckt.lr + ckt.lm)
    v0 = nsolve(lambda v: eq(v)[0], v0min * .99 + v0max * .01, v0max * .99 + v0min * .01)
    residue, ss = eq(v0)

    if abs(residue) > MINIMUM_VOLTAGE:
        print(f'Warning: cannot find steady-state for tfwd = {fmt(tfwd)}s.')
    return ss


def evaluate_operating_point(pout, ckt, t12min=500e-9, fswmax=100e3):
    # tfwdmax = math.pi / 2 / ((ckt.lr + ckt.lm) * ckt.cr)**-.5
    vdion = ckt.vout / ckt.lm * (ckt.lr + ckt.lm)
    tfwdmax = (math.pi / 2 + math.asin(vdion / (ckt.vbus + vdion))) / ((ckt.lr + ckt.lm) * ckt.cr)**-.5
    sspmax = find_steady_state(tfwdmax, ckt, t12min, fswmax)
    pmax = evaluate_switching_period(sspmax).iout * ckt.vout
    if 0 < pout <= pmax:
        tf = nsolve(lambda t: evaluate_switching_period(find_steady_state(t, ckt, t12min, fswmax)).iout * ckt.vout - pout,
                    MINIMUM_FORWARD_TIME, tfwdmax)
        ss = find_steady_state(tf, ckt, t12min, fswmax)
        return tf, ss, evaluate_switching_period(ss), pmax
    else:
        return 0, [State()], Evaluation(), pmax
