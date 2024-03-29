import math
from scipy.optimize import brentq as nsolve
from plothelper import fmt


MINIMUM_TIME = 1e-9
MINIMUM_VOLTAGE = 1e-5
MINIMUM_CURRENT = 1e-7
MAXIMUM_COMMUTATION_TIME = 600e-9
MINIMUM_HIGH_SIDE_ON_TIME = 100e-9


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
        self.imax = kwargs.get('imax', 0)
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
            else:
                s += f'{"Idc,out":>28} = {fmt(self.iout, 3)}V\n'
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


def state_l1(i0, v0, vhb0, im0, ckt, con, state_func_list):
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
    next_state = state_func_list['state_l0']
    return next_state, State(state='l1', dt=dt,
                             i0=i0, v0=v0, vhb0=vhb0, im0=im0,
                             i1=i1, v1=v1, vhb1=vhb0, im1=i1,
                             r=r, phi=phi, w=w, z=z,
                             cr=cr, chb=chb,
                             vcen=vcen, vout=vout, km=-km) # yapf: disable


def state_l0(i0, v0, vhb0, im0, ckt, con, state_func_list):
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
        next_state = state_func_list['state_l1']
    elif i_t12 <= 0:
        dt = t12
        next_state = state_func_list['state_c0']
    else:
        dt = dt_izc
        next_state = state_func_list['state_c0']

    i1 = r * math.cos(w * dt + phi) / z
    v1 = r * math.sin(w * dt + phi) + vcen
    return next_state, State(state='l0', dt=dt,
                             i0=i0, v0=v0, vhb0=vhb0, im0=i0,
                             i1=i1, v1=v1, vhb1=vhb0, im1=i1,
                             r=r, phi=phi, w=w, z=z,
                             cr=cr, chb=chb,
                             vcen=vcen, vout=vout, km=0) # yapf: disable


def state_h0_voff(i0, v0, vhb0, im0, ckt, con, state_func_list):
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
    # voff == r * sin(w * t + phi) + vcen
    if -r <= voff - vcen <= r:
        dt = min((math.asin((voff - vcen) / r) - phi) % (2 * math.pi),
                 (math.pi - math.asin((voff - vcen) / r) - phi) % (2 * math.pi)) / w
    else:  # then let's switch off the high-side at maximum voltage
        dt = ((math.pi / 2 - phi) % (2 * math.pi)) / w
    i1 = r * math.cos(w * dt + phi) / z
    v1 = r * math.sin(w * dt + phi) + vcen

    next_state = state_func_list['state_c0']
    return next_state, State(state='h0', dt=dt,
                             i0=i0, v0=v0, vhb0=vhb0, im0=i0,
                             i1=i1, v1=v1, vhb1=vhb0, im1=i1,
                             r=r, phi=phi, w=w, z=z,
                             cr=cr, chb=chb,
                             vcen=vcen, vout=vout, km=0) # yapf: disable


def state_h0_cpm(i0, v0, vhb0, im0, ckt, con, state_func_list):
    """Sub-state of:
    h - high-side device on,
    0 - output rectifier off.

    Exit criteria: high-side turned off when v reaches voff

    cpm - current-programming-mode, or peak current control
    """
    del vhb0, im0
    ioff = con

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

    if i0 == 0:
        # v0 + dvoff == r * sin(w * t + phi) + vcen
        # ioff = r * cos(w * t + phi) / z
        assert -r <= ioff * z <= r
        dt = min((math.acos(ioff * z / r) - phi) % (2 * math.pi),
                 (-math.acos(ioff * z / r) - phi) % (2 * math.pi)) / w
    else:
        # 0 = r * cos(w * t + phi) / z, w * t + phi = -pi / 2
        dt = ((-math.pi / 2 - phi) % (2 * math.pi)) / w
    i1 = r * math.cos(w * dt + phi) / z
    v1 = r * math.sin(w * dt + phi) + vcen

    next_state = state_func_list['state_c0']
    return next_state, State(state='h0', dt=dt,
                             i0=i0, v0=v0, vhb0=vhb0, im0=i0,
                             i1=i1, v1=v1, vhb1=vhb0, im1=i1,
                             r=r, phi=phi, w=w, z=z,
                             cr=cr, chb=chb,
                             vcen=vcen, vout=vout, km=0) # yapf: disable


def state_h0_vinc(i0, v0, vhb0, im0, ckt, con, state_func_list):
    """Sub-state of:
    h - high-side device on,
    0 - output rectifier off.

    Exit criteria: high-side turned off when v reaches voff

    vinc - voltage increment control
    """
    del vhb0, im0
    dvoff = con

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

    if i0 == 0:
        # v0 + dvoff == r * sin(w * t + phi) + vcen
        assert -r <= v0 + dvoff - vcen <= r
        dt = min((math.asin((v0 + dvoff - vcen) / r) - phi) % (2 * math.pi),
                 (math.pi - math.asin((v0 + dvoff - vcen) / r) - phi) % (2 * math.pi)) / w
    else:
        # 0 = r * cos(w * t + phi) / z, w * t + phi = -pi / 2
        dt = ((-math.pi / 2 - phi) % (2 * math.pi)) / w
    i1 = r * math.cos(w * dt + phi) / z
    v1 = r * math.sin(w * dt + phi) + vcen

    next_state = state_func_list['state_c0']
    return next_state, State(state='h0', dt=dt,
                             i0=i0, v0=v0, vhb0=vhb0, im0=i0,
                             i1=i1, v1=v1, vhb1=vhb0, im1=i1,
                             r=r, phi=phi, w=w, z=z,
                             cr=cr, chb=chb,
                             vcen=vcen, vout=vout, km=0) # yapf: disable


def state_h0_tfwd(i0, v0, vhb0, im0, ckt, con, state_func_list):
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

    next_state = state_func_list['state_c0']
    return next_state, State(state='h0', dt=dt,
                             i0=i0, v0=v0, vhb0=vhb0, im0=i0,
                             i1=i1, v1=v1, vhb1=vhb0, im1=i1,
                             r=r, phi=phi, w=w, z=z,
                             cr=cr, chb=chb,
                             vcen=vcen, vout=vout, km=0) # yapf: disable


def state_h0_vincscomp(i0, v0, vhb0, im0, ckt, con, state_func_list):
    """Sub-state of:
    h - high-side device on,
    0 - output rectifier off.

    Exit criteria: high-side turned off when v reaches voff

    vincscomp - voltage increment control with slope compensation
    """
    del vhb0, im0
    dvoff = con

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

    slope = 40 / 10e-6
    if i0 == 0:
        # v0 + dvoff - slope * t == r * sin(w * t + phi) + vcen
        dt = nsolve(lambda t: v0 + dvoff - slope * t - r * math.sin(w * t + phi) - vcen, 0, math.pi / w)
        # dt = min((math.asin((v0 + dvoff - vcen) / r) - phi) % (2 * math.pi),
        #          (math.pi - math.asin((v0 + dvoff - vcen) / r) - phi) % (2 * math.pi)) / w
    else:
        # 0 = r * cos(w * t + phi) / z, w * t + phi = -pi / 2
        dt = ((-math.pi / 2 - phi) % (2 * math.pi)) / w
    i1 = r * math.cos(w * dt + phi) / z
    v1 = r * math.sin(w * dt + phi) + vcen

    next_state = state_func_list['state_c0']
    return next_state, State(state='h0', dt=dt,
                             i0=i0, v0=v0, vhb0=vhb0, im0=i0,
                             i1=i1, v1=v1, vhb1=vhb0, im1=i1,
                             r=r, phi=phi, w=w, z=z,
                             cr=cr, chb=chb,
                             vcen=vcen, vout=vout, km=0) # yapf: disable


def state_h0_dvoff(i0, v0, vhb0, im0, ckt, con, state_func_list):
    """Sub-state of:
    h - high-side device on,
    0 - output rectifier off.

    Exit criteria: high-side turned off when v reaches voff

    dvoff - delta-voltage turning-off control
    """
    del vhb0, im0
    dvoff = con

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

    assert -r <= v0 + dvoff - vcen <= r
    dt = min((math.asin((v0 + dvoff - vcen) / r) - phi) % (2 * math.pi),
                (math.pi - math.asin((v0 + dvoff - vcen) / r) - phi) % (2 * math.pi)) / w
    i1 = r * math.cos(w * dt + phi) / z
    v1 = r * math.sin(w * dt + phi) + vcen

    next_state = state_func_list['state_c0']
    return next_state, State(state='h0', dt=dt,
                             i0=i0, v0=v0, vhb0=vhb0, im0=i0,
                             i1=i1, v1=v1, vhb1=vhb0, im1=i1,
                             r=r, phi=phi, w=w, z=z,
                             cr=cr, chb=chb,
                             vcen=vcen, vout=vout, km=0) # yapf: disable


def state_c0(i0, v0, vhb0, im0, ckt, con, state_func_list):
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
        next_state = state_func_list['state_h0']
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
            next_state = state_func_list['state_c1']
        else:
            dt = dtcomm
            next_state = state_func_list['state_l0']

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


def state_c1(i0, v0, vhb0, im0, ckt, con, state_func_list):
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
        next_state = state_func_list['state_c0']
    else:
        dt = dtcomm
        next_state = state_func_list['state_l1']

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


def _sim_phase(isf, cond, i0, v0, vhb0, im0, ckt, con, sfl):
    nsf, s = isf(i0, v0, vhb0, im0, ckt, con, sfl)
    states = [s]
    while cond(nsf):
        nsf, s = nsf(s.i1, s.v1, s.vhb1, s.im1, ckt, con, sfl)
        states.append(s)
    return nsf, states


def sim(v0, ckt, con, constr=None):
    sfl = dict(state_l1=state_l1, state_l0=state_l0, state_h0=state_h0_tfwd, state_c0=state_c0, state_c1=state_c1)
    if constr is not None:
        sfl['state_h0'] = constr
    assert v0 < ckt.vbus
    commutating = lambda s: s in {sfl['state_c0'], sfl['state_c1']}
    active = lambda s: s not in {sfl['state_c0'], sfl['state_c1']}

    def max_dv(i0, v0):
        l = ckt.lr + ckt.lm
        c = 1 / (1 / ckt.cr + 1 / ckt.chb)
        z = (l / c)**.5
        r = math.hypot(v0, i0 * z)
        phi = math.atan2(v0, i0 * z)
        return r * (1 + math.sin(phi))

    def eq_zvon(t):
        _, ls_on = _sim_phase(ls_on_inf, active, hs_off[-1].i1, hs_off[-1].v1, hs_off[-1].vhb1, hs_off[-1].im1, ckt, t, sfl)
        return max_dv(ls_on[-1].i1, ls_on[-1].v1) - (1 + ckt.chb / ckt.cr) * ckt.vbus

    conv, t12min = con
    i0 = 0
    vhb0 = ckt.vbus

    # high-side on phase
    nsf, hs_on = _sim_phase(sfl['state_h0'], active, i0, v0, vhb0, i0, ckt, conv, sfl)

    # high-side off phase
    ls_on_inf, hs_off = _sim_phase(nsf, commutating, hs_on[-1].i1, hs_on[-1].v1, hs_on[-1].vhb1, hs_on[-1].im1, ckt, None, sfl)

    # low-side on phase
    _, ls_on = _sim_phase(ls_on_inf, active, hs_off[-1].i1, hs_off[-1].v1, hs_off[-1].vhb1, hs_off[-1].im1, ckt, t12min, sfl)
    if max_dv(ls_on[-1].i1, ls_on[-1].v1) < (1 + ckt.chb / ckt.cr) * ckt.vbus:
        t12max = (math.pi - ls_on[-1].phi) / ls_on[-1].w
        if t12min < t12max:
            if 0 < eq_zvon(t12max):
                t12_zvon = nsolve(eq_zvon, t12min, t12max)
            else:
                t12_zvon = t12max
            _, ls_on = _sim_phase(ls_on_inf, active, hs_off[-1].i1, hs_off[-1].v1, hs_off[-1].vhb1, hs_off[-1].im1, ckt, t12_zvon, sfl)

    # low-side off phase
    _, ls_off = _sim_phase(sfl['state_c0'], commutating, ls_on[-1].i1, ls_on[-1].v1, ls_on[-1].vhb1, ls_on[-1].im1, ckt, None, sfl)

    # high-side on phase
    _, hs_on2 = _sim_phase(sfl['state_h0'], active, ls_off[-1].i1, ls_off[-1].v1, ls_off[-1].vhb1, ls_off[-1].im1, ckt, None, sfl)

    return hs_on2[-1].v1 - v0, hs_on + hs_off + ls_on + ls_off + hs_on2


def sim_dvoff(i0, v0, ckt, con):
    """Simulator used only for delta-V turning-off control.
    The first return variable is a tuple of two floats
    representing the residue current and voltage.

    :param i0: initial current
    :type i0: float
    :param v0: initial voltage
    :type v0: float
    :param ckt: circuit
    :type ckt: AHBLLC
    :param con: controlling variable
    :type con: tuple of two floats
    :return: residues and states
    :rtype: Tuple[residues, states] where residues::Tuple[float, float], states::State
    """
    sfl = dict(state_l1=state_l1, state_l0=state_l0, state_h0=state_h0_dvoff, state_c0=state_c0, state_c1=state_c1)
    commutating = lambda s: s in {sfl['state_c0'], sfl['state_c1']}
    active = lambda s: s not in {sfl['state_c0'], sfl['state_c1']}

    def max_dv(i0, v0):
        l = ckt.lr + ckt.lm
        c = 1 / (1 / ckt.cr + 1 / ckt.chb)
        z = (l / c)**.5
        r = math.hypot(v0, i0 * z)
        phi = math.atan2(v0, i0 * z)
        return r * (1 + math.sin(phi))

    def eq_zvon(t):
        _, ls_on = _sim_phase(ls_on_inf, active, hs_off[-1].i1, hs_off[-1].v1, hs_off[-1].vhb1, hs_off[-1].im1, ckt, t, sfl)
        return max_dv(ls_on[-1].i1, ls_on[-1].v1) - (1 + ckt.chb / ckt.cr) * ckt.vbus

    dvoff, t12min = con
    vhb0 = ckt.vbus

    # high-side on phase
    nsf, hs_on = _sim_phase(state_h0_dvoff, active, i0, v0, vhb0, i0, ckt, dvoff, sfl)

    # high-side off phase
    ls_on_inf, hs_off = _sim_phase(nsf, commutating, hs_on[-1].i1, hs_on[-1].v1, hs_on[-1].vhb1, hs_on[-1].im1, ckt, None, sfl)

    # low-side on phase
    _, ls_on = _sim_phase(ls_on_inf, active, hs_off[-1].i1, hs_off[-1].v1, hs_off[-1].vhb1, hs_off[-1].im1, ckt, t12min, sfl)
    if max_dv(ls_on[-1].i1, ls_on[-1].v1) < (1 + ckt.chb / ckt.cr) * ckt.vbus:
        t12max = (math.pi - ls_on[-1].phi) / ls_on[-1].w
        if t12min < t12max:
            if 0 < eq_zvon(t12max):
                t12_zvon = nsolve(eq_zvon, t12min, t12max)
            else:
                t12_zvon = t12max
            _, ls_on = _sim_phase(ls_on_inf, active, hs_off[-1].i1, hs_off[-1].v1, hs_off[-1].vhb1, hs_off[-1].im1, ckt, t12_zvon, sfl)

    # low-side off phase
    _, ls_off = _sim_phase(state_c0, commutating, ls_on[-1].i1, ls_on[-1].v1, ls_on[-1].vhb1, ls_on[-1].im1, ckt, None, sfl)

    return (ls_off[-1].i1 - i0,  ls_off[-1].v1 - v0), hs_on + hs_off + ls_on + ls_off


def sim_voff(i0, ckt, con):
    sfl = dict(state_l1=state_l1, state_l0=state_l0, state_h0=state_h0_voff, state_c0=state_c0, state_c1=state_c1)
    commutating = lambda s: s in {sfl['state_c0'], sfl['state_c1']}
    active = lambda s: s not in {sfl['state_c0'], sfl['state_c1']}

    def eq_zvon(t, hs_off_fins):
        # from the end moment of the high-side turning-off state (high- to low-side commutation finishes),
        # find out what t12 gives ZV-on of the high-side device.
        _, ss = _sim_phase(isf_ls_on, active, *hs_off_fins, ckt, t, sfl)
        ls_on_fins = (ss[-1].i1, ss[-1].v1, ss[-1].vhb1, ss[-1].im1)
        _, ss = _sim_phase(state_c0, commutating, *ls_on_fins, ckt, t, sfl)
        c0 = ss[-1]
        r, chb, cr = c0.r, c0.chb, c0.cr
        vout, v0, vhb0 = c0['vout'], c0.v0, c0.vhb0
        vhbmax = r / (1 + chb/cr) - (vout-v0) / (1 + chb/cr) + vhb0 / (1 + cr/chb)
        return vhbmax - ckt.vbus

    # def eq_vcont(t, hs_off_fins, ratio_v_ls_on=.95):
    #     # resonant voltage continuity equation:
    #     # t should be so long that when we turn off the low-side,
    #     # v the resonant voltage is at (or below) ratio_v_ls_on * voff
    #     _, ss = _sim1(isf_ls_on, active, *hs_off_fins, ckt, t)
    #     v_ls_on_fin = ss[-1].v1
    #     return v_ls_on_fin - voff * ratio_v_ls_on

    # def eq_vcont(t, hs_off_fins, ratio_v_ls_on=.95):
    #     # resonant voltage continuity equation:
    #     # t should be so long that when we turn off the low-side and go into the high-side on state,
    #     # the minimum value of v, the resonant voltage, is at (or below) ratio_v_ls_on * voff
    #     # 这样做看起来很理想，但是
    #     #   1) 无法实现，除非控制器有办法得知 hs_on 状态中 v 的最小值
    #     #   2) 会改变 sim 返回值 i_fin - i0 随 i0 的单调性，使方程无法可靠求解
    #     _, ss = _sim1(isf_ls_on, active, *hs_off_fins, ckt, t)
    #     ls_on_fins = (ss[-1].i1, ss[-1].v1, ss[-1].vhb1, ss[-1].im1)
    #     _, ss = _sim1(state_c0, commutating, *ls_on_fins, ckt, t)
    #     v0 = ss[-1].v1
    #     i0 = ss[-1].i1
    #     z = ((ckt.lr + ckt.lm) / ckt.cr)**.5
    #     vcen = ckt.vbus
    #     r = math.hypot(v0 - vcen, i0 * z)
    #     vmin = -r + vcen
    #     return vmin - voff * ratio_v_ls_on

    voff, t12min = con
    # starting at the turning-off moment of the high-side device
    v0 = voff
    vhb0 = ckt.vbus

    # high-side turning-off phase
    nsf, hs_off = _sim_phase(state_c0, commutating, i0, v0, vhb0, i0, ckt, None, sfl)
    hs_off_fins = (hs_off[-1].i1, hs_off[-1].v1, hs_off[-1].vhb1, hs_off[-1].im1)

    isf_ls_on = nsf  # remember the initial state of the low-side on phase,
                     # we may needed repetitively in determine t12 for high-side ZV-on
    # low-side on phase
    nsf, ls_on = _sim_phase(isf_ls_on, active, *hs_off_fins, ckt, t12min, sfl)
    ls_on_fins = (ls_on[-1].i1, ls_on[-1].v1, ls_on[-1].vhb1, ls_on[-1].im1)
    t12capm = ls_on[-1].dt  # state_l0 automatically increases t12 and avoids capactive mode
    # low-side turning-off phase
    nsf, ls_off = _sim_phase(state_c0, commutating, *ls_on_fins, ckt, None, sfl)
    ls_off_fins = (ls_off[-1].i1, ls_off[-1].v1, ls_off[-1].vhb1, ls_off[-1].im1)

    # correct t12
    t12zvon = 0
    t12max = (math.pi - ls_on[-1].phi) / ls_on[-1].w  # the t12 that generates the most negative current before low-side turns off
    if ls_off[-1].vhb1 < ckt.vbus - MINIMUM_VOLTAGE:  # hard-switching occurs, see if we can avoid it
        if eq_zvon(t12max, hs_off_fins) >= 0:  # yes, we can
            t12zvon = nsolve(lambda t: eq_zvon(t, hs_off_fins), t12min / 2, t12max)
        else:  # hard-switching we cannot avoid, in conditions, e.g., too much chb
            t12zvon = t12max
    t12 = max(t12min, t12capm, t12zvon)
    # t12vcont = 0  # in practice, the switching frequency is limited
    # if ls_off[-1].v1 > voff:
    #     t12vcont = nsolve(lambda t: eq_vcont(t, hs_off_fins), t12min / 2, t12max)
    # t12 = max(t12, t12vcont)

    # t12 correction needed
    if t12 > t12min:  # re-run the low-side on and low- to high-side commutation
        nsf, ls_on = _sim_phase(isf_ls_on, active, *hs_off_fins, ckt, t12, sfl)
        ls_on_fins = (ls_on[-1].i1, ls_on[-1].v1, ls_on[-1].vhb1, ls_on[-1].im1)
        nsf, ls_off = _sim_phase(state_c0, commutating, *ls_on_fins, ckt, None, sfl)
        ls_off_fins = (ls_off[-1].i1, ls_off[-1].v1, ls_off[-1].vhb1, ls_off[-1].im1)

    # high-side on phase
    _, hs_on = _sim_phase(state_h0_voff, active, *ls_off_fins, ckt, voff, sfl)
    res = hs_on[-1].i1 - i0
    if hs_on[-1].v1 > voff + MINIMUM_VOLTAGE:
        dt_ls_off = sum(s.dt for s in ls_off)
        dt_hs_on = MAXIMUM_COMMUTATION_TIME - dt_ls_off
        if dt_hs_on < MINIMUM_HIGH_SIDE_ON_TIME:
            dt_hs_on = MINIMUM_HIGH_SIDE_ON_TIME
        hs_on_v0 = ls_off[-1].v1
        hs_on_i0 = ls_off[-1].i1
        hs_on_z = ((ckt.lr + ckt.lm) / ckt.cr)**.5
        hs_on_w = ((ckt.lr + ckt.lm) * ckt.cr)**-.5
        hs_on_vcen = ckt.vbus
        hs_on_r = math.hypot(hs_on_v0 - hs_on_vcen, hs_on_i0 * hs_on_z)
        hs_on_phi = math.atan2(hs_on_v0 - hs_on_vcen, hs_on_i0 * hs_on_z)
        hs_on_v1 = hs_on_r * math.sin(hs_on_w * dt_hs_on + hs_on_phi) + hs_on_vcen
        _, hs_on = _sim_phase(state_h0_voff, active, *ls_off_fins, ckt, hs_on_v1)
        # voltage continuity violates if we reach here, which indicates that
        # the i_fin - i0 function could have multiple zeros, and the solver may fail.
        # consider changing the res to keep it monotone
        res = hs_on[-1].i1 - i0

    states = hs_off + ls_on + ls_off + hs_on
    return res, t12, states  # TODO: 统一所有 simulator 的返回格式


def evaluate_switching_period(states):
    if (abs(states[0].i0 - states[-1].i1) < MINIMUM_CURRENT and
        abs(states[0].v0 - states[-1].v1) < MINIMUM_VOLTAGE):
        operation = 'steady-state'
    else:
        operation = 'dynamics'

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
        if state.phi <= 0 <= state.w * state.dt + state.phi:
            imax = j

    irms = (i2acc / tsw)**.5
    vavg = vacc / tsw
    vrms = (v2acc / tsw)**.5
    iout = qout / tsw

    t12 = [s for s in states if s.state == 'l0'][-1].dt
    return Evaluation(operation=operation,
                      tsw=tsw,
                      irms=irms,
                      imax=imax,
                      vavg=vavg,
                      vrms=vrms,
                      iout=iout,
                      t_diode_on=t_diode_on,
                      t_hs_off=t_hs_off,
                      t_ls_on=t_ls_on,
                      t_ls_off=t_ls_off,
                      t_hs_on=t_hs_on,
                      t12=t12)


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
        dv = ss[-1].v1 - voff
        return di, dv, ss

    i0max = ckt.vbus / ((ckt.lr + ckt.lm) / ckt.cr)**.5
    i0 = nsolve(lambda i: eq(i)[0], .001, i0max * 2)
    _, _, ss = eq(i0)
    return ss


def evaluate_operating_point(pout, ckt, t12min=500e-9, fswmax=100e3):
    pmax = evaluate_switching_period(find_steady_state(ckt.vbus, ckt, t12min, fswmax)).iout * ckt.vout
    if 0 < pout <= pmax:
        voff = nsolve(lambda v: evaluate_switching_period(find_steady_state(v, ckt, t12min, fswmax)).iout * ckt.vout - pout,
                      ckt.vout/10, ckt.vbus)
        ss = find_steady_state(voff, ckt, t12min, fswmax)
        return voff, ss, evaluate_switching_period(ss)
    else:
        return 0, [State()], Evaluation()
