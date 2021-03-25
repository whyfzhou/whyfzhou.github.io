import math
from scipy.optimize import ridder as nsolve
import ahbllc
import plothelper


def evaluate_conditions(conditions, lr, lm, cr, use_tls_after=True):
    vins = []
    vouts = []
    pouts = []
    tlss = []
    fsws = []
    duties = []
    dduties = []
    iouts = []
    irms_outs = []
    ipk_outs = []
    irms_pris = []
    ippk_pris = []
    inpk_pris = []
    vrms_rs = []
    vrms_pris = []
    psi_pris = []
    qcomms = []
    thsrevs = []
    thsfwds = []
    tlsrevs = []
    tlsfwds = []
    tls_diodes = []
    tls_befores = []
    tls_afters = []
    for vin, vout, pout, t in conditions:
        if not use_tls_after:
            tls = t
        else:
            def eq(tl):
                *_, ta = evaluate_design(vin, vout, pout, tl, lr, lm, cr)
                return ta - t
            tls = nsolve(eq, 1.1 * math.pi * (lr * cr)**.5, 2.1 * math.pi * (lr * cr)**.5)
        perf, thsrev, thsfwd, tlsrev, tlsfwd, tls_diode, tls_before, tls_after = evaluate_design(vin, vout, pout, tls, lr, lm, cr)
        vins.append(vin)
        vouts.append(vout)
        pouts.append(pout)
        tlss.append(tls)
        fsws.append(perf['fsw'])
        duties.append(perf['duty'])
        dduties.append(perf['dduty'])
        iouts.append(perf['iout'])
        irms_outs.append(perf['irms_out'])
        ipk_outs.append(perf['ipk_out'])
        irms_pris.append(perf['irms_pri'])
        ippk_pris.append(perf['ippk_pri'])
        inpk_pris.append(perf['inpk_pri'])
        vrms_rs.append(perf['vrms_r'])
        vrms_pris.append(perf['vrms_pri'])
        psi_pris.append(perf['psi_pri'])
        qcomms.append(perf['qcomm'])
        thsrevs.append(thsrev)
        thsfwds.append(thsfwd)
        tlsrevs.append(tlsrev)
        tlsfwds.append(tlsfwd)
        tls_diodes.append(tls_diode)
        tls_befores.append(tls_before)
        tls_afters.append(tls_after)

    s = []
    s.append('{:>10} = {}H'.format('Lr', plothelper.fmt(lr)))
    s.append('{:>10} = {}H'.format('Lm', plothelper.fmt(lm)))
    s.append('{:>10} = {}F'.format('Cr', plothelper.fmt(cr)))
    lmlr = lm / lr
    s.append('{:>10} = {}'.format('1/λ', plothelper.fmt(lmlr)))
    s.append('{:>10} = {}rad/s'.format('ω0', plothelper.fmt((lr * cr)**-.5)))
    s.append('{:>10} = {}Hz'.format('f0', plothelper.fmt((lr * cr)**-.5 / 2 / math.pi)))
    s.append('{:>10} = {}Ω'.format('Z0', plothelper.fmt((lr / cr)**.5)))
    s.append(','.join(['vin', *[plothelper.fmt(x) for x in vins]]))
    s.append(','.join(['vout', *[plothelper.fmt(x) for x in vouts]]))
    s.append(','.join(['pout', *[plothelper.fmt(x) for x in pouts]]))
    s.append(','.join(['tls', *[plothelper.fmt(x) for x in tlss]]))
    s.append(','.join(['fsw', *[plothelper.fmt(x) for x in fsws]]))
    s.append(','.join(['duty', *[plothelper.fmt(x) for x in duties]]))
    s.append(','.join(['dduty', *[plothelper.fmt(x) for x in dduties]]))
    s.append(','.join(['iout', *[plothelper.fmt(x) for x in iouts]]))
    s.append(','.join(['irms_out', *[plothelper.fmt(x) for x in irms_outs]]))
    s.append(','.join(['ipk_out', *[plothelper.fmt(x) for x in ipk_outs]]))
    s.append(','.join(['irms_pri', *[plothelper.fmt(x) for x in irms_pris]]))
    s.append(','.join(['ippk_pri', *[plothelper.fmt(x) for x in ippk_pris]]))
    s.append(','.join(['inpk_pri', *[plothelper.fmt(x) for x in inpk_pris]]))
    s.append(','.join(['vrms_r', *[plothelper.fmt(x) for x in vrms_rs]]))
    s.append(','.join(['vrms_pri', *[plothelper.fmt(x) for x in vrms_pris]]))
    s.append(','.join(['psi_pri', *[plothelper.fmt(x) for x in psi_pris]]))
    s.append(','.join(['qcomm', *[plothelper.fmt(x) for x in qcomms]]))
    s.append(','.join(['thsrev', *[plothelper.fmt(x) for x in thsrevs]]))
    s.append(','.join(['thsfwd', *[plothelper.fmt(x) for x in thsfwds]]))
    s.append(','.join(['tlsrev', *[plothelper.fmt(x) for x in tlsrevs]]))
    s.append(','.join(['tlsfwd', *[plothelper.fmt(x) for x in tlsfwds]]))
    s.append(','.join(['tls_diode', *[plothelper.fmt(x) for x in tls_diodes]]))
    s.append(','.join(['tls_before', *[plothelper.fmt(x) for x in tls_befores]]))
    s.append(','.join(['tls_after', *[plothelper.fmt(x) for x in tls_afters]]))

    return '\n'.join(s)


def evaluate_design(vin, vout, pout, tls, lr, lm, cr, show=False):
    ckt = asyllcm1.AsymmetricalHalfBridgeLLC(lr, lm, cr, vin, vout)
    i0, *_ = asyllcm1.find_i0(pout, tls, ckt)
    ss = asyllcm1.find_ss(i0, tls, ckt)
    sampled_data = asyllcm1.sample(ss)
    perf = asyllcm1.evaluate_sampled(*sampled_data, ss)

    lmlr = lm / lr
    thsrev, thsfwd, tlsrev, tlsfwd, tls_diode, tls_before, tls_after = asyllcm1.calc_dts(*sampled_data, ss)
    ths = 1 / perf['fsw'] - tls

    if show:
        print('## Evaluating Design ##')
        print('{:>10} = {}H'.format('Lr', plothelper.fmt(lr)))
        print('{:>10} = {}H'.format('Lm', plothelper.fmt(lm)))
        print('{:>10} = {}F'.format('Cr', plothelper.fmt(cr)))
        print('{:>10} = {}'.format('1/λ', plothelper.fmt(lmlr)))
        print('{:>10} = {}rad/s'.format('ω0', plothelper.fmt((lr * cr)**-.5)))
        print('{:>10} = {}Hz'.format('f0', plothelper.fmt((lr * cr)**-.5 / 2 / math.pi)))
        print('{:>10} = {}Ω'.format('Z0', plothelper.fmt((lr / cr)**.5)))
        print('## At Oprating Point ##')
        print('{:>10} = {}V'.format('Vin', plothelper.fmt(vin)))
        print('{:>10} = {}V'.format('Vout', plothelper.fmt(vout)))
        print('{:>10} = {}W'.format('Pout', plothelper.fmt(pout)))
        print('## With Control Constants ##')
        print('{:>10} = {}s'.format('Tls', plothelper.fmt(tls)))
        print('## Results ##')
        print('{:>10} = {}Hz'.format('fsw', plothelper.fmt(perf['fsw'])))
        print('{:>10} = {:.2%}'.format('δhs', perf['duty']))
        print('{:>10} = {:.2%}'.format('δd,ls', perf['dduty']))
        print('{:>10} = {}A'.format('Iout', plothelper.fmt(perf['iout'])))
        print('{:>10} = {}A'.format('Iout,rms', plothelper.fmt(perf['irms_out'])))
        print('{:>10} = {}A'.format('Iout,pk', plothelper.fmt(perf['ipk_out'])))
        print('{:>10} = {}A'.format('Ipri,rms', plothelper.fmt(perf['irms_pri'])))
        print('{:>10} = {}A'.format('Ipri,ppk', plothelper.fmt(perf['ippk_pri'])))
        print('{:>10} = {}A'.format('Ipri,npk', plothelper.fmt(perf['inpk_pri'])))
        print('{:>10} = {}V'.format('Vr,rms', plothelper.fmt(perf['vrms_r'])))
        print('{:>10} = {}V'.format('Vpri,rms', plothelper.fmt(perf['vrms_pri'])))
        print('{:>10} = {}Wb'.format('Ψpri', plothelper.fmt(perf['psi_pri'])))
        print('{:>10} = {}C'.format('Qcomm', plothelper.fmt(perf['qcomm'])))
        print('## Details ##')
        print('{:>10} = {}s'.format('Ths', plothelper.fmt(ths)))
        print('{:>10} = {}s'.format('Ths,rev', plothelper.fmt(thsrev)))
        print('{:>10} = {}s'.format('Ths,fwd', plothelper.fmt(thsfwd)))
        print('{:>10} = {}s'.format('Tls', plothelper.fmt(tls)))
        print('{:>10} = {}s'.format('Tls,rev', plothelper.fmt(tlsrev)))
        print('{:>10} = {}s'.format('Tls,fwd', plothelper.fmt(tlsfwd)))
        print('{:>10} = {}s'.format('Tls,before', plothelper.fmt(tls_before)))
        print('{:>10} = {}s'.format('Tls,diode', plothelper.fmt(tls_diode)))
        print('{:>10} = {}s'.format('Tls,after', plothelper.fmt(tls_after)))


    return perf, thsrev, thsfwd, tlsrev, tlsfwd, tls_diode, tls_before, tls_after


def plot_design(vin, vout, pout, tls, lr, lm, cr, show=True, filename=None):
    ckt = asyllcm1.AsymmetricalHalfBridgeLLC(lr, lm, cr, vin, vout)
    i0, *_ = asyllcm1.find_i0(pout, tls, ckt)
    ss = asyllcm1.find_ss(i0, tls, ckt)
    sampled_data = asyllcm1.sample(ss)
    perf = asyllcm1.evaluate_sampled(*sampled_data, ss)
    plothelper.plot(ss, show=show, filename=filename)
    plothelper.plot_sampled(*sampled_data, ss, perf, tls, ckt, show=show, filename=filename)
