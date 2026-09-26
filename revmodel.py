# -*- coding: utf-8 -*-
"""
Revised thermodynamic model of the data-center cooling system (revision of NEXUS-D-26-00499).

Changes with respect to the submitted model (cooling_tower_tradeoff.py / cooling_tower_vectorized_corr.py):
  1. Moist-air entropy: dry air and water vapor are evaluated at their partial pressures
     (ideal-gas mixture), s = s_a(T) - R_a ln(p_a/p0) + w [s_v(T) - R_v ln(p_v/p_sat(T))]   (Reviewer 5.1)
  2. Chiller mode: the heat rejected by the tower includes the compressor work,
     Q_T = P_ITE (1 + 1/COP)                                                                 (Reviewer 5.2)
  3. Outlet state of the tower: counterflow closure based on the local second law (Merkel operating
     line kept below the saturation enthalpy of the water, shifted by the design approach a_ev). The air
     flow is the minimum compatible with that constraint, and the outlet air is saturated. The global
     entropy balance is kept as a check (S_gen >= 0).                                        (Reviewer 5.3, 7.4)
  4. Feasibility of the modes: dry free cooling if T2 <= T4 - a_dry; evaporative cooling if the inlet
     wet-bulb temperature Twb2 <= T4 - a_ev; otherwise chiller plus evaporative cooling with the design
     condenser-water set point T4 = 26 C, raised to Twb2 + a_ev when needed (condenser-water reset),
     with T1 = T4 + 24 C.                                                                    (Reviewers 2.4, 5.5, 7.3)
  5. Atmospheric pressure of each station from the standard atmosphere at its elevation.   (Reviewer 7.10)
  6. Auxiliary power: P = dp * volumetric flow with effective design pressure rises (efficiency folded in),
     250 Pa on the air stream and 0.21 MPa on the water loop.                                 (Reviewer 5.4)
Property fits (Table 1) are the REFPROP polynomials of the submitted model, whose independent variable is
X = T[C] + 273.16 (0.01 K above the thermodynamic temperature; the shift is absorbed by the fitted coefficients).
The same offset K is used in the ideal-gas volume and in the Carnot expression of the lift-based COP (effect < 0.004 %).
Optional sensitivity case: tower_lam = water value (kWh/L) for a water-saving fan control that raises the tower air
flow above the minimum when this lowers fan energy + tower_lam * evaporation.
All quantities are per kW of ITE power.
"""
import numpy as np

K = 273.16
R_U = 8.314            # kJ/(kmol K)
MW_W, MW_A = 18.015, 28.97
R_A, R_V = R_U/MW_A, R_U/MW_W      # kJ/(kg K)
EPS = MW_W/MW_A
P0 = 1.0e5             # Pa, pressure of the property fits

# Table 1 coefficients (a0, a1, a2); T in K
C_hL = [-1143.79850778, 4.18765442]
C_hV = [1.89674720e+03, 2.54423376e+00, -1.22199599e-03]
C_hA = [124.0129842, 1.00763588]
C_sL = [-5.59476137e+00, 2.59636662e-02, -2.00282656e-05]
C_sV = [2.13062766e+01, -6.41403034e-02, 7.17739861e-05]
C_sA = [2.45260622e+00, 6.22899746e-03, -4.78974682e-06]

DEFAULT = dict(dp_fan=250.0, dp_pump=0.21e6, a_ev=4.0, a_dry=6.0, T4_design=26.0, dT_chiller=24.0,
               T_ite_supply=18.0, T_ite_return=27.0, cop_a0=6.88, cop_a1=0.18, cop_a2=0.12, cop_model='eq15',
               cop_scale=1.0, f_cpL=1.0, f_hfg=1.0, f_cpa=1.0, RH5=100.0, comp_heat=True, lam=np.inf, ngrid=81,
               tower_lam=None)


def poly(c, T):
    Tk = np.asarray(T, float) + K
    return sum(ci*Tk**i for i, ci in enumerate(c))


class Props:
    """Property functions with optional reference-consistent perturbation factors (pinned at 20 C)."""
    def __init__(self, f_cpL=1.0, f_hfg=1.0, f_cpa=1.0):
        self.f_cpL, self.f_hfg, self.f_cpa = f_cpL, f_hfg, f_cpa
        self.hL20, self.hA20 = poly(C_hL, 20.0), poly(C_hA, 20.0)

    def hL(self, T):
        return self.hL20 + self.f_cpL*(poly(C_hL, T) - self.hL20)

    def hV(self, T):
        return self.hL(T) + self.f_hfg*(poly(C_hV, T) - poly(C_hL, T))

    def hA(self, T):
        return self.hA20 + self.f_cpa*(poly(C_hA, T) - self.hA20)

    @staticmethod
    def sL(T):
        return poly(C_sL, T)

    @staticmethod
    def sV(T):
        return poly(C_sV, T)

    @staticmethod
    def sA(T):
        return poly(C_sA, T)


def psat(T):
    """Eq. (2), Pa."""
    T = np.asarray(T, float)
    return 610.78*np.exp(17.27*T/(237.3+T))


def w_from(pv, p):
    return EPS*pv/(p-pv)


def p_std(z):
    """Standard-atmosphere pressure (Pa) at elevation z (m)."""
    return 101325.0*(1.0-2.25577e-5*np.asarray(z, float))**5.25588


def cop_fun(T4, T1, prm):
    if prm['cop_model'] == 'eq15':
        c = prm['cop_a0'] + prm['cop_a1']*prm['T_ite_supply'] - prm['cop_a2']*T4
    elif prm['cop_model'] == 'lift':
        # Carnot fraction on the refrigerant lift: evaporator 2 K below the chilled-water outlet,
        # condenser 1.5 K above the leaving condenser water (T1); fraction eta calibrated so that the
        # lift model equals Eq. (16) at a standard rating point (7 C chilled water, 32 C in, 37 C out).
        Tev, Tcd = prm['T_ite_supply']-2.0+K, T1+1.5+K
        eta = LIFT_ETA
        c = eta*Tev/(Tcd-Tev)
    return prm['cop_scale']*c


def _lift_eta():
    Tchw, Tin, Tout = 7.0, 32.0, 37.0   # standard rating conditions of the source chiller data [Kim & Kim 2025]
    cop_ref = 6.88 + 0.18*Tchw - 0.12*Tin
    Tev, Tcd = Tchw-2.0+K, Tout+1.5+K
    return cop_ref/(Tev/(Tcd-Tev))


LIFT_ETA = _lift_eta()


def wet_bulb(T2, RH2, p, pr):
    T2 = np.asarray(T2, float)
    pv = RH2/100*psat(T2); w2 = w_from(pv, p); H2 = pr.hA(T2) + w2*pr.hV(T2)
    lo = np.full_like(T2, -45.0); hi = T2.copy()
    for _ in range(48):
        m = 0.5*(lo+hi); ws = w_from(psat(m), p)
        f = pr.hA(m) + ws*pr.hV(m) - (ws-w2)*pr.hL(m) - H2
        hi = np.where(f > 0, m, hi); lo = np.where(f > 0, lo, m)
    return 0.5*(lo+hi)


def entropy_gen(mw, T1, T4, T2, w2, pv2, T5, w5, pv5, ma, T3, p, pr):
    """Global entropy generation, kW/K per kW_ITE, with the ideal-gas mixture at partial pressures."""
    s2 = pr.sA(T2) - R_A*np.log((p-pv2)/P0) + w2*(pr.sV(T2) - R_V*np.log(pv2/psat(T2)))
    s5 = pr.sA(T5) - R_A*np.log((p-pv5)/P0) + w5*(pr.sV(T5) - R_V*np.log(pv5/psat(T5)))
    return mw*(pr.sL(T4)-pr.sL(T1)) + ma*(s5-s2) - (w5-w2)*ma*pr.sL(T3)


def tower_states(T1, T4, T2, RH2, Q, p, prm, pr, mode_evap):
    """Tower (evaporative) or dry-cooler (mode_evap False) outlet state for given water temperatures.
    Returns dict of arrays: ma, T5, w5, WUE, VA5, mw, S."""
    T1, T4, T2, RH2, Q = [np.asarray(x, float) for x in (T1, T4, T2, RH2, Q)]
    p = np.broadcast_to(np.asarray(p, float), T2.shape).astype(float)
    pv2 = RH2/100*psat(T2); w2 = w_from(pv2, p); H2 = pr.hA(T2) + w2*pr.hV(T2)
    hL1, hL4, hL3 = pr.hL(T1), pr.hL(T4), pr.hL(T2)
    mw = Q/(hL1-hL4)
    n = T2.shape[0]
    out = dict(ma=np.full(n, np.nan), T5=np.full(n, np.nan), w5=np.full(n, np.nan), WUE=np.full(n, np.nan),
               VA5=np.full(n, np.nan), mw=mw, S=np.full(n, np.nan))
    if not mode_evap:
        T5 = T1 - prm['a_dry']
        ma = Q/(pr.hA(T5)+w2*pr.hV(T5) - H2)
        VA5 = ma*R_A*1e3*(K+T5)/(p-pv2)
        out.update(ma=ma, T5=T5, w5=w2, WUE=0.0*ma, VA5=VA5,
                   S=entropy_gen(mw, T1, T4, T2, w2, pv2, T5, w2, pv2, ma, T2, p, pr))
        return out
    a = prm['a_ev']
    f = np.linspace(0.0, 1.0, prm['ngrid'])[1:]
    for s in range(0, n, 20000):
        sl = slice(s, min(n, s+20000))
        Tw = T4[sl, None] + f[None, :]*(T1[sl]-T4[sl])[:, None]
        Hs = pr.hA(Tw-a) + w_from(psat(Tw-a), p[sl, None])*pr.hV(Tw-a)
        need = (pr.hL(Tw)-hL4[sl, None])*mw[sl, None]
        ma = np.max(need/np.maximum(Hs-H2[sl, None], 1e-12), axis=1)
        # top of the tower: saturated outlet air not warmer than T1 - a (full balance incl. makeup enthalpy)
        Ttop = T1[sl]-a; wtop = w_from(psat(Ttop), p[sl])
        ma_top = Q[sl]/np.maximum(pr.hA(Ttop)+wtop*pr.hV(Ttop)-H2[sl]-(wtop-w2[sl])*hL3[sl], 1e-12)
        ma = np.maximum(ma, ma_top)
        # outlet state at RH5 (default saturated) from the energy balance, Eq. (7)
        rh5 = prm['RH5']/100.0
        lo = np.full(ma.shape, -20.0); hi = T1[sl] + 10.0
        for _ in range(55):
            m = 0.5*(lo+hi); pv5 = rh5*psat(m); w5 = w_from(pv5, p[sl])
            g = ma*(pr.hA(m)+w5*pr.hV(m) - H2[sl] - (w5-w2[sl])*hL3[sl]) - Q[sl]
            hi = np.where(g > 0, m, hi); lo = np.where(g > 0, lo, m)
        T5 = 0.5*(lo+hi); pv5 = rh5*psat(T5); w5 = w_from(pv5, p[sl])
        if prm.get('tower_lam') is not None:
            # water-saving fan control (sensitivity case): the air flow is raised above the minimum when this
            # lowers fan energy + tower_lam * evaporation; the outlet air stays at RH5 (energy balance, Eq. (7))
            best = (ma*R_A*1e3*(K+T5)/(p[sl]-pv5))*prm['dp_fan']/1000.0 + prm['tower_lam']*(w5-w2[sl])*ma*3600.0
            ma0 = ma.copy()
            for kf in (1.1, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0):
                mk = kf*ma0
                lo = np.full(mk.shape, -20.0); hi = T1[sl] + 10.0
                for _ in range(55):
                    m = 0.5*(lo+hi); pvk = rh5*psat(m); wk = w_from(pvk, p[sl])
                    g = mk*(pr.hA(m)+wk*pr.hV(m) - H2[sl] - (wk-w2[sl])*hL3[sl]) - Q[sl]
                    hi = np.where(g > 0, m, hi); lo = np.where(g > 0, lo, m)
                Tk = 0.5*(lo+hi); pvk = rh5*psat(Tk); wk = w_from(pvk, p[sl])
                Jk = (mk*R_A*1e3*(K+Tk)/(p[sl]-pvk))*prm['dp_fan']/1000.0 + prm['tower_lam']*(wk-w2[sl])*mk*3600.0
                better = Jk < best
                best = np.where(better, Jk, best); ma = np.where(better, mk, ma)
                T5 = np.where(better, Tk, T5); pv5 = np.where(better, pvk, pv5); w5 = np.where(better, wk, w5)
        out['ma'][sl] = ma; out['T5'][sl] = T5; out['w5'][sl] = w5
        out['WUE'][sl] = (w5-w2[sl])*ma*3600.0
        out['VA5'][sl] = ma*R_A*1e3*(K+T5)/(p[sl]-pv5)
        out['S'][sl] = entropy_gen(mw[sl], T1[sl], T4[sl], T2[sl], w2[sl], pv2[sl], T5, w5, pv5, ma, T2[sl], p[sl], pr)
    return out


def run(T2, RH2, p=P0, **kw):
    """Run the cooling-system model on a series of air states (T2 in C, RH2 in %).
    Mode selection minimizes J = PUE + lam*WUE among the feasible modes; lam = inf reproduces the
    water-priority rule (free cooling whenever feasible), lam = 0 is the energy-priority rule."""
    prm = dict(DEFAULT); prm.update(kw)
    pr = Props(prm['f_cpL'], prm['f_hfg'], prm['f_cpa'])
    T2 = np.asarray(T2, float); RH2 = np.clip(np.asarray(RH2, float), 0.5, 100.0)
    n = T2.shape[0]
    p = np.broadcast_to(np.asarray(p, float), T2.shape).astype(float)
    Ts, Tr = prm['T_ite_supply'], prm['T_ite_return']
    Twb = wet_bulb(T2, RH2, p, pr)
    free_ok = T2 <= Ts - prm['a_dry']
    evap_ok = Twb <= Ts - prm['a_ev']
    res = {}
    # free cooling
    one = np.ones(n)
    fr = tower_states(Tr*one, Ts*one, T2, RH2, one, p, prm, pr, False)
    ev = tower_states(Tr*one, Ts*one, T2, RH2, one, p, prm, pr, True)
    # chiller + evaporative, condenser-water reset
    T4c = np.maximum(prm['T4_design'], Twb + prm['a_ev'])
    T1c = T4c + prm['dT_chiller']
    cop = cop_fun(T4c, T1c, prm)
    Qc = 1.0 + (1.0/cop if prm['comp_heat'] else 0.0*cop)
    ch = tower_states(T1c, T4c, T2, RH2, Qc, p, prm, pr, True)

    def aux(st):
        return st['VA5']*prm['dp_fan']/1000.0, (st['mw']/1000.0)*prm['dp_pump']/1000.0
    ffan, fpump = aux(fr); efan, epump = aux(ev); cfan, cpump = aux(ch)
    PUE_f = 1 + ffan + fpump; PUE_e = 1 + efan + epump; PUE_c = 1 + cfan + cpump + 1.0/cop
    WUE_f = np.zeros(n); WUE_e = ev['WUE']; WUE_c = ch['WUE']
    lam = prm['lam']
    if np.isinf(lam):
        mode = np.where(free_ok, 0, np.where(evap_ok, 1, 2))
    else:
        Jf = np.where(free_ok, PUE_f + lam*WUE_f, np.inf)
        Je = np.where(evap_ok, PUE_e + lam*WUE_e, np.inf)
        Jc = PUE_c + lam*WUE_c
        mode = np.argmin(np.vstack([Jf, Je, Jc]), axis=0)
    sel = lambda a, b, c: np.choose(mode, [a, b, c])
    res['mode'] = mode
    res['reset'] = (mode == 2) & (T4c > prm['T4_design'] + 1e-9)
    res['PUE'] = sel(PUE_f, PUE_e, PUE_c)
    res['WUE'] = sel(WUE_f, WUE_e, WUE_c)
    res['fan'] = sel(ffan, efan, cfan); res['pump'] = sel(fpump, epump, cpump)
    res['chill'] = np.where(mode == 2, 1.0/cop, 0.0)
    res['T1'] = sel(Tr*one, Tr*one, T1c); res['T4'] = sel(Ts*one, Ts*one, T4c)
    res['T5'] = sel(fr['T5'], ev['T5'], ch['T5']); res['ma'] = sel(fr['ma'], ev['ma'], ch['ma'])
    res['VA5'] = sel(fr['VA5'], ev['VA5'], ch['VA5']); res['S'] = sel(fr['S'], ev['S'], ch['S'])
    res['Qtower'] = sel(one, one, Qc); res['cop'] = np.where(mode == 2, cop, np.nan)
    res['Twb'] = Twb
    return res
