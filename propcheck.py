# Checks of the property functions of revmodel.py (Section 2.1 and Table 2 of the article):
# property fits against IAPWS-95 and the reference equation of state of dry air (CoolProp), the saturation pressure,
# humidity ratio, and wet-bulb temperature against ASHRAE (psychrolib), and the moist-air entropy of Eq. (11) against
# the ASHRAE RP-1485 real moist-air formulation (CoolProp HumidAir).
import numpy as np
import CoolProp.CoolProp as CP
import psychrolib as PL
import revmodel as R

PL.SetUnitSystem(PL.SI)
pr = R.Props()

# ---------------- Table 2: fits of h and s, 0 to 100 C (0 to 99 C for liquid water), 0.1 MPa and saturated vapor
T = np.arange(0.0, 100.0 + 1e-9, 1.0)
TK = T + 273.15; TK[0] = 273.16
lq = T <= 99
ref = {'hL': np.array([CP.PropsSI('H', 'T', t, 'P', 1e5, 'Water') for t in TK[lq]])/1e3,
       'sL': np.array([CP.PropsSI('S', 'T', t, 'P', 1e5, 'Water') for t in TK[lq]])/1e3,
       'hv': np.array([CP.PropsSI('H', 'T', t, 'Q', 1, 'Water') for t in TK])/1e3,
       'sv': np.array([CP.PropsSI('S', 'T', t, 'Q', 1, 'Water') for t in TK])/1e3,
       'ha': np.array([CP.PropsSI('H', 'T', t, 'P', 1e5, 'Air') for t in TK])/1e3,
       'sa': np.array([CP.PropsSI('S', 'T', t, 'P', 1e5, 'Air') for t in TK])/1e3}
print('Property fits of Table 2 against CoolProp (IAPWS-95 water, reference equation of dry air):')
for lab, c, k, msk in (('hL', R.C_hL, 'hL', lq), ('hv', R.C_hV, 'hv', None), ('ha', R.C_hA, 'ha', None),
                       ('sL', R.C_sL, 'sL', lq), ('sv', R.C_sV, 'sv', None), ('sa', R.C_sA, 'sa', None)):
    Tt = T[msk] if msk is not None else T
    d = R.poly(c, Tt) - ref[k]
    print(f'  {lab}: RMSE {np.sqrt(np.mean(d**2)):.4f}, max {np.max(np.abs(d)):.4f} (units of the property)')
hfg = R.poly(R.C_hV, T[lq]) - R.poly(R.C_hL, T[lq]); hfg_ref = ref['hv'][lq] - ref['hL']
print(f'  latent heat from the fits: max relative error {100*np.max(np.abs(hfg/hfg_ref - 1)):.3f} %')

# ---------------- saturation pressure, Eq. (2), against IAPWS-95, 0 to 50 C
Ts = np.arange(0.0, 50.0 + 1e-9, 0.5)
ps_ref = np.array([CP.PropsSI('P', 'T', t + 273.15, 'Q', 0, 'Water') for t in Ts])
print(f'Saturation pressure, 0 to 50 C: max relative error {100*np.max(np.abs(R.psat(Ts)/ps_ref - 1)):.3f} %')

# ---------------- humidity ratio and wet-bulb temperature against ASHRAE at 101.325 kPa, 5 to 40 C
p = 101325.0
dw, dtw = [], []
for t in np.arange(5.0, 40.0 + 1e-9, 0.5):
    for rh in np.arange(5.0, 100.0 + 1e-9, 2.5):
        wm = float(R.w_from(rh/100*R.psat(t), p)); wa = PL.GetHumRatioFromRelHum(float(t), rh/100, p)
        dw.append(abs(wm/wa - 1))
        twa = PL.GetTWetBulbFromRelHum(float(t), rh/100, p)
        if twa > 0:       # below 0 C, psychrolib uses ice; the model uses liquid water
            twm = float(R.wet_bulb(np.array([t]), np.array([rh]), np.array([p]), pr)[0]); dtw.append(abs(twm - twa))
print(f'Humidity ratio, 5 to 40 C: max relative difference {100*max(dw):.3f} %')
print(f'Wet-bulb temperature, 5 to 40 C where it is above 0 C: max difference {max(dtw):.3f} K')

# ---------------- moist-air entropy differences, Eq. (11), against ASHRAE RP-1485 (CoolProp HumidAir)
def s_model(t, rh, mix=True):
    pv = rh/100*R.psat(t); w = float(R.w_from(pv, p))
    sa, sv = float(pr.sA(t)), float(pr.sV(t))
    if mix:
        sa -= R.R_A*np.log((p - pv)/R.P0); sv -= R.R_V*np.log(pv/R.psat(t))
    return sa + w*sv


def s_ref(t, rh):
    return CP.HAPropsSI('S', 'T', t + 273.15, 'P', p, 'R', rh/100)/1e3


base = (25.0, 30.0)
err_mix, err_nomix = [], []
for st in [(30.0, 100.0), (36.8, 100.0), (40.0, 100.0), (45.0, 100.0), (27.0, 100.0), (32.0, 90.0), (35.0, 95.0)]:
    dref = s_ref(*st) - s_ref(*base)
    err_mix.append(abs((s_model(*st) - s_model(*base))/dref - 1))
    err_nomix.append((s_model(*st, mix=False) - s_model(*base, mix=False))/dref - 1)
print(f'Moist-air entropy rise from 25 C, 30 % to saturated outlet states: mixture within {100*max(err_mix):.2f} %; '
      f'without the entropy of mixing {100*min(err_nomix):.1f} to {100*max(err_nomix):.1f} %')
