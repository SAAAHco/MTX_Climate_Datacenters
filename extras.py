# Additional diagnostics for the manuscript (written to extras_out.json):
# data completeness, drift share, Merkel number and L/G of the design states, air-flow sizing ratios,
# marginal water value of extra tower air, Carnot fraction of the chiller curves, all-pair climate transfer,
# ET0 totals, evaporator-side pumping.
import numpy as np, pandas as pd, json, pickle
import revmodel as R
import stations as S
from config import MADRID_DATA

OUT = {}
d = pd.read_pickle('res_design.pkl')
pr = R.Props(); prm = dict(R.DEFAULT)

# ---------------------------------------------------------------- data completeness
raw = pd.read_csv(MADRID_DATA, sep='\t', usecols=[1, 2, 3, 4, 5, 6, 7], encoding='latin-1')
raw.columns = ['St', 'Fecha', 'Year', 'Doy', 'HM', 'T', 'RH']
exp_steps = (pd.Timestamp('2019-12-31') - pd.Timestamp('2006-01-01')).days + 1
exp_steps *= 48
comp = {}
for s, g in raw.groupby('St'):
    key = (pd.to_datetime(g['Fecha'], format='%d/%m/%Y') + pd.to_timedelta(g['HM']//100, 'h') + pd.to_timedelta(g['HM'] % 100, 'min'))
    comp[str(s)] = dict(rows=int(len(g)), expected=int(exp_steps), absent=int(exp_steps-key.nunique()), duplicated=int(key.duplicated().sum()),
                        noT=int(g['T'].isna().sum()), noRH=int(g['RH'].isna().sum()), both=int((g['T'].isna() & g['RH'].isna()).sum()),
                        first=str(key.min()), last=str(key.max()))
OUT['madrid_completeness'] = comp
OUT['madrid_totals'] = dict(rows=int(len(raw)), expected=int(7*exp_steps), absent=int(sum(v['absent'] for v in comp.values())),
                            noT=int(raw['T'].isna().sum()), noRH=int(raw['RH'].isna().sum()), both=int((raw['T'].isna() & raw['RH'].isna()).sum()),
                            onlyT=int((raw['T'].isna() & raw['RH'].notna()).sum()), onlyRH=int((raw['RH'].isna() & raw['T'].notna()).sum()))
tx = pd.read_pickle('texas.pkl')
txc = {}
for s, g in tx.groupby('St'):
    ts = pd.to_datetime(dict(year=g.Year, month=g.Month, day=g.Day, hour=g.Hour, minute=g.Minute))
    dt = np.diff(np.sort(ts.values)).astype('timedelta64[m]').astype(int)
    txc[str(s)] = dict(n=int(len(g)), first=str(ts.min()), last=str(ts.max()), steps=sorted(set(dt.tolist())), dup=int(ts.duplicated().sum()))
OUT['texas_completeness'] = txc

# ---------------------------------------------------------------- drift share (0.005 % of the circulating flow)
hL = pr.hL
mw = d['Qtower'].values/(hL(d['T1'].values)-hL(d['T4'].values))           # kg/s per kW
evap = d['WUE'].values/3600.0                                               # kg/s per kW
drift = 5e-5*mw
ds = {}
for key, m in (('evap', d['mode'].values == 1), ('chill', d['mode'].values == 2)):
    sh = 100*drift[m]/evap[m]
    ds[key] = dict(p1=float(np.percentile(sh, 1)), p50=float(np.percentile(sh, 50)), p99=float(np.percentile(sh, 99)), max=float(sh.max()), mw=float(np.median(mw[m])))
for site, g in [('Madrid', d.Region == 'Madrid')] + [(S.LABEL[s], d.St == s) for s in S.TEXAS]:
    m = g.values & (d['mode'].values > 0)
    ds['annual_'+site] = float(100*drift[m].sum()/evap[m].sum())
OUT['drift'] = ds

# ---------------------------------------------------------------- Merkel number and L/G of the design states (subsample)
rng = np.random.default_rng(1)
mer = {}
for s in S.MADRID+S.TEXAS:
    g = d[(d.St == s) & (d['mode'] > 0)]
    idx = rng.choice(len(g), size=min(4000, len(g)), replace=False)
    g = g.iloc[idx]
    p = S.pressure(s)
    for mode in (1, 2):
        gg = g[g['mode'] == mode]
        if len(gg) == 0:
            continue
        T1, T4, T2, RH2, Q = [gg[c].values.astype(float) for c in ('T1', 'T4', 'T', 'RH', 'Qtower')]
        st = R.tower_states(T1, T4, T2, RH2, Q, p, prm, pr, True)
        pv2 = RH2/100*R.psat(T2); w2 = R.w_from(pv2, p); H2 = pr.hA(T2)+w2*pr.hV(T2)
        f = np.linspace(0, 1, 401)
        Tw = T4[:, None] + f[None, :]*(T1-T4)[:, None]
        hs = pr.hA(Tw) + R.w_from(R.psat(Tw), p)*pr.hV(Tw)
        ha = H2[:, None] + (st['mw'][:, None]*(hL(Tw)-hL(T4)[:, None]))/st['ma'][:, None]
        integrand = (hL(Tw[:, 1:])-hL(Tw[:, :-1]))/(0.5*((hs-ha)[:, 1:]+(hs-ha)[:, :-1]))
        Me = integrand.sum(axis=1)
        LG = st['mw']/st['ma']
        mer.setdefault(mode, {'Me': [], 'LG': []})
        mer[mode]['Me'] += Me.tolist(); mer[mode]['LG'] += LG.tolist()
OUT['merkel'] = {('evap' if k == 1 else 'chill'): dict(Me_p1=float(np.percentile(v['Me'], 1)), Me_p50=float(np.percentile(v['Me'], 50)),
                                                     Me_p99=float(np.percentile(v['Me'], 99)), Me_max=float(np.max(v['Me'])),
                                                     LG_p1=float(np.percentile(v['LG'], 1)), LG_p50=float(np.percentile(v['LG'], 50)), LG_p99=float(np.percentile(v['LG'], 99)))
                 for k, v in mer.items()}
# Merkel number of the Figure 4 design state (T1 50, T4 26, 25 C / 30 %)
c0 = float(R.cop_fun(26.0, 50.0, prm)); Q0 = 1+1/c0
st = R.tower_states(np.array([50.0]), np.array([26.0]), np.array([25.0]), np.array([30.0]), np.array([Q0]), 101325.0, prm, pr, True)
Tw = np.linspace(26, 50, 2001); pv2 = 0.3*R.psat(25.0); w2 = R.w_from(pv2, 101325.0); H2 = pr.hA(25.0)+w2*pr.hV(25.0)
hs = pr.hA(Tw)+R.w_from(R.psat(Tw), 101325.0)*pr.hV(Tw); ha = H2+st['mw'][0]*(hL(Tw)-hL(26.0))/st['ma'][0]
OUT['merkel_fig4'] = dict(Me=float(np.sum((hL(Tw[1:])-hL(Tw[:-1]))/(0.5*((hs-ha)[1:]+(hs-ha)[:-1])))), LG=float(st['mw'][0]/st['ma'][0]))

# ---------------------------------------------------------------- air-flow sizing ratios (volumetric air flow per kW ITE)
af = {}
for site, g in [('Madrid', d[d.Region == 'Madrid'])] + [(S.LABEL[s], d[d.St == s]) for s in S.TEXAS]:
    a = {}
    for key, m in (('free', g['mode'] == 0), ('evap', g['mode'] == 1), ('chill', g['mode'] == 2), ('all', g['mode'] >= 0)):
        v = g['VA5'].values[m.values]
        if len(v):
            a[key] = dict(mean=float(v.mean()), p50=float(np.median(v)), p99=float(np.percentile(v, 99)), max=float(v.max()))
    a['max_over_mean'] = float(g['VA5'].max()/g['VA5'].mean()); a['p99_over_mean'] = float(np.percentile(g['VA5'], 99)/g['VA5'].mean())
    af[site] = a
OUT['airflow'] = af

# ---------------------------------------------------------------- marginal water value of extra tower air along saturated outlet states
mv = {}
for s in S.MADRID+S.TEXAS:
    g = d[(d.St == s) & (d['mode'] > 0)].sample(n=3000, random_state=2)
    p = S.pressure(s)
    for mode in (1, 2):
        gg = g[g['mode'] == mode]
        if len(gg) == 0:
            continue
        T1, T4, T2, RH2, Q = [gg[c].values.astype(float) for c in ('T1', 'T4', 'T', 'RH', 'Qtower')]
        st = R.tower_states(T1, T4, T2, RH2, Q, p, prm, pr, True)
        lamstar = []
        for kf in (1.1, 1.5, 2.0):
            ma = kf*st['ma']
            lo = np.full(ma.shape, -20.0); hi = T1+10
            pv2 = RH2/100*R.psat(T2); w2 = R.w_from(pv2, p); H2 = pr.hA(T2)+w2*pr.hV(T2); hL3 = hL(T2)
            for _ in range(55):
                m = 0.5*(lo+hi); w5 = R.w_from(R.psat(m), p)
                gg_ = ma*(pr.hA(m)+w5*pr.hV(m)-H2-(w5-w2)*hL3)-Q
                hi = np.where(gg_ > 0, m, hi); lo = np.where(gg_ > 0, lo, m)
            T5 = 0.5*(lo+hi); pv5 = R.psat(T5); w5 = R.w_from(pv5, p)
            WUEk = (w5-w2)*ma*3600; VAk = ma*R.R_A*1e3*(R.K+T5)/(p-pv5)
            dW = st['WUE']-WUEk; dP = (VAk-st['VA5'])*prm['dp_fan']/1000
            lamstar.append(np.where(dW > 0, dP/np.maximum(dW, 1e-12), np.inf))
        mv.setdefault(mode, []).append(np.vstack(lamstar))
OUT['marginal_lambda'] = {}
for mode, arrs in mv.items():
    A = np.hstack(arrs)
    OUT['marginal_lambda']['evap' if mode == 1 else 'chill'] = {f'k{k}': dict(p1=float(np.percentile(A[i], 1)), p5=float(np.percentile(A[i], 5)), p50=float(np.percentile(A[i], 50)),
                                                                              min=float(np.min(A[i])))
                                                                 for i, k in enumerate((1.1, 1.5, 2.0))}
    # water saved per kWh of extra fan energy (L/kWh) for a 10 % air increase
    OUT['marginal_lambda']['evap' if mode == 1 else 'chill']['L_per_kWh_k1.1'] = dict(p50=float(1/np.percentile(A[0], 50)), p99=float(1/np.percentile(A[0], 1)))

# ---------------------------------------------------------------- Carnot fraction of Eq. (16) with the 24 K and 10 K condenser ranges
cf = {}
for rng_ in (24.0, 10.0):
    for T4 in (26.0, 29.0, 33.0):
        T1 = T4+rng_
        c16 = float(R.cop_fun(T4, T1, prm)); Tev = 18-2+R.K; Tcd = T1+1.5+R.K
        cf[f'range{int(rng_)}_T4_{int(T4)}'] = dict(cop16=c16, carnot=float(Tev/(Tcd-Tev)), frac=float(c16/(Tev/(Tcd-Tev))),
                                                   cop_lift=float(R.cop_fun(T4, T1, dict(prm, cop_model='lift'))))
cf['eta_lift'] = float(R.LIFT_ETA)
OUT['carnot'] = cf

# ---------------------------------------------------------------- evaporator-side pumping in the chilling mode
mw_ite = 1.0/(hL(27.0)-hL(18.0))
OUT['evap_side_pump'] = {f'{dp:.2f}MPa': float(mw_ite/1000*dp*1e6/1000) for dp in (0.03, 0.05, 0.1)}
OUT['evap_side_pump']['mw_ite'] = float(mw_ite)

# ---------------------------------------------------------------- climate transfer, all donor-target pairs
T_ = pickle.load(open('transfer_maps.pkl', 'rb'))
dens, maps = T_['dens'], T_['maps']
emp = {s: {y: float(d[d.St == s][y].mean()) for y in ('PUE', 'WUE')} for s in S.MADRID+S.TEXAS}
hell = lambda a, b: float(np.sqrt(max(0.0, 1-np.sum(np.sqrt(a*b)))))
pred = lambda t, f: {y: float(np.sum(maps[t][y]*f)/np.sum(f)) for y in ('PUE', 'WUE')}
rows = []
for reg, sts in (('Madrid', S.MADRID), ('Texas', S.TEXAS)):
    for t in sts:
        donors = [(S.LABEL[o], dens[o]) for o in sts if o != t]
        pooled = sum(dens[o] for o in sts if o != t); pooled = pooled/pooled.sum()
        donors.append(('pooled', pooled))
        for nm, f in donors:
            pp = pred(t, f)
            rows.append(dict(region=reg, target=S.LABEL[t], donor=nm, H=hell(dens[t], f),
                             ePUE=100*abs(pp['PUE']-emp[t]['PUE'])/emp[t]['PUE'], eWUE=100*abs(pp['WUE']-emp[t]['WUE'])/emp[t]['WUE']))
OUT['transfer_all'] = rows
tr = pd.DataFrame(rows)
for reg in ('Madrid', 'Texas'):
    g = tr[tr.region == reg]
    OUT[f'transfer_corr_{reg}'] = dict(r_PUE=float(np.corrcoef(g['H'], g['ePUE'])[0, 1]), r_WUE=float(np.corrcoef(g['H'], g['eWUE'])[0, 1]),
                                       H_range=[float(g['H'].min()), float(g['H'].max())])

# ---------------------------------------------------------------- annual ET0 and equivalent reference-crop area
A = json.load(open('analysis_out.json'))['annual']
OUT['et0'] = {k: dict(ET0_mm=float(np.sum(v['ET_m'])), W_ML=float(v['W_y']), area_ha=float(v['W_y']*1000/(np.sum(v['ET_m'])*10.0))) for k, v in A.items()}
json.dump(OUT, open('extras_out.json', 'w'), indent=1, default=float)
print(json.dumps({k: v for k, v in OUT.items() if k not in ('transfer_all', 'madrid_completeness')}, indent=1, default=float)[:6000])
print(tr.pivot_table(index=['region', 'target'], columns='donor', values=['H', 'ePUE', 'eWUE']).round(3).to_string()[:4000])
