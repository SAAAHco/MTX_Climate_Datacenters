# Sensitivity cases, Pareto front over the water value lambda, and the lambda interval of the selection rule,
# all evaluated on the full weather records (the same records as Table 4), so that the design case of every
# table reproduces the long-term means of Table 4 exactly.
import numpy as np, pandas as pd, json, time
import revmodel as R
import stations as S

d = pd.read_pickle('res_design.pkl')[['Region', 'St', 'T', 'RH']]
REC = {s: (d[d.St == s]['T'].values.astype(float), d[d.St == s]['RH'].values.astype(float)) for s in S.MADRID+S.TEXAS}


def rh_const_w(T, RH, dT, p):
    """Relative humidity after a uniform warming dT at constant vapor pressure (constant specific humidity)."""
    pv = RH/100*R.psat(T)
    return np.clip(100*pv/R.psat(T+dT), 0.5, 100.0)


def evaluate(dT=0.0, const_w=False, **kw):
    out = {}
    tot = {}
    for s in S.MADRID+S.TEXAS:
        T, H = REC[s]
        p = kw.get('p_const', S.pressure(s))
        kk = {k: v for k, v in kw.items() if k != 'p_const'}
        H2 = rh_const_w(T, H, dT, p) if const_w else H
        r = R.run(T+dT, H2, p=p, **kk)
        n = len(T)
        occ = [float(np.mean(r['mode'] == m)) for m in range(3)]
        out[S.LABEL[s]] = dict(PUE=float(np.mean(r['PUE'])), WUE=float(np.mean(r['WUE'])), occ=occ,
                               reset=float(np.mean(r['reset'])), n=n,
                               fan=float(np.mean(r['fan'])), pump=float(np.mean(r['pump'])), chill=float(np.mean(r['chill'])))
    # regional pooled (record weighted), as the Madrid and Texas rows of Table 4
    for reg, sts in (('Madrid', S.MADRID), ('Texas', S.TEXAS)):
        n = np.array([out[S.LABEL[s]]['n'] for s in sts], float)
        out[reg] = {y: float(np.sum([out[S.LABEL[s]][y]*ni for s, ni in zip(sts, n)])/n.sum()) for y in ('PUE', 'WUE', 'reset', 'fan', 'pump', 'chill')}
        out[reg]['occ'] = [float(np.sum([out[S.LABEL[s]]['occ'][m]*ni for s, ni in zip(sts, n)])/n.sum()) for m in range(3)]
    return out


CASES = [('baseline', {}), ('limit', dict(a_ev=0.0, a_dry=0.0)),
         ('a_ev2', dict(a_ev=2.0)), ('a_ev6', dict(a_ev=6.0)), ('a_dry3', dict(a_dry=3.0)), ('a_dry9', dict(a_dry=9.0)),
         ('RH5_95', dict(RH5=95.0)), ('RH5_90', dict(RH5=90.0)), ('p_0.1MPa', dict(p_const=1.0e5)),
         ('no_comp_heat', dict(comp_heat=False)), ('cop_lift', dict(cop_model='lift')),
         ('cop_x0.8', dict(cop_scale=0.8)), ('cop_x1.2', dict(cop_scale=1.2)),
         ('range10', dict(dT_chiller=10.0)), ('range10_lift', dict(dT_chiller=10.0, cop_model='lift')),
         ('dp_fan150', dict(dp_fan=150.0)), ('dp_fan350', dict(dp_fan=350.0)), ('dp_pump0.35', dict(dp_pump=0.35e6)),
         ('fanJ_0.02', dict(lam=0.02, tower_lam=0.02)), ('fanJ_0.2', dict(lam=0.2, tower_lam=0.2)),
         ('warm+1', dict(dT=1.0)), ('warm+2', dict(dT=2.0)), ('warm+2_w', dict(dT=2.0, const_w=True))]
RES = {}
t0 = time.time()
for name, kw in CASES:
    kw = dict(kw); dT = kw.pop('dT', 0.0); cw = kw.pop('const_w', False)
    RES[name] = evaluate(dT=dT, const_w=cw, **kw)
    b = RES[name]
    print(f"{name:13s} " + ' '.join(f"{k}:{v['PUE']:.4f}/{v['WUE']:.3f}" for k, v in b.items() if k in ('Madrid', 'El Paso', 'Lubbock', 'Dallas', 'Houston')) + f"  [{time.time()-t0:.0f}s]", flush=True)
    json.dump(dict(cases=RES), open('sens_out_partial.json', 'w'), indent=1)
LAMS = [0.0, 0.001, 1/370, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 10.0, 1000.0]
PARETO = {}
for lam in LAMS:
    PARETO[repr(lam)] = evaluate(lam=lam)
    b = PARETO[repr(lam)]
    print(f"lam={lam:<8.4g} " + ' '.join(f"{k}:{v['PUE']:.4f}/{v['WUE']:.3f}" for k, v in b.items() if k in ('Madrid', 'El Paso', 'Lubbock', 'Dallas', 'Houston')) + f"  [{time.time()-t0:.0f}s]", flush=True)

# lambda interval that reproduces the rule, per record
iv = {}
for s in S.MADRID+S.TEXAS:
    T, H = REC[s]
    prm = dict(R.DEFAULT); pr = R.Props()
    p = S.pressure(s)
    Twb = R.wet_bulb(T, np.clip(H, 0.5, 100), np.full(T.shape, p), pr)
    Hc = np.clip(H, 0.5, 100)
    free_ok = T <= 18-prm['a_dry']; ev_ok = Twb <= 18-prm['a_ev']
    one = np.ones_like(T)
    fr = R.tower_states(27*one, 18*one, T, Hc, one, p, prm, pr, False)
    ev = R.tower_states(27*one, 18*one, T, Hc, one, p, prm, pr, True)
    T4c = np.maximum(26.0, Twb+4.0); cop = R.cop_fun(T4c, T4c+24, prm)
    ch = R.tower_states(T4c+24, T4c, T, Hc, 1+1/cop, p, prm, pr, True)
    aux = lambda st: st['VA5']*prm['dp_fan']/1000 + (st['mw']/1000)*prm['dp_pump']/1000
    Pf = 1+aux(fr); Pe = 1+aux(ev); Pc = 1+aux(ch)+1/cop
    both = free_ok & ev_ok
    lam_lo = (Pf-Pe)/ev['WUE']                    # free preferred iff lam > lam_lo
    m2 = (~free_ok) & ev_ok & (ev['WUE'] > ch['WUE'])
    lam_hi = (Pc-Pe)/(ev['WUE']-ch['WUE'])        # evap preferred iff lam < lam_hi
    iv[S.LABEL[s]] = dict(lam_lo_max=float(np.max(lam_lo[both])) if both.any() else None,
                          lam_lo_p99=float(np.percentile(lam_lo[both], 99)) if both.any() else None,
                          lam_hi_min=float(np.min(lam_hi[m2])) if m2.any() else None,
                          share_evap_wetter=float(np.sum(m2)/np.sum(ev_ok & ~free_ok)) if (ev_ok & ~free_ok).any() else 0.0)
print('lambda interval:', json.dumps(iv, indent=0))
json.dump(dict(cases=RES, pareto=PARETO, interval=iv, lams=LAMS, basis='records'), open('sens_out.json', 'w'), indent=1)
print('done', time.time()-t0)
