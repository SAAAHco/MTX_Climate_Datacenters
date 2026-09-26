# Baseline (design case) and thermodynamic-limit runs on the complete series, per record.
import numpy as np, pandas as pd, time, sys
import revmodel as R
from loaddata import madrid, texas
import stations as S

KEEP = ['mode', 'reset', 'PUE', 'WUE', 'fan', 'pump', 'chill', 'T1', 'T4', 'T5', 'Twb', 'Qtower', 'VA5']


def frame():
    m = madrid().copy(); t = texas().copy()
    mins = (m['HM']//100)*60 + m['HM'] % 100
    m['ts'] = pd.to_datetime(m['Fecha'], format='%d/%m/%Y') + pd.to_timedelta(mins, unit='m')
    m['Slot'] = ((mins % 1440)//30).astype(int)
    m['Region'] = 'Madrid'
    t['ts'] = pd.to_datetime(dict(year=t.Year, month=t.Month, day=t.Day, hour=t.Hour, minute=t.Minute))
    t['Slot'] = t['Hour'].astype(int); t['Region'] = 'Texas'
    cols = ['Region', 'St', 'ts', 'Year', 'Month', 'Slot', 'T', 'RH']
    m['St'] = m['St'].astype(object); t['St'] = t['St'].astype(object)
    d = pd.concat([m[cols], t[cols]], ignore_index=True)
    return d


def run_case(d, tag, **kw):
    out = {k: np.empty(len(d)) for k in KEEP}
    t0 = time.time()
    for st in S.MADRID + S.TEXAS:
        idx = np.where(d['St'].values == st)[0]
        r = R.run(d['T'].values[idx], d['RH'].values[idx], p=S.pressure(st), **kw)
        for k in KEEP:
            out[k][idx] = r[k]
        print(f'  {tag} {S.LABEL[st]:8s} PUE {np.mean(r["PUE"]):.4f} WUE {np.mean(r["WUE"]):.4f} '
              f'modes {np.bincount(r["mode"], minlength=3)/len(idx)*100} reset {100*np.mean(r["reset"]):.1f}% '
              f'[{time.time()-t0:.0f}s]', flush=True)
    res = d.copy()
    for k in KEEP:
        res[k] = out[k]
    return res


if __name__ == '__main__':
    d = frame()
    base = run_case(d, 'design')
    base.to_pickle('res_design.pkl')
    lim = run_case(d, 'limit', a_ev=0.0, a_dry=0.0)
    lim[['PUE', 'WUE', 'mode']].to_pickle('res_limit.pkl')
