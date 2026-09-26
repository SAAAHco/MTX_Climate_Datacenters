# Variance-based global sensitivity analysis (Sobol, Saltelli sampling) and Latin-hypercube uncertainty propagation
# of the long-term mean PUE and WUE, revised model, climate-weighted over (T, RH) histograms.
import numpy as np, pandas as pd, json, time
from SALib.sample import saltelli
from SALib.analyze import sobol
from scipy.stats import qmc
import revmodel as R
import stations as S

d = pd.read_pickle('res_design.pkl')[['Region', 'St', 'T', 'RH']]
HIST = {}
for s in S.MADRID+S.TEXAS:
    g = d[d.St == s]
    Tb = np.round(g['T'].values/1.0); Hb = np.round(g['RH'].values/2.5)
    key = pd.DataFrame(dict(a=Tb, b=Hb, T=g['T'].values, H=g['RH'].values)).groupby(['a', 'b']).agg(T=('T', 'mean'), H=('H', 'mean'), n=('T', 'size'))
    HIST[s] = (key['T'].values, key['H'].values, key['n'].values.astype(float))
b = R.DEFAULT
problem = dict(num_vars=10,
               names=['cop_a0', 'cop_a1', 'cop_a2', 'dp_fan', 'dp_pump', 'a_ev', 'a_dry', 'f_cpL', 'f_hfg', 'f_cpa'],
               bounds=[[b['cop_a0']*0.9, b['cop_a0']*1.1], [b['cop_a1']*0.8, b['cop_a1']*1.2], [b['cop_a2']*0.8, b['cop_a2']*1.2],
                       [b['dp_fan']*0.5, b['dp_fan']*1.5], [b['dp_pump']*0.7, b['dp_pump']*1.3], [2.0, 6.0], [3.0, 9.0],
                       [0.995, 1.005], [0.995, 1.005], [0.995, 1.005]])
OUTS = ['Madrid', 'El Paso', 'Lubbock', 'Dallas', 'Houston']


def model(x):
    kw = dict(zip(problem['names'], x)); kw['ngrid'] = 41
    res = {}
    mad = []
    for s in S.MADRID+S.TEXAS:
        T, H, w = HIST[s]
        r = R.run(T, H, p=S.pressure(s), **kw)
        pu, wu = np.sum(w*r['PUE'])/w.sum(), np.sum(w*r['WUE'])/w.sum()
        if s in S.MADRID:
            mad.append((pu, wu, w.sum()))
        else:
            res[S.LABEL[s]] = (pu, wu)
    n = np.array([m[2] for m in mad])
    res['Madrid'] = (np.sum([m[0] for m in mad]*n)/n.sum(), np.sum([m[1] for m in mad]*n)/n.sum())
    return np.array([res[k][0] for k in OUTS] + [res[k][1] for k in OUTS])


def base_values():
    """The design case evaluated with the same histogram method (used to scale the intervals to the record-level values)."""
    x0 = [b['cop_a0'], b['cop_a1'], b['cop_a2'], b['dp_fan'], b['dp_pump'], b['a_ev'], b['a_dry'], 1.0, 1.0, 1.0]
    y = model(x0)
    return {lab: float(v) for lab, v in zip([f'PUE_{k}' for k in OUTS] + [f'WUE_{k}' for k in OUTS], y)}


if __name__ == '__main__':
    t0 = time.time()
    json.dump(base_values(), open('gsa_base.json', 'w'), indent=1)
    X = saltelli.sample(problem, 1024, calc_second_order=False)
    print('Saltelli samples', X.shape, flush=True)
    Y = np.empty((X.shape[0], 2*len(OUTS)))
    for i, x in enumerate(X):
        Y[i] = model(x)
        if i % 1000 == 0:
            print(i, f'{time.time()-t0:.0f}s', flush=True)
    np.save('gsa_Y.npy', Y); np.save('gsa_X.npy', X)
    SI = {}
    for j, lab in enumerate([f'PUE_{k}' for k in OUTS] + [f'WUE_{k}' for k in OUTS]):
        si = sobol.analyze(problem, Y[:, j], calc_second_order=False, print_to_console=False, seed=11)
        SI[lab] = {n: dict(S1=float(si['S1'][i]), S1c=float(si['S1_conf'][i]), ST=float(si['ST'][i]), STc=float(si['ST_conf'][i]))
                   for i, n in enumerate(problem['names'])}
    # Latin hypercube propagation
    lhs = qmc.LatinHypercube(d=10, seed=7).random(4096)
    lo = np.array([bb[0] for bb in problem['bounds']]); hi = np.array([bb[1] for bb in problem['bounds']])
    XL = lo + lhs*(hi-lo)
    YL = np.array([model(x) for x in XL])
    np.save('gsa_YL.npy', YL)
    U = {}
    for j, lab in enumerate([f'PUE_{k}' for k in OUTS] + [f'WUE_{k}' for k in OUTS]):
        U[lab] = dict(mean=float(YL[:, j].mean()), std=float(YL[:, j].std()),
                      p2_5=float(np.percentile(YL[:, j], 2.5)), p97_5=float(np.percentile(YL[:, j], 97.5)))
    json.dump(dict(problem=problem, SI=SI, U=U), open('gsa_out.json', 'w'), indent=1)
    print('done', f'{time.time()-t0:.0f}s')
    for lab in SI:
        top = sorted(SI[lab].items(), key=lambda kv: -kv[1]['ST'])[:4]
        print(lab, ' '.join(f"{n}:ST={v['ST']:.2f}/S1={v['S1']:.2f}" for n, v in top), '| 95%:', f"{U[lab]['p2_5']:.4f}-{U[lab]['p97_5']:.4f}")
