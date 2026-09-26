# Statistical analyses of the revised PUE/WUE series (mirrors the procedures of the V2 pipeline).
import numpy as np, pandas as pd, pickle, json
from scipy.ndimage import gaussian_filter
import revmodel as R
import stations as S

d = pd.read_pickle('res_design.pkl')
lim = pd.read_pickle('res_limit.pkl')
d['lim_PUE'] = lim['PUE'].values; d['lim_WUE'] = lim['WUE'].values; d['lim_mode'] = lim['mode'].values
d['cls'] = np.where(d['mode'] < 2, d['mode'], np.where(d['reset'] > 0, 3, 2)).astype(int)
OUT = {}
P_IT = 100.0  # MW


def stats_block(g):
    r = {}
    for y in ('PUE', 'WUE'):
        x = g[y].values
        q = np.percentile(x, [0, 25, 50, 75, 95, 99, 100])
        r[y] = dict(mean=x.mean(), std=x.std(), p0=q[0], p25=q[1], p50=q[2], p75=q[3], p95=q[4], p99=q[5], p100=q[6])
    occ = np.bincount(g['cls'].values, minlength=4)/len(g)*100
    r['occ'] = dict(free=occ[0], evap=occ[1], chill=occ[2], reset=occ[3])
    r['lim'] = dict(PUE=g['lim_PUE'].mean(), WUE=g['lim_WUE'].mean(),
                    occ=(np.bincount(g['lim_mode'].values.astype(int), minlength=3)/len(g)*100).tolist())
    r['n'] = len(g)
    r['maxT4'] = float(g['T4'].max()); r['maxT1'] = float(g['T1'].max())
    return r


OUT['stats'] = {S.LABEL[s]: stats_block(d[d.St == s]) for s in S.MADRID+S.TEXAS}
OUT['stats']['Madrid'] = stats_block(d[d.Region == 'Madrid'])
OUT['stats']['Texas'] = stats_block(d[d.Region == 'Texas'])


def nested_var(g, y, spatial=True):
    x = g[y].values; mu = x.mean(); tot = x.var()
    keys = (['St'] if spatial else []) + ['Year', 'Month', 'Slot']
    lv = []
    for k in range(1, len(keys)+1):
        lv.append(g.groupby(keys[:k])[y].transform('mean').values)
    comps = [np.mean((lv[0]-mu)**2)] + [np.mean((lv[i]-lv[i-1])**2) for i in range(1, len(lv))] + [np.mean((x-lv[-1])**2)]
    names = (['spatial'] if spatial else []) + ['interannual', 'seasonal', 'intraday', 'interday']
    return {n: 100*c/tot for n, c in zip(names, comps)}, sum(comps)/tot


OUT['vardec'] = {}
for reg in ('Madrid', 'Texas'):
    g = d[d.Region == reg]
    OUT['vardec'][reg] = {y: nested_var(g, y)[0] for y in ('PUE', 'WUE')}
for s in S.TEXAS:
    g = d[d.St == s]
    OUT['vardec'][S.LABEL[s]] = {y: nested_var(g, y, spatial=False)[0] for y in ('PUE', 'WUE')}

# averaging-horizon convergence
conv = {}
for s in S.MADRID+S.TEXAS:
    g = d[d.St == s].set_index('ts').sort_index()
    g = g[~g.index.duplicated(keep='first')]
    rs = {}
    for y in ('PUE', 'WUE'):
        v0 = g[y].var(); rr = {}
        for lab, rule in [('1h', 'h'), ('1d', 'D'), ('1wk', '7D'), ('1mo', 'MS'), ('1yr', 'YS')]:
            agg = g[y].resample(rule).agg(['mean', 'count'])
            med = agg['count'].median(); agg = agg[agg['count'] >= 0.5*med]
            rr[lab] = agg['mean'].var()/v0
        rs[y] = rr
    ann = g.groupby(g.index.year)[['PUE', 'WUE']].mean()
    ann = ann[g.groupby(g.index.year).size() > 0.9*g.groupby(g.index.year).size().median()]
    rs['annual'] = dict(PUE_std=ann['PUE'].std(), WUE_std=ann['WUE'].std(), PUE_min=ann['PUE'].min(), PUE_max=ann['PUE'].max(),
                        WUE_min=ann['WUE'].min(), WUE_max=ann['WUE'].max(), nyears=len(ann))
    conv[S.LABEL[s]] = rs
OUT['conv'] = conv

# mode occupancy and conditional footprint per region and station
modes = {}
for key, g in [('Madrid', d[d.Region == 'Madrid']), ('Texas', d[d.Region == 'Texas'])] + [(S.LABEL[s], d[d.St == s]) for s in S.TEXAS]:
    mm = {}
    for c, lab in enumerate(['free', 'evap', 'chill', 'reset']):
        gm = g[g.cls == c]
        mm[lab] = dict(occ=100*len(gm)/len(g), PUE=gm['PUE'].mean() if len(gm) else np.nan, WUE=gm['WUE'].mean() if len(gm) else np.nan)
    modes[key] = mm
OUT['modes'] = modes

# climate transfer: kernel densities and deterministic maps
Tg = np.arange(-15.0, 48.0+1e-9, 0.5); Hg = np.arange(0.5, 100.0, 1.0)
TT, HH = np.meshgrid(Tg, Hg, indexing='ij')
dens = {}; maps = {}
bw = {}
for reg, sts in (('Madrid', S.MADRID), ('Texas', S.TEXAS)):
    hT = np.mean([d[d.St == s]['T'].std()*len(d[d.St == s])**(-1/6) for s in sts])
    hR = np.mean([d[d.St == s]['RH'].std()*len(d[d.St == s])**(-1/6) for s in sts])
    bw[reg] = (hT, hR)
    for s in sts:
        g = d[d.St == s]
        Hc, _, _ = np.histogram2d(g['T'].values, g['RH'].values, bins=[np.append(Tg-0.25, Tg[-1]+0.25), np.arange(0, 101, 1.0)])
        f = gaussian_filter(Hc, sigma=(hT/0.5, hR/1.0), mode='constant')
        dens[s] = f/f.sum()
        r = R.run(TT.ravel(), HH.ravel(), p=S.pressure(s))
        maps[s] = {y: r[y].reshape(TT.shape) for y in ('PUE', 'WUE')}
OUT['bandwidth'] = {k: [float(v[0]), float(v[1])] for k, v in bw.items()}


def hell(a, b):
    return float(np.sqrt(max(0.0, 1-np.sum(np.sqrt(a*b)))))


def pred(target, donor_f):
    return {y: float(np.sum(maps[target][y]*donor_f)/np.sum(donor_f)) for y in ('PUE', 'WUE')}


emp = {s: {y: d[d.St == s][y].mean() for y in ('PUE', 'WUE')} for s in S.MADRID+S.TEXAS}
err = lambda p, e: 100*abs(p-e)/e
trans = {}
for s in S.MADRID:
    others = sum(dens[o] for o in S.MADRID if o != s); others = others/others.sum()
    po, pp = pred(s, dens[s]), pred(s, others)
    trans[S.LABEL[s]] = dict(own={y: err(po[y], emp[s][y]) for y in po}, loo={y: err(pp[y], emp[s][y]) for y in pp})
pairs = {'ElPaso': 'Lubbock', 'Lubbock': 'ElPaso', 'Dallas': 'Houston', 'Houston': 'Dallas'}
for s in S.TEXAS:
    others = sum(dens[o] for o in S.TEXAS if o != s); others = others/others.sum()
    po, pa, pp = pred(s, dens[s]), pred(s, dens[pairs[s]]), pred(s, others)
    trans[S.LABEL[s]] = dict(own={y: err(po[y], emp[s][y]) for y in po}, pair={y: err(pa[y], emp[s][y]) for y in pa},
                             pooled={y: err(pp[y], emp[s][y]) for y in pp})
OUT['transfer'] = trans
HM = {}
for reg, sts in (('Madrid', S.MADRID), ('Texas', S.TEXAS)):
    vals = [hell(dens[a], dens[b]) for i, a in enumerate(sts) for b in sts[i+1:]]
    HM[reg] = (min(vals), max(vals))
HM['pairs'] = {f'{S.LABEL[a]}-{S.LABEL[b]}': hell(dens[a], dens[b]) for a, b in
               [('ElPaso', 'Lubbock'), ('Dallas', 'Houston'), ('ElPaso', 'Dallas'), ('ElPaso', 'Houston'), ('Lubbock', 'Dallas'), ('Lubbock', 'Houston')]}
OUT['hellinger'] = HM
pickle.dump(dict(Tg=Tg, Hg=Hg, dens=dens, maps=maps), open('transfer_maps.pkl', 'wb'))

# annual demands, monthly calendar, ET0 (FAO-56 Hargreaves), peaks, interannual bounds
days = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]); hm = days*24.0; H = hm.sum()
ann = {}
for key, g, lat in [('Madrid', d[d.Region == 'Madrid'], np.mean([S.LAT[s] for s in S.MADRID]))] + \
        [(S.LABEL[s], d[d.St == s], S.LAT[s]) for s in S.TEXAS]:
    mm = g.groupby('Month')[['PUE', 'WUE']].mean().reindex(range(1, 13))
    E_m = P_IT*hm*mm['PUE'].values/1e3            # GWh
    W_m = P_IT*1e3*hm*mm['WUE'].values/1e6        # ML (1e6 L)
    E_y = P_IT*H*g['PUE'].mean()/1e3; W_y = P_IT*1e3*H*g['WUE'].mean()/1e6
    # ET0 from daily Tmax/Tmin per station, then averaged over stations of the region
    gg = g[['St', 'ts', 'T']].copy(); gg['date'] = gg['ts'].dt.normalize()
    dly = gg.groupby(['St', 'date'])['T'].agg(['max', 'min', 'mean', 'count']).reset_index()
    dly = dly[dly['count'] >= 0.75*dly['count'].median()]
    doy = dly['date'].dt.dayofyear.values; phi = np.radians(np.array([S.LAT[s] for s in dly['St'].values]))
    dr = 1+0.033*np.cos(2*np.pi*doy/365); dec = 0.409*np.sin(2*np.pi*doy/365-1.39)
    ws = np.arccos(np.clip(-np.tan(phi)*np.tan(dec), -1, 1))
    Ra = (24*60/np.pi)*0.0820*dr*(ws*np.sin(phi)*np.sin(dec)+np.cos(phi)*np.cos(dec)*np.sin(ws))
    et = 0.0023*0.408*Ra*(dly['mean'].values+17.8)*np.sqrt(np.clip(dly['max'].values-dly['min'].values, 0, None))
    dly['ET0'] = et; dly['month'] = dly['date'].dt.month
    ET_m = dly.groupby('month')['ET0'].mean().reindex(range(1, 13)).values*days
    r = float(np.corrcoef(W_m, ET_m)[0, 1])
    # annual totals per calendar year (interannual range)
    yy = g.groupby('Year')[['PUE', 'WUE']].mean(); cnt = g.groupby('Year').size(); yy = yy[cnt > 0.9*cnt.median()]
    E_years = P_IT*H*yy['PUE'].values/1e3; W_years = P_IT*1e3*H*yy['WUE'].values/1e6
    # peaks (record level): cooling electricity (MW) and water rate (m3/h) for 100 MW ITE
    cool = (g['PUE'].values-1)*P_IT; wr = g['WUE'].values*P_IT*1e3/1e3
    # daily totals
    gd = g.copy(); gd['date'] = gd['ts'].dt.normalize()
    daily = gd.groupby(['St', 'date'])[['PUE', 'WUE']].mean()
    Wd = daily['WUE'].values*P_IT*1e3*24/1e3            # m3/day
    Ed = daily['PUE'].values*P_IT*24/1e3                # GWh/day
    ann[key] = dict(E_y=E_y, W_y=W_y, E_over=E_y-P_IT*H/1e3, E_m=E_m.tolist(), W_m=W_m.tolist(), ET_m=ET_m.tolist(), r=r,
                    summer_winter=float(W_m[[5, 6, 7]].mean()/W_m[[11, 0, 1]].mean()), aprsep=float(100*W_m[3:9].sum()/W_m.sum()),
                    peakW=int(np.argmax(W_m))+1, peakET=int(np.argmax(ET_m))+1,
                    E_years=(float(E_years.min()), float(E_years.max()), float(E_years.std())),
                    W_years=(float(W_years.min()), float(W_years.max()), float(W_years.std())),
                    cool_MW=dict(mean=float(cool.mean()), p95=float(np.percentile(cool, 95)), p99=float(np.percentile(cool, 99)), max=float(cool.max())),
                    water_m3h=dict(mean=float(wr.mean()), p95=float(np.percentile(wr, 95)), p99=float(np.percentile(wr, 99)), max=float(wr.max())),
                    water_m3d=dict(mean=float(Wd.mean()), p95=float(np.percentile(Wd, 95)), p99=float(np.percentile(Wd, 99)), max=float(Wd.max())),
                    PUE_p=dict(p95=float(np.percentile(g['PUE'], 95)), p99=float(np.percentile(g['PUE'], 99))),
                    WUE_p=dict(p95=float(np.percentile(g['WUE'], 95)), p99=float(np.percentile(g['WUE'], 99))),
                    ha=float(W_y*1000/6500))
OUT['annual'] = ann
json.dump(OUT, open('analysis_out.json', 'w'), indent=1, default=float)


def pr(x, n=4):
    return f'{x:.{n}f}'


print('== Descriptive (design) ==')
for k, v in OUT['stats'].items():
    print(f"{k:8s} PUE {pr(v['PUE']['mean'])} sd {pr(v['PUE']['std'],3)} p50 {pr(v['PUE']['p50'])} p99 {pr(v['PUE']['p99'])} max {pr(v['PUE']['p100'])} | "
          f"WUE {pr(v['WUE']['mean'],3)} sd {pr(v['WUE']['std'],2)} p50 {pr(v['WUE']['p50'],2)} p99 {pr(v['WUE']['p99'],2)} max {pr(v['WUE']['p100'],2)} | "
          f"occ F/E/C/Cr {v['occ']['free']:.1f}/{v['occ']['evap']:.1f}/{v['occ']['chill']:.1f}/{v['occ']['reset']:.1f} | limit PUE {pr(v['lim']['PUE'])} WUE {pr(v['lim']['WUE'],3)} | maxT4 {v['maxT4']:.1f}")
print('== Variance decomposition ==')
for k, v in OUT['vardec'].items():
    print(k, {y: {a: round(b, 2) for a, b in v[y].items()} for y in v})
print('== Convergence (Madrid mean, Texas per station) ==')
for k, v in OUT['conv'].items():
    print(k, {y: {a: round(b, 3) for a, b in v[y].items()} for y in ('PUE', 'WUE')}, {a: round(b, 4) for a, b in v['annual'].items()})
print('== Modes ==')
for k, v in OUT['modes'].items():
    print(k, {a: (round(b['occ'], 1), round(b['PUE'], 4), round(b['WUE'], 3)) for a, b in v.items()})
print('== Transfer (errors %) ==', OUT['bandwidth'])
for k, v in OUT['transfer'].items():
    print(k, {a: {y: round(z, 2) for y, z in b.items()} for a, b in v.items()})
print('== Hellinger ==', OUT['hellinger'])
print('== Annual (100 MW) ==')
for k, v in OUT['annual'].items():
    print(k, {a: (round(b, 3) if isinstance(b, float) else b) for a, b in v.items() if a not in ('E_m', 'W_m', 'ET_m')})
