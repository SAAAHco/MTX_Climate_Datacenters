# Access to the numerical results for the manuscript and the response letter (all numbers come from here).
import json, numpy as np
A = json.load(open('analysis_out.json'))
S = json.load(open('sens_out.json'))
try:
    G = json.load(open('gsa_out.json'))
    GB = json.load(open('gsa_base.json'))      # design case evaluated with the histogram of the Sobol runs
except FileNotFoundError:
    G = None; GB = None
try:
    E = json.load(open('extras_out.json'))
except FileNotFoundError:
    E = None
MAD = ['M01', 'M02', 'M03', 'M04', 'M05', 'M06', 'M102']
TX = ['El Paso', 'Lubbock', 'Dallas', 'Houston']
SITES = ['Madrid'] + TX


def f(x, n=2):
    return f'{x:,.{n}f}'


def st(site):
    return A['stats'][site]


def mean(site, y):
    return st(site)[y]['mean']


def rng(vals, n=2, sep=' to '):
    lo, hi = min(vals), max(vals)
    a, b = f'{lo:.{n}f}', f'{hi:.{n}f}'
    return a if a == b else f'{a}{sep}{b}'


def sgn(x, n=1):
    r = round(x, n)
    if r == 0:
        return f'{0:.{n}f}'          # no signed zero
    return f'{r:+.{n}f}'.replace('-', '−')


def case(name, site, y):
    return S['cases'][name][site][y]


def gsa_interval(y, site):
    """95 % parameter interval of the long-term mean, scaled to the record-level design value (the Sobol runs use
    histograms that reproduce the record means within 0.5 %)."""
    base = st(site)[y]['mean']; gb = GB[f'{y}_{site}']; u = G['U'][f'{y}_{site}']
    if y == 'PUE':      # scale the non-ITE part
        return 1 + (u['p2_5']-1)*(base-1)/(gb-1), 1 + (u['p97_5']-1)*(base-1)/(gb-1)
    return u['p2_5']*base/gb, u['p97_5']*base/gb


def lam_interval():
    iv = S['interval']
    lo = max(v['lam_lo_max'] for v in iv.values()); hi = min(v['lam_hi_min'] for v in iv.values())
    los = [v['lam_lo_max'] for v in iv.values()]; his = [v['lam_hi_min'] for v in iv.values()]
    wet = [v['share_evap_wetter'] for v in iv.values()]
    return lo, hi, los, his, wet


def lam_fmt():
    """The water-value interval of the selection rule, rounded inward so that every value quoted lies inside it:
    lower end rounded up to 3 decimals, upper end rounded down to 2 decimals."""
    lo, hi = lam_interval()[:2]
    import math
    return f'{math.ceil(lo*1000-1e-9)/1000:.3f}', f'{math.floor(hi*100+1e-9)/100:.2f}'


def up(x, n=1):
    """Upper bound rounded up at n decimals (for 'within' and 'at most' statements)."""
    import math
    return f'{math.ceil(x*10**n-1e-9)/10**n:.{n}f}'


def down(x, n=1):
    import math
    return f'{math.floor(x*10**n+1e-9)/10**n:.{n}f}'


def pareto(site, lam):
    keys = list(S['pareto'].keys())
    k = min(keys, key=lambda kk: abs((float(kk) if kk != 'inf' else 1e9) - lam))
    return S['pareto'][k][site]


def pct(a, b):
    return 100*(a-b)/b
