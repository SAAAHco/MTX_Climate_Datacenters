# Figures of the revised manuscript.
import numpy as np, pandas as pd, json, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon
from matplotlib.colors import BoundaryNorm, ListedColormap
import revmodel as R
import stations as S

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 9.5, 'axes.labelsize': 9.5,
                     'xtick.labelsize': 8.5, 'ytick.labelsize': 8.5, 'legend.fontsize': 8, 'figure.dpi': 300,
                     'savefig.dpi': 300, 'savefig.bbox': 'tight', 'axes.spines.top': False, 'axes.spines.right': False})
OUT = 'figs/'
EPS = 'figs/eps/'            # vector copies for the journal: EPS (submitted) and PDF (used only to check the vector output)
import os, warnings
os.makedirs(OUT, exist_ok=True); os.makedirs(EPS, exist_ok=True)
plt.rcParams.update({'ps.fonttype': 42, 'pdf.fonttype': 42})


def save(fig, name):
    """Save the PNG used in the manuscript and the EPS and PDF copies of the same figure. EPS has no transparency, so
    any non-rasterized partially transparent artist stops the export instead of being drawn opaque."""
    fig.savefig(OUT + name)
    base = os.path.splitext(name)[0]
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        fig.savefig(EPS + base + '.eps', format='eps')
    bad = [str(x.message) for x in w if 'transparen' in str(x.message).lower()]
    if bad:
        raise RuntimeError(f'{name}: {bad}')
    fig.savefig(EPS + base + '.pdf', format='pdf')
C_FREE, C_EVAP, C_CHIL, C_RESET = '#4C9BE8', '#3BAF7C', '#E8A33D', '#C8553D'
P_REF = 101325.0


def lab(ax, s, x=-0.14, y=1.03):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom')


# ------------------------------------------------------------------ Figure 3: flowchart
def fig3():
    E = json.load(open('numbering.json'))['eq']
    fig, ax = plt.subplots(figsize=(7.2, 6.3)); ax.set_xlim(-2, 114); ax.set_ylim(-4, 96); ax.axis('off')

    def box(x, y, w, h, txt, fc, fs=8.2):
        ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle='round,pad=0.6', fc=fc, ec='#333333', lw=1.0))
        ax.text(x, y, txt, ha='center', va='center', fontsize=fs, linespacing=1.45)
        return y+h/2+0.6, y-h/2-0.6            # outer top and bottom of the rounded frame

    def dia(x, y, w, h, txt):
        ax.add_patch(Polygon([[x-w/2, y], [x, y+h/2], [x+w/2, y], [x, y-h/2]], closed=True, fc='#F2F2F2', ec='#333333', lw=1.0))
        ax.text(x, y, txt, ha='center', va='center', fontsize=8.4)

    def arr(pts, t=None, tp=None):
        xs = [q[0] for q in pts]; ys = [q[1] for q in pts]
        ax.plot(xs, ys, color='#333333', lw=1.0)
        ax.annotate('', xy=pts[-1], xytext=pts[-2], arrowprops=dict(arrowstyle='-|>', lw=1.0, color='#333333', shrinkA=0, shrinkB=0))
        if t:
            ax.text(tp[0], tp[1], t, fontsize=8, ha='center', va='center', color='#333333')
    _, b_top = box(56, 88, 106, 8, 'Weather record: air temperature $T_2$ and relative humidity $RH_2$\n'
                                  'Station pressure $p$, Eq. (%d); wet-bulb temperature $T_{wb2}$ (adiabatic saturation), Eq. (%d)' % (E['patm'], E['twb']), '#FFFFFF')
    dia(22, 68, 36, 13, '$T_2 \\leq T_4 - a_{dry}$ ?')
    dia(63, 68, 36, 13, '$T_{wb2} \\leq T_4 - a_{ev}$ ?')
    ax.text(42.5, 63.2, 'Eq. (%d)' % E['feas'], fontsize=7.6, ha='center', style='italic')
    arr([(22, b_top), (22, 74.5)])
    arr([(40, 68), (45, 68)], 'No', (42.5, 70.4))
    t_f, b_f = box(22, 45, 36, 18, 'Free cooling\n$T_1$ = 27 °C, $T_4$ = 18 °C\n$\\omega_5 = \\omega_2$, $WUE$ = 0\n'
                                   'outlet air $T_5 = T_1 - a_{dry}$', '#DCEBFA')
    t_e, b_e = box(59, 45, 26, 18, 'Evaporative cooling\n$T_1$ = 27 °C\n$T_4$ = 18 °C\n$RH_5$ = 100 %', '#DDF2E7')
    t_c, b_c = box(95, 45, 36, 23, 'Chiller plus evaporative\ncooling\n$T_4$ = max(26 °C, $T_{wb2} + a_{ev}$)\n'
                                   '$T_1 = T_4$ + 24 °C, Eq. (%d)\n$RH_5$ = 100 %%\n$Q_T = P_{ITE}\\,(1 + 1/COP)$' % E['reset'], '#FBEBD3', fs=7.9)
    arr([(81, 68), (95, 68), (95, t_c)], 'No', (88, 70.4))
    arr([(22, 61.5), (22, t_f)], 'Yes', (26.5, 58.3))
    arr([(63, 61.5), (63, t_e)], 'Yes', (67.5, 58.3))
    t_b, _ = box(56, 13, 110, 17, 'Tower or dry cooler: minimum air flow allowed by the counterflow limit with approach $a$, Eq. (%d)\n'
                                  '(air operating line below the saturation enthalpy at $T_w - a_{ev}$; outlet air not warmer than $T_1 - a$)\n'
                                  'Mass and energy balances, Eqs. (%d) to (%d): $WUE$\n'
                                  'Fan and pump power, Eq. (%d), and $PUE$, Eq. (%d); check $\\dot{S}_{gen} \\geq 0$, Eq. (%d)'
                                  % (E['closure'], E['qt'], E['qma5'], E['pow'], E['pue'], E['sgen']), '#FFFFFF', fs=7.4)
    for x, y0 in ((22, b_f), (59, b_e), (95, b_c)):
        arr([(x, y0), (x, t_b)])
    ax.text(56, -2.6, 'The sequence is the minimum of $J = PUE + \\lambda\\,WUE$, Eq. (%d), over the feasible modes for any water value $\\lambda$ in the interval of Section 3.2' % E['J'],
            ha='center', fontsize=7.6, style='italic')
    save(fig, 'fig03_flowchart.png'); plt.close(fig)


# ------------------------------------------------------------------ tower operation (Figs 4 and 5)
def tower_grid(T2=25.0, RH2=30.0, T1=50.0, T4=26.0, p=P_REF):
    prm = dict(R.DEFAULT); pr = R.Props()
    cop = float(R.cop_fun(T4, T1, prm)); Q = 1.0 + 1.0/cop
    T5 = np.linspace(T2+0.05, T1, 300); RH5 = np.linspace(RH2+0.5, 100.0, 300)
    TT, HH = np.meshgrid(T5, RH5)
    pv2 = RH2/100*R.psat(T2); w2 = R.w_from(pv2, p); H2 = pr.hA(T2)+w2*pr.hV(T2)
    pv5 = HH/100*R.psat(TT); w5 = R.w_from(pv5, p)
    hL1, hL4, hL3 = pr.hL(T1), pr.hL(T4), pr.hL(T2)
    mw = Q/(hL1-hL4)
    den = pr.hA(TT)+w5*pr.hV(TT)-H2-(w5-w2)*hL3
    ma = np.where(den > 0, Q/den, np.nan)
    WUE = (w5-w2)*ma*3600
    VA5 = ma*R.R_A*1e3*(R.K+TT)/(p-pv5)
    Ptower = VA5*prm['dp_fan']/1000 + (mw/1000)*prm['dp_pump']/1000
    S = R.entropy_gen(mw, T1, T4, T2, w2, pv2, TT, w5, pv5, ma, T2, p, pr)
    # counterflow admissibility for a = 0 and a = 4: ma >= ma_min(a)
    ok = {}
    pts = {}
    for a in (0.0, 4.0):
        st = R.tower_states(np.array([T1]), np.array([T4]), np.array([T2]), np.array([RH2]), np.array([Q]), p,
                            dict(prm, a_ev=a), pr, True)
        ok[a] = ma >= st['ma'][0]*(1-1e-9)
        wm = float(st['WUE'][0]); pm = float(st['VA5'][0]*prm['dp_fan']/1000 + (st['mw'][0]/1000)*prm['dp_pump']/1000)
        pts[a] = (float(st['T5'][0]), 100.0, wm, pm)
    return TT, HH, WUE, Ptower, S, ok, pts, Q


def _label_positions(cs, target, wx=25.0, wy=70.0):
    """One manual label position per contour level: the vertex of the level closest to target(x, y) in axis-scaled units."""
    pos = []
    xs, ys = target
    for lev, segs in zip(cs.levels, cs.allsegs):
        best = None
        for seg in segs:
            if len(seg) < 2:
                continue
            d = ((seg[:, 0]-xs)/wx)**2 + ((seg[:, 1]-ys)/wy)**2
            i = int(np.argmin(d))
            if best is None or d[i] < best[0]:
                best = (d[i], (float(seg[i, 0]), float(seg[i, 1])))
        if best is not None:
            pos.append(best[1])
    return pos


def _clear_labels(ax, fig, cs, others, fmt, fs, target, placed, bounds, wx=25.0, wy=70.0):
    """For each level of cs, the vertex of that level closest to target whose distance, in display units, to every
    curve of the other family and to every label already placed exceeds the half-length of the label plus a margin."""
    fig.canvas.draw()
    tr = ax.transData
    oth = np.vstack([sg for segs in others for sg in segs if len(sg) > 1])
    # resample the other family densely so that point distances approximate curve distances
    dense = []
    for segs in others:
        for sg in segs:
            if len(sg) < 2:
                continue
            for a_, b_ in zip(sg[:-1], sg[1:]):
                n_ = max(2, int(np.hypot(*(tr.transform(b_) - tr.transform(a_)))/1.5))
                dense.append(np.linspace(a_, b_, n_))
    oth_px = tr.transform(np.vstack(dense))
    x0, x1, y0, y1 = bounds
    pos = []
    px_per_pt = fig.dpi/72.0
    for lev, segs in zip(cs.levels, cs.allsegs):
        txt = fmt(lev); hl = 0.5*len(txt)*fs*0.62*px_per_pt; hh = 0.5*fs*px_per_pt
        need = hl + 3.0*px_per_pt
        best = None
        for sg in segs:
            if len(sg) < 2:
                continue
            for a_, b_ in zip(sg[:-1], sg[1:]):
                for q in np.linspace(a_, b_, 6):
                    if not (x0 < q[0] < x1 and y0 < q[1] < y1):
                        continue
                    qp = tr.transform(q)
                    d_o = np.min(np.hypot(*(oth_px - qp).T))
                    d_l = min([np.hypot(*(qp - pp)) for pp in placed] + [1e9])
                    if d_o < need or d_l < 2*need:
                        continue
                    dt = ((q[0]-target[0])/wx)**2 + ((q[1]-target[1])/wy)**2
                    if best is None or dt < best[0]:
                        best = (dt, (float(q[0]), float(q[1])), qp)
        if best is None:
            raise RuntimeError(f'no clear position for label {txt}')
        pos.append(best[1]); placed.append(best[2])
    return pos


def _check_labels(fig, ax, texts, others):
    """Verify on the rendered figure that no label box is crossed by a curve of the other family."""
    from matplotlib.path import Path
    fig.canvas.draw()
    bad = []
    paths = [Path(ax.transData.transform(sg)) for segs in others for sg in segs if len(sg) > 1]
    for t in texts:
        bb = t.get_window_extent().expanded(1.05, 1.15)
        if any(pth.intersects_bbox(bb, filled=False) for pth in paths):
            bad.append(t.get_text())
    return bad


def fig4_5():
    TT, HH, WUE, Pt, S, ok, pts, Q = tower_grid()
    fig, ax = plt.subplots(figsize=(6.3, 4.6))
    inad = ~ok[0.0] | (S < 0)
    ax.contourf(TT, HH, np.where(inad, 1, np.nan), levels=[0.5, 1.5], colors=['#E6E6E6'])
    ax.contourf(TT, HH, np.where(ok[0.0] & ~ok[4.0], 1, np.nan), levels=[0.5, 1.5], colors=['#F7E9C6'])
    ax.set_xlim(25, 50); ax.set_ylim(30, 100)
    cw = ax.contour(TT, HH, WUE, levels=[1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6], colors='#1F4E79', linewidths=1.1)
    cp = ax.contour(TT, HH, Pt, levels=[0.004, 0.005, 0.0075, 0.01], colors='#A23B2A', linewidths=1.0, linestyles='--')
    fw = lambda v: f'{v:.1f}'; fpp = lambda v: f'{v:.4g}'
    placed = []
    # WUE labels near RH5 of about 62 %, away from the corner (T5 = T2, RH5 = RH2) where all contours meet;
    # every label is placed where no dashed curve passes within its half-length
    wpos = _clear_labels(ax, fig, cw, cp.allsegs, fw, 8.0, (41.0, 62.0), placed, (25.8, 49.2, 32.0, 97.0))
    ppos = _clear_labels(ax, fig, cp, cw.allsegs, fpp, 7.5, (38.0, 88.0), placed, (25.8, 49.2, 32.0, 97.0))
    tw = ax.clabel(cw, fmt=fw, fontsize=8, inline_spacing=4, manual=wpos)
    tp = ax.clabel(cp, fmt=fpp, fontsize=7.5, inline_spacing=4, manual=ppos)
    bad = _check_labels(fig, ax, tw, cp.allsegs) + _check_labels(fig, ax, tp, cw.allsegs)
    assert not bad, f'labels crossed by the other family: {bad}'
    for a, mk in ((0.0, 'o'), (4.0, 's')):
        t5, rh5, w, pwr = pts[a]
        ax.plot(t5, rh5-0.8, mk, ms=7, mfc='k', mec='white', clip_on=False, zorder=5)
    ax.set_xlabel('Outlet air temperature $T_5$ (°C)'); ax.set_ylabel('Outlet relative humidity $RH_5$ (%)')
    ax.set_xlim(25, 50); ax.set_ylim(30, 100)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    (t0, _, w0, p0), (t4, _, w4, p4) = pts[0.0], pts[4.0]
    ax.legend(handles=[Line2D([], [], color='#1F4E79', lw=1.1, label='$WUE$ (L·kWh$_{ITE}^{-1}$)'),
                       Line2D([], [], color='#A23B2A', lw=1.0, ls='--', label='$PUE_{tower}$ (kWh·kWh$_{ITE}^{-1}$)'),
                       Patch(fc='#E6E6E6', label='not admissible (counterflow limit or $\\dot S_{gen}<0$)'),
                       Patch(fc='#F7E9C6', label='admissible only with zero approach'),
                       Patch(fc='white', ec='#BFBFBF', lw=0.6, label='admissible with the design approach of 4 K'),
                       Line2D([], [], ls='', marker='o', mfc='k', mec='white', ms=7,
                              label=f'limit, $a_{{ev}}$ = 0: $T_5$ = {t0:.1f} °C, $WUE$ = {w0:.2f}, $PUE_{{tower}}$ = {p0:.4f}'),
                       Line2D([], [], ls='', marker='s', mfc='k', mec='white', ms=7,
                              label=f'design, $a_{{ev}}$ = 4 K: $T_5$ = {t4:.1f} °C, $WUE$ = {w4:.2f}, $PUE_{{tower}}$ = {p4:.4f}')],
              loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False, fontsize=7.4)
    save(fig, 'fig04_tower_map.png'); plt.close(fig)
    # Figure 5: admissible (PUE_tower, WUE) states, lower frontier, saturated states reachable with more air, iso-cost lines
    fig, ax = plt.subplots(figsize=(6.3, 4.6))
    adm = ok[0.0] & (S >= 0) & np.isfinite(WUE) & (Pt < 0.021)
    ax.scatter(Pt[adm], WUE[adm], s=0.4, color='#C9C9C9', rasterized=True, lw=0)
    xb = np.linspace(0.0036, 0.02, 56); fr = []
    for x0_, x1_ in zip(xb[:-1], xb[1:]):
        m = adm & (Pt >= x0_) & (Pt < x1_)
        fr.append((0.5*(x0_+x1_), np.nanmin(WUE[m]) if m.any() else np.nan))
    fr = np.array(fr)
    ax.plot(fr[:, 0], fr[:, 1], color='#1F4E79', lw=2.0)
    j = int(np.argmin(np.abs(HH[:, 0]-100.0))); row = ok[4.0][j] & (S[j] >= 0) & np.isfinite(WUE[j])
    ax.plot(Pt[j][row], WUE[j][row], color='#2CA02C', lw=2.0)
    for rh, ls in ((80, '--'), (60, ':')):
        jj = int(np.argmin(np.abs(HH[:, 0]-rh))); m = adm[jj]
        ax.plot(Pt[jj][m], WUE[jj][m], color='#555555', lw=1.0, ls=ls)
    (t0, _, w0, p0), (t4, _, w4, p4) = pts[0.0], pts[4.0]
    ax.plot(p0, w0, 'o', ms=8, mfc='white', mec='k', mew=1.3, zorder=5)
    ax.plot(p4, w4, 's', ms=6, mfc='k', mec='k', zorder=6)
    for lam, ls in ((0.02, '-'), (0.2, '--')):
        xs = np.array([p4, 0.02]); ax.plot(xs, w4-(xs-p4)/lam, color='#D62728', lw=1.2, ls=ls)
    ax.set_xlim(0.003, 0.02); ax.set_ylim(0.8, 1.75)
    ax.set_xlabel('Fan and pump energy $PUE_{tower}$ (kWh·kWh$_{ITE}^{-1}$)'); ax.set_ylabel('$WUE$ (L·kWh$_{ITE}^{-1}$)')
    from matplotlib.lines import Line2D
    h = [Line2D([], [], ls='', marker='s', mfc='#C9C9C9', mec='#C9C9C9', ms=6, label='admissible outlet states ($a_{ev}$ = 0)'),
         Line2D([], [], color='#1F4E79', lw=2.0, label='lower frontier (minimum $WUE$, unsaturated outlet)'),
         Line2D([], [], color='#2CA02C', lw=2.0, label='saturated outlet, air flow above the minimum ($a_{ev}$ = 4 K)'),
         Line2D([], [], color='#555555', lw=1.0, ls='--', label='$RH_5$ = 80 %'),
         Line2D([], [], color='#555555', lw=1.0, ls=':', label='$RH_5$ = 60 %'),
         Line2D([], [], ls='', marker='o', mfc='white', mec='k', mew=1.3, ms=8, label='minimum air flow, limit ($a_{ev}$ = 0)'),
         Line2D([], [], ls='', marker='s', mfc='k', mec='k', ms=6, label='minimum air flow, design ($a_{ev}$ = 4 K)'),
         Line2D([], [], color='#D62728', lw=1.2, label='iso-cost line through the design state, $\\lambda$ = 0.02 kWh·L$^{-1}$'),
         Line2D([], [], color='#D62728', lw=1.2, ls='--', label='iso-cost line, $\\lambda$ = 0.2 kWh·L$^{-1}$')]
    ax.legend(handles=h, loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False, fontsize=7.4)
    save(fig, 'fig05_tradeoff.png'); plt.close(fig)
    return pts, Q


# ------------------------------------------------------------------ Figure 6: mode maps over (T2, RH2)
LV_P = [1.00, 1.01, 1.02, 1.03, 1.05, 1.10, 1.15, 1.16, 1.17]
LV_W = [0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.3]


def _cmaps():
    cp = plt.get_cmap('YlOrRd', len(LV_P)-1).copy(); cp.set_over('#3F007D')
    cw = plt.get_cmap('Blues', len(LV_W)-1).copy(); cw.set_under('#E3E3E3'); cw.set_over('#08306B')
    return cp, cw


def fig6():
    Tg = np.linspace(-10, 46, 281); Hg = np.linspace(2, 100, 197)
    TT, HH = np.meshgrid(Tg, Hg)
    r = R.run(TT.ravel(), HH.ravel(), p=P_REF)
    PUE = r['PUE'].reshape(TT.shape); WUE = r['WUE'].reshape(TT.shape); Twb = r['Twb'].reshape(TT.shape)
    unobs = Twb > 30.0
    PUE = np.where(unobs, np.nan, PUE); WUE = np.where(unobs, np.nan, WUE)
    cmp_, cmw = _cmaps()
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 3.9), sharey=True)
    fig.subplots_adjust(wspace=0.1)
    cf = axs[0].contourf(TT, HH, PUE, levels=LV_P, cmap=cmp_, norm=BoundaryNorm(LV_P, cmp_.N), extend='max')
    cb = fig.colorbar(cf, ax=axs[0], orientation='horizontal', pad=0.2, fraction=0.06, aspect=28, ticks=LV_P)
    cb.set_label('$PUE$ (kWh·kWh$_{ITE}^{-1}$)', fontsize=9); cb.ax.tick_params(labelsize=8.5)
    cb.ax.set_xticklabels([f'{v:.2f}' for v in LV_P], rotation=45)
    cf2 = axs[1].contourf(TT, HH, WUE, levels=LV_W, cmap=cmw, norm=BoundaryNorm(LV_W, cmw.N), extend='both')
    cb2 = fig.colorbar(cf2, ax=axs[1], orientation='horizontal', pad=0.2, fraction=0.06, aspect=28, ticks=LV_W)
    cb2.set_label('$WUE$ (L·kWh$_{ITE}^{-1}$); gray: $WUE$ = 0', fontsize=9); cb2.ax.tick_params(labelsize=8.5)
    cb2.ax.set_xticklabels([f'{v:.1f}' for v in LV_W])
    box = dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=1.0)
    arw = dict(arrowstyle='-|>', lw=0.8, color='k', shrinkA=1, shrinkB=0)
    for ax in axs:
        ax.axvline(12.0, color='k', lw=1.2)
        ax.contour(TT, HH, Twb, levels=[14.0, 22.0], colors='k', linewidths=1.2, linestyles=['-', '--'])
        ax.axvline(18.0, color='k', lw=0.8, ls=':')
        ax.contour(TT, HH, Twb, levels=[18.0, 26.0], colors='k', linewidths=0.8, linestyles=':')
        ax.set_xlabel('Air temperature $T_2$ (°C)'); ax.set_xlim(-10, 46); ax.set_ylim(2, 100)
        ax.text(0, 60, 'free', fontsize=8.5, ha='center', va='center', fontweight='bold', bbox=box, zorder=6)
        ax.text(15.0, 40, 'evaporative', fontsize=8.5, ha='center', va='center', rotation=90, fontweight='bold', bbox=box, zorder=6)
        # the chilling and raised-set-point bands are crossed by border lines everywhere: their labels are placed in
        # the blank area (wet-bulb temperature above 30 C) with an arrow into each band
        ax.annotate('chilling', xy=(24.5, 40), xytext=(40.8, 73), fontsize=8.5, ha='center', va='center', fontweight='bold',
                    bbox=box, arrowprops=arw, zorder=6)
        ax.annotate('raised $T_4$', xy=(30.5, 78), xytext=(40.0, 91), fontsize=8.5, ha='center', va='center', style='italic',
                    bbox=box, arrowprops=arw, zorder=6)
    axs[0].set_ylabel('Relative humidity $RH_2$ (%)')
    lab(axs[0], '(a)', -0.16); lab(axs[1], '(b)', -0.06)
    save(fig, 'fig06_mode_maps.png'); plt.close(fig)


# ------------------------------------------------------------------ Figure 7: Pareto front (water value lambda)
def fig7():
    import resq as Q
    sens = json.load(open('sens_out.json'))
    P = sens['pareto']
    keys = sorted(P.keys(), key=lambda k: float(k))
    lo_s, hi_s = Q.lam_fmt()
    fig, ax = plt.subplots(figsize=(6.3, 4.4))
    sts = [('Madrid', '#1F77B4', '-'), ('El Paso', '#9467BD', '--'), ('Lubbock', '#8C564B', '-.'), ('Dallas', '#FF7F0E', ':'), ('Houston', '#D62728', (0, (5, 1, 1, 1, 1, 1)))]
    for st, c, ls in sts:
        xs = [P[k][st]['WUE'] for k in keys]; ys = [P[k][st]['PUE'] for k in keys]
        ax.plot(xs, ys, color=c, lw=1.4, ls=ls, marker='.', ms=4, zorder=2)
        b = sens['cases']['baseline'][st]
        ax.plot(b['WUE'], b['PUE'], 'o', ms=7, mfc=c, mec='k', zorder=4)
        k0 = keys[0]; ax.plot(P[k0][st]['WUE'], P[k0][st]['PUE'], 's', ms=6, mfc='white', mec=c, mew=1.3, zorder=5)
        kd = [k for k in keys if abs(float(k)-1/370) < 1e-9][0]
        ax.plot(P[kd][st]['WUE'], P[kd][st]['PUE'], '^', ms=4.5, mfc=c, mec=c, zorder=6)
        km = keys[-1]; ax.plot(P[km][st]['WUE'], P[km][st]['PUE'], 'D', ms=5.5, mfc='white', mec=c, mew=1.3, zorder=5)
        ax.annotate(st, (P[k0][st]['WUE'], P[k0][st]['PUE']), xytext=(6, 4), textcoords='offset points', fontsize=8, color=c)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=c, lw=1.4, ls=ls, label=s) for s, c, ls in sts]
    h += [Line2D([], [], ls='', marker='o', mfc='grey', mec='k', ms=7, label=f'selection rule ({lo_s} ≤ $\\lambda$ ≤ {hi_s} kWh·L$^{{-1}}$)'),
          Line2D([], [], ls='', marker='s', mfc='white', mec='grey', mew=1.3, ms=6, label='minimum energy, $\\lambda$ = 0'),
          Line2D([], [], ls='', marker='^', mfc='grey', mec='grey', ms=4.5, label='$\\lambda$ = 1/370 kWh·L$^{-1}$ (desalination)'),
          Line2D([], [], ls='', marker='D', mfc='white', mec='grey', mew=1.3, ms=5.5, label='minimum water, $\\lambda$ = 1,000 kWh·L$^{-1}$'),
          Line2D([], [], ls='', marker='.', color='grey', ms=4, label='other values of $\\lambda$ computed')]
    ax.legend(handles=h, loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=3, frameon=False, fontsize=7.4)
    ax.set_xlabel('Long-term mean $WUE$ (L·kWh$_{ITE}^{-1}$)'); ax.set_ylabel('Long-term mean $PUE$ (kWh·kWh$_{ITE}^{-1}$)')
    save(fig, 'fig07_pareto.png'); plt.close(fig)


# ------------------------------------------------------------------ Figures 8 to 11: series of M01 and Dallas
def series_figs():
    d = pd.read_pickle('res_design.pkl')
    sel = [(1, 'Madrid M01'), ('Dallas', 'Dallas')]
    # Figure 8
    fig, axs = plt.subplots(2, 2, figsize=(6.6, 5.2), sharex=True, sharey=True)
    bp, bw = LV_P, LV_W
    cm_p, cm_w = _cmaps()
    for j, (st, nm) in enumerate(sel):
        g = d[d.St == st]
        for i, (y, b, cm) in enumerate((('PUE', bp, cm_p), ('WUE', bw, cm_w))):
            ax = axs[i, j]
            ext = 'max' if y == 'PUE' else 'both'
            norm = BoundaryNorm(b, cm.N)
            o = np.argsort(g[y].values)
            sc = ax.scatter(g['T'].values[o], g['RH'].values[o], c=g[y].values[o], s=0.6, cmap=cm, norm=norm, rasterized=True, lw=0)
            ax.set_title(f'${y}$, {nm}', fontsize=9)
            lab(ax, '(' + 'abcd'[2*i+j] + ')', -0.16 if j == 0 else -0.06)
            if j == 1:
                cb = fig.colorbar(sc, ax=axs[i, :], pad=0.015, fraction=0.03, ticks=b, extend=ext)
                cb.set_label(f'${y}$ ' + ('(kWh·kWh$_{ITE}^{-1}$)' if y == 'PUE' else '(L·kWh$_{ITE}^{-1}$); gray: 0'), fontsize=9)
                cb.ax.tick_params(labelsize=8.5)
                cb.ax.set_yticklabels([f'{v:.2f}' if y == 'PUE' else f'{v:.1f}' for v in b])
            ax.set_xlim(-12, 46); ax.set_ylim(0, 100)
            ax.axvline(12, color='k', lw=0.6)
    for ax in axs[1]:
        ax.set_xlabel('$T_2$ (°C)')
    for ax in axs[:, 0]:
        ax.set_ylabel('$RH_2$ (%)')
    save(fig, 'fig08_series_maps.png'); plt.close(fig)
    # Figure 9: distributions
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.9))
    fig.subplots_adjust(wspace=0.32)
    for j, y in enumerate(('PUE', 'WUE')):
        ax = axs[j]
        for (st, nm), c, ls, hat in zip(sel, ('#1F77B4', '#FF7F0E'), ('-', '--'), (None, '////')):
            g = d[d.St == st][y].values
            bins = np.linspace(1.0, 1.17, 69) if y == 'PUE' else np.linspace(0, 2.3, 47)
            ax.hist(g, bins=bins, density=True, histtype='stepfilled', alpha=0.30, color=c, hatch=hat, rasterized=True)
            ax.hist(g, bins=bins, density=True, histtype='step', color=c, lw=1.0, ls=ls)
            ax.axvline(g.mean(), color=c, ls=':', lw=1.3)
        ax.set_xlabel(f'${y}$ ' + ('(kWh·kWh$_{ITE}^{-1}$)' if y == 'PUE' else '(L·kWh$_{ITE}^{-1}$)'))
        ax.set_ylabel('Probability density')
        if j == 0:
            from matplotlib.patches import Patch
            from matplotlib.colors import to_rgb
            mix = lambda c: tuple(0.30*v + 0.70 for v in to_rgb(c))
            ax.legend(handles=[Patch(fc=mix(c), ec=c, hatch=hat, lw=1.0, ls=ls, label=nm)
                               for (st, nm), c, ls, hat in zip(sel, ('#1F77B4', '#FF7F0E'), ('-', '--'), (None, '////'))],
                      frameon=False, loc='upper center', bbox_to_anchor=(1.12, -0.3), ncol=2)
        if y == 'PUE':
            ax.set_yscale('log')
    lab(axs[0], '(a)'); lab(axs[1], '(b)')
    save(fig, 'fig09_distributions.png'); plt.close(fig)
    # Figures 10 and 11: intraday evolution on the equinoxes and solstices
    days = [(3, 21, 'Spring equinox'), (6, 21, 'Summer solstice'), (9, 21, 'Autumn equinox'), (12, 21, 'Winter solstice')]
    for y, fn, ylim in (('PUE', 'fig10_intraday_pue.png', (0.99, 1.18)), ('WUE', 'fig11_intraday_wue.png', (-0.05, 2.2))):
        fig, axs = plt.subplots(4, 2, figsize=(6.6, 7.6), sharex=True, sharey=True)
        for j, (st, nm) in enumerate(sel):
            g = d[d.St == st]
            hh = g['ts'].dt.hour.values + g['ts'].dt.minute.values/60.0
            for i, (mo, dy, dn) in enumerate(days):
                m = (g['ts'].dt.month.values == mo) & (g['ts'].dt.day.values == dy)
                ax = axs[i, j]
                sr = pd.Series(g[y].values[m]).groupby(hh[m])
                q = sr.quantile([0.1, 0.9]).unstack(); mm = sr.mean()
                ax.fill_between(q.index, q[0.1], q[0.9], color='#BDBDBD', lw=0, zorder=1)
                ax.plot(mm.index, mm.values, color='#D62728', lw=1.6, zorder=3)
                ax.set_title(f'{dn}, {nm}', fontsize=8.5)
                ax.set_xlim(0, 24); ax.set_ylim(*ylim); ax.set_xticks([0, 6, 12, 18, 24])
        for ax in axs[-1]:
            ax.set_xlabel('Time of day (h)')
        for ax in axs[:, 0]:
            ax.set_ylabel(f'${y}$ ' + ('(kWh·kWh$_{ITE}^{-1}$)' if y == 'PUE' else '(L·kWh$_{ITE}^{-1}$)'), fontsize=8)
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
        fig.legend(handles=[Line2D([], [], color='#D62728', lw=1.6, label='mean over the years'),
                            Patch(fc='#BDBDBD', label='10th to 90th percentiles of the years')],
                   loc='lower center', bbox_to_anchor=(0.5, -0.02), ncol=2, frameon=False, fontsize=8)
        fig.subplots_adjust(bottom=0.09)
        save(fig, fn); plt.close(fig)


# ------------------------------------------------------------------ Figure 12: monthly water vs ET0
def fig12():
    A = json.load(open('analysis_out.json'))['annual']
    fig, axs = plt.subplots(1, 5, figsize=(7.0, 2.5), sharey=True)
    mon = np.arange(1, 13)
    for ax, st in zip(axs, ['Madrid', 'El Paso', 'Lubbock', 'Dallas', 'Houston']):
        a = A[st]
        ax.axvspan(3.5, 9.5, color='#EEF5E9', zorder=0)
        ax.bar(mon, np.array(a['W_m']), color='#4C9BE8', width=0.75)
        ax2 = ax.twinx(); ax2.plot(mon, a['ET_m'], color='#2CA02C', marker='o', ms=2.5, lw=1.2)
        ax2.set_ylim(0, 300)
        if st != 'Houston':
            ax2.set_yticks([]); ax2.spines['right'].set_visible(False)
        else:
            ax2.set_ylabel('ET$_0$ (mm·month$^{-1}$)', color='#2CA02C')
        ax2.spines['top'].set_visible(False)
        ax.set_title(f"{st}\n$r$ = {a['r']:.2f}", fontsize=8.5)
        ax.set_xticks([1, 4, 7, 10]); ax.set_xticklabels(['Jan', 'Apr', 'Jul', 'Oct'], fontsize=7.5)
    axs[0].set_ylabel('Water demand (ML·month$^{-1}$)')
    save(fig, 'fig12_water_calendar.png'); plt.close(fig)


if __name__ == '__main__':
    which = sys.argv[1:] or ['3', '45', '6', '7', 'series', '12']
    if '3' in which: fig3()
    if '45' in which: print(fig4_5())
    if '6' in which: fig6()
    if '7' in which: fig7()
    if 'series' in which: series_figs()
    if '12' in which: fig12()
    print('figures done')
