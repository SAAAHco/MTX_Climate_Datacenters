# Figure 2: heat transport from the data center to the atmosphere through free cooling, evaporative cooling, or chiller
# plus evaporative cooling, depending on the position of the 3-way valves. Redrawn from the authors' schematic with the
# same components, pipe connections, and valve ports (the side port of every 3-way valve is its common port, as in the
# original), at a size that keeps every label legible.
import os, warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Polygon, Circle, Ellipse, PathPatch, Wedge
from matplotlib.path import Path

HOT, COLD, SPRAY, INNER, AIR = '#E0312B', '#1F4FD1', '#12B5CB', '#C21BC2', '#9ACD32'
plt.rcParams.update({'ps.fonttype': 42, 'pdf.fonttype': 42})
LW = 1.5
FS = 8.5


def fig2(out='figs/fig02_schematic.png'):
    fig, ax = plt.subplots(figsize=(7.0, 6.8))
    ax.set_xlim(0, 102); ax.set_ylim(-1, 96); ax.set_aspect('equal'); ax.axis('off')

    def pipe(pts, c, lw=LW, z=2, ls='-'):
        xs, ys = zip(*pts); ax.plot(xs, ys, color=c, lw=lw, zorder=z, ls=ls, solid_capstyle='butt')

    def head(x, y, d, c, s=0.8, z=3):
        dx, dy = {'r': (1, 0), 'l': (-1, 0), 'u': (0, 1), 'd': (0, -1)}[d]
        px, py = -dy, dx
        ax.add_patch(Polygon([(x + dx*s, y + dy*s), (x - dx*s*0.4 + px*s*0.6, y - dy*s*0.4 + py*s*0.6),
                              (x - dx*s*0.4 - px*s*0.6, y - dy*s*0.4 - py*s*0.6)], closed=True, fc=c, ec=c, lw=0.5, zorder=z))

    def pump(x, y, d, r=1.4):
        ax.add_patch(Circle((x, y), r, fc='white', ec='k', lw=0.9, ls=(0, (2, 1.2)), zorder=4))
        dx, dy = {'r': (1, 0), 'l': (-1, 0), 'u': (0, 1), 'd': (0, -1)}[d]; px, py = -dy, dx
        ax.add_patch(Polygon([(x + dx*r*0.7, y + dy*r*0.7), (x - dx*r*0.45 + px*r*0.55, y - dy*r*0.45 + py*r*0.55),
                              (x - dx*r*0.45 - px*r*0.55, y - dy*r*0.45 - py*r*0.55)], closed=True, fc='k', ec='k', lw=0.4, zorder=5))

    def valve3(x, y, main='v', stub='l', s=1.3):
        # two triangles along the main pipe and a smaller third triangle on the side port (the common port)
        k, ks = 0.55*s, 0.45*s
        if main == 'v':
            tri = [[(x-k, y+s), (x+k, y+s), (x, y)], [(x-k, y-s), (x+k, y-s), (x, y)]]
        else:
            tri = [[(x-s, y-k), (x-s, y+k), (x, y)], [(x+s, y-k), (x+s, y+k), (x, y)]]
        st = {'r': [(x+s, y-ks), (x+s, y+ks), (x, y)], 'l': [(x-s, y-ks), (x-s, y+ks), (x, y)],
              'u': [(x-ks, y+s), (x+ks, y+s), (x, y)], 'd': [(x-ks, y-s), (x+ks, y-s), (x, y)]}[stub]
        for t in tri + [st]:
            ax.add_patch(Polygon(t, closed=True, fc='white', ec='k', lw=0.9, zorder=4))

    def needle(x, y, s=1.3):
        k = 0.7*s
        for t in ([(x-s, y-k), (x-s, y+k), (x, y)], [(x+s, y-k), (x+s, y+k), (x, y)]):
            ax.add_patch(Polygon(t, closed=True, fc='white', ec='k', lw=0.9, zorder=4))
        ax.add_patch(Circle((x, y), 0.28, fc='k', zorder=5))

    def plate_hx(x0, y0, x1, y1, fw=2.2, lab=None, lab_side='r'):
        # end frames on the left and right, plates parallel to them
        ax.add_patch(Rectangle((x0, y0), fw, y1-y0, fc='white', ec='k', lw=1.0, zorder=3))
        ax.add_patch(Rectangle((x1-fw, y0), fw, y1-y0, fc='white', ec='k', lw=1.0, zorder=3))
        for yy in (y0+0.7, y1-0.7):
            ax.plot([x0+fw, x1-fw], [yy, yy], color='k', lw=0.8, zorder=3)
        for xx in np.linspace(x0+fw+1.0, x1-fw-1.0, 6):
            ax.plot([xx, xx], [y0+0.7, y1-0.7], color='k', lw=0.6, zorder=3)
        if lab:
            xl = x1-fw/2 if lab_side == 'r' else x0+fw/2
            ax.text(xl, (y0+y1)/2, lab, rotation=90, ha='center', va='center', fontsize=FS-1.5, zorder=6)

    def plate_hx_h(x0, y0, x1, y1, fw=2.0):
        # end frames at the top and bottom, plates parallel to them
        ax.add_patch(Rectangle((x0, y1-fw), x1-x0, fw, fc='white', ec='k', lw=1.0, zorder=3))
        ax.add_patch(Rectangle((x0, y0), x1-x0, fw, fc='white', ec='k', lw=1.0, zorder=3))
        for xx in (x0+0.7, x1-0.7):
            ax.plot([xx, xx], [y0+fw, y1-fw], color='k', lw=0.8, zorder=3)
        for yy in np.linspace(y0+fw+0.8, y1-fw-0.8, 5):
            ax.plot([x0+0.7, x1-0.7], [yy, yy], color='k', lw=0.6, zorder=3)

    def curve(p0, c1, c2, p1, c, d):
        pth = Path([p0, c1, c2, p1], [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
        ax.add_patch(PathPatch(pth, fc='none', ec=c, lw=1.2, ls=(0, (3, 2)), zorder=5))
        head(p1[0], p1[1], d, c, s=0.8, z=6)

    def fan(x, y, w=7.0, h=1.6, z=4):
        ax.add_patch(Ellipse((x-w/4, y), w/2, h, fc='white', ec='k', lw=0.9, zorder=z))
        ax.add_patch(Ellipse((x+w/4, y), w/2, h, fc='white', ec='k', lw=0.9, zorder=z))
        ax.add_patch(Circle((x, y), 0.35, fc='k', zorder=z+1))

    def motor(x, y, w=4.0, h=3.2):
        ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle='round,pad=0.25', fc='white', ec='k', lw=0.9, zorder=4))
        for yy in (y-h/2+0.5, y+h/2-0.5):
            ax.plot([x-w/2, x+w/2], [yy, yy], color='k', lw=0.6, zorder=5)

    def dbox(x0, y0, x1, y1):
        ax.add_patch(Rectangle((x0, y0), x1-x0, y1-y0, fc='none', ec='#666666', lw=1.1, ls=(0, (4, 2.5)), zorder=1))

    def air(pts):
        xs, ys = np.array(pts).T
        t = np.linspace(0, 1, len(xs)); tt = np.linspace(0, 1, 80)
        cx = np.polyfit(t, xs, min(3, len(xs)-1)); cy = np.polyfit(t, ys, min(3, len(xs)-1))
        X = np.polyval(cx, tt); Y = np.polyval(cy, tt)
        ax.plot(X[:-3], Y[:-3], color=AIR, lw=3.0, ls=(0, (0.6, 0.6)), zorder=6, solid_capstyle='butt')
        dx, dy = X[-1]-X[-4], Y[-1]-Y[-4]; n = np.hypot(dx, dy); dx, dy = dx/n, dy/n; px, py = -dy, dx
        ax.add_patch(Polygon([(X[-1]+dx*1.2, Y[-1]+dy*1.2), (X[-4]+px*1.3, Y[-4]+py*1.3), (X[-4]-px*1.3, Y[-4]-py*1.3)],
                             closed=True, fc=AIR, ec=AIR, zorder=6))
        return X[-1], Y[-1]

    # ---------------- key (symbols, then line styles)
    R1, R2, R3 = 93.0, 89.2, 85.4
    pump(3.5, R1, 'r', r=1.3); ax.text(6.0, R1, 'Pump / compressor', va='center', fontsize=FS)
    needle(3.5, R2); ax.text(6.0, R2, 'Needle valve', va='center', fontsize=FS)
    valve3(3.5, R3, 'v', 'l'); ax.text(6.0, R3, '3-way valve', va='center', fontsize=FS)
    plate_hx(30.5, R1-1.4, 36.5, R1+1.4, fw=0.9); ax.text(38.0, R1, 'Heat exchanger', va='center', fontsize=FS)
    fan(33.5, R2, w=5.5, h=1.3); ax.text(38.0, R2, 'Fan', va='center', fontsize=FS)
    motor(33.5, R3, w=3.0, h=2.2); ax.text(38.0, R3, 'Fan motor', va='center', fontsize=FS)
    pipe([(59.5, R1), (63.0, R1)], HOT); ax.text(64.5, R1, 'hot fluid', va='center', fontsize=FS)
    pipe([(59.5, R2), (63.0, R2)], COLD); ax.text(64.5, R2, 'cold fluid', va='center', fontsize=FS)
    pipe([(59.5, R3+0.6), (63.0, R3+0.6)], SPRAY, lw=1.2)
    ax.plot([59.5, 63.0], [R3-0.7, R3-0.7], color=SPRAY, lw=2.2, ls=(0, (1.6, 1.0))); ax.text(64.5, R3, 'spray water', va='center', fontsize=FS)
    ax.plot([80.5, 84.0], [R1, R1], color=INNER, lw=1.2, ls=(0, (3, 2))); ax.text(85.5, R1, 'flow in heat', va='center', fontsize=FS)
    ax.text(85.5, R1-2.3, 'exchanger', va='center', fontsize=FS)
    ax.plot([80.5, 84.0], [R3, R3], color=AIR, lw=3.0, ls=(0, (0.6, 0.6))); ax.text(85.5, R3, 'air', va='center', fontsize=FS)
    TT = 80.9      # box titles, above the boxes as in the original

    # ---------------- ITE side (secondary loop) and its heat exchanger
    dbox(1.5, 55.5, 24.5, 79.2)
    ax.add_patch(Rectangle((3.0, 73.0), 20.0, 4.6, fc='white', ec='k', lw=1.1, zorder=3))
    ax.text(13.0, 75.3, 'ITE heat load', ha='center', va='center', fontsize=FS+0.5, fontweight='bold', zorder=4)
    pipe([(7, 73.0), (7, 68.0)], HOT); pump(7, 71.2, 'd', r=1.0); head(7, 69.0, 'd', HOT, 0.7)
    pipe([(19, 68.0), (19, 73.0)], COLD); head(19, 71.2, 'u', COLD, 0.7)
    plate_hx_h(5.0, 59.0, 21.0, 68.0, fw=2.0)
    curve((8.0, 65.6), (9.5, 64.3), (16.5, 64.3), (18.0, 65.6), INNER, 'u')
    curve((18.0, 61.4), (16.5, 62.7), (9.5, 62.7), (8.0, 61.4), INNER, 'd')
    ax.text(7.6, 53.3, '27 °C', ha='left', va='center', fontsize=FS)      # between the two pipes, each beside its own pipe
    ax.text(18.4, 50.3, '18 °C', ha='right', va='center', fontsize=FS)

    # ---------------- free cooling (dry cooler)
    dbox(28.0, 44.5, 60.0, 79.2)
    ax.text(44.0, TT, 'Free cooling', ha='center', fontsize=FS+1, fontweight='bold', va='center')
    body = [(33.6, 58.0), (36.0, 56.0), (53.5, 56.0), (56.0, 58.0), (56.0, 67.0), (53.5, 69.0), (36.0, 69.0), (33.6, 67.0)]
    ax.add_patch(Polygon(body, closed=True, fc='white', ec='k', lw=1.1, zorder=3))
    ax.add_patch(Rectangle((34.2, 61.6), 21.2, 1.8, fc='white', ec='k', lw=0.8, hatch='////', zorder=3))
    ax.add_patch(Rectangle((34.2, 59.8), 21.2, 1.8, fc='white', ec='k', lw=0.8, zorder=3))
    ax.plot([56.0, 57.2], [62.5, 62.5], color='k', lw=0.8, zorder=3)
    ax.add_patch(Wedge((57.2, 62.5), 1.1, -90, 90, fc='white', ec='k', lw=0.9, zorder=4))
    fan(45.0, 66.4); fan(45.0, 57.9, h=1.4)
    ax.add_patch(Rectangle((31.2, 58.2), 2.4, 8.6, fc='white', ec='k', lw=1.0, zorder=3))
    ax.plot([45, 45], [56.0, 53.4], color='k', lw=0.9, zorder=3); motor(45.0, 51.8, w=4.0, h=2.8)
    ax1 = air([(49.8, 47.5), (49.0, 58.0), (48.0, 66.0), (48.8, 75.6)])
    ax.text(ax1[0]+1.6, ax1[1]-0.6, 'Air', fontsize=FS, va='center')

    # ---------------- evaporative cooling (closed-circuit tower)
    dbox(62.5, 44.5, 101.0, 79.2)
    ax.text(81.5, TT, 'Evaporative cooling', ha='center', fontsize=FS+1, fontweight='bold', va='center')
    # casing with an open fan stack and an inward step to the basin; the right wall is open only at the air inlet
    ax.plot([70.0, 68.0, 68.0, 75.5, 75.5], [50.5, 50.5, 64.0, 68.8, 71.2], color='k', lw=1.1, zorder=3)
    ax.plot([86.5, 86.5, 94.0, 94.0], [71.2, 68.8, 64.0, 58.4], color='k', lw=1.1, zorder=3)
    ax.plot([92.0, 94.0, 94.0], [50.5, 50.5, 54.2], color='k', lw=1.1, zorder=3)
    ax.add_patch(Rectangle((70.0, 47.8), 22.0, 2.7, fc=SPRAY, ec='k', lw=1.0, zorder=3))
    ax.text(81.0, 49.15, 'Water', ha='center', va='center', fontsize=FS-0.5, zorder=4)
    fan(81.0, 70.0, w=10, h=1.7); ax.plot([81, 81], [70.85, 71.9], color='k', lw=0.9, zorder=3); motor(81.0, 73.4, w=3.6, h=2.8)
    ax.plot([71.0, 91.0], [62.0, 62.0], color='k', lw=1.3, zorder=4)
    for xn in (76.0, 86.0):
        ax.plot([xn-1.6, xn, xn+1.6], [60.4, 62.0, 60.4], color='k', lw=1.1, zorder=4)
        ax.plot([xn, xn], [62.0, 60.4], color='k', lw=1.1, zorder=4)
        for dx_ in (-2.8, 0.0, 2.8):          # spray falls over the coil into the basin
            ax.plot([xn+dx_*0.4, xn+dx_, xn+dx_], [60.2, 58.4, 50.6], color=SPRAY, lw=2.2, ls=(0, (1.6, 1.0)), zorder=3)
    coil = [(93.0, 53.0), (71.2, 53.0), (71.2, 55.0), (90.8, 55.0), (90.8, 57.0), (69.3, 57.0)]
    pipe(coil, 'k', lw=1.3, z=4)
    pipe([(69.3, 57.0), (67.0, 57.0)], COLD, z=4)
    pipe([(70.0, 49.15), (63.6, 49.15), (63.6, 62.0), (71.0, 62.0)], SPRAY, lw=1.2); pump(64.9, 49.15, 'l', r=1.0); head(69.9, 62.0, 'r', SPRAY, 0.8)
    ax2 = air([(100.2, 56.0), (92.0, 56.1), (82.0, 60.8), (86.6, 75.8)])
    ax.text(ax2[0]+1.8, ax2[1]-0.4, 'Air', fontsize=FS, va='center')

    # ---------------- chiller (mechanical cooling)
    dbox(9.5, 4.5, 72.5, 26.5)
    ax.text(41.0, 1.6, 'Chiller (mechanical cooling)', ha='center', va='center', fontsize=FS+1, fontweight='bold')
    plate_hx(18.0, 8.0, 32.0, 22.0, fw=2.4, lab='Evaporator', lab_side='r')
    plate_hx(50.0, 8.0, 64.0, 22.0, fw=2.4, lab='Condenser', lab_side='l')
    curve((20.2, 9.3), (24.5, 11.5), (24.5, 18.5), (20.2, 20.7), INNER, 'u')
    curve((28.9, 20.7), (25.3, 18.5), (25.3, 11.5), (28.9, 9.3), INNER, 'd')
    curve((53.1, 9.3), (56.6, 11.5), (56.6, 18.5), (53.1, 20.7), INNER, 'u')
    curve((61.8, 20.7), (57.6, 18.5), (57.6, 11.5), (61.8, 9.3), INNER, 'd')
    pipe([(50.0, 20.0), (32.0, 20.0)], COLD); needle(41.0, 20.0); head(35.0, 20.0, 'l', COLD)
    pipe([(32.0, 10.0), (50.0, 10.0)], HOT); pump(41.0, 10.0, 'r'); head(47.2, 10.0, 'r', HOT)
    ax.text(41.0, 15.2, 'Vapor\ncompression\ncycle', ha='center', va='center', fontsize=FS-1.5, linespacing=1.15)

    # ---------------- primary circuit. As in the original, the side port of every 3-way valve is its common port:
    # V1 and V2 divide the hot water (inlet on the side port), V3 and V4 join the cold water (outlet on the side port),
    # V5 divides the water leaving the tower coil, and V6 joins the two lines that feed it.
    s = 1.3
    V1, V2, V3, V4, V5, V6 = (9.5, 45.0), (12.0, 31.5), (21.5, 51.0), (24.0, 39.0), (67.0, 41.5), (97.0, 34.5)
    # hot 27 C water from the ITE heat exchanger into the side port of V1
    pipe([(7.0, 59.0), (7.0, V1[1]), (V1[0]-s, V1[1])], HOT); head(7.0, 51.5, 'd', HOT)
    # V1 top port: free cooling (pump, dry-cooler header)
    pipe([(V1[0], V1[1]+s), (V1[0], 47.5), (32.4, 47.5), (32.4, 58.2)], HOT); pump(30.0, 47.5, 'r', r=1.3); head(32.4, 54.5, 'u', HOT)
    # V1 bottom port into the side port of V2
    pipe([(V1[0], V1[1]-s), (V1[0], V2[1]), (V2[0]-s, V2[1])], HOT); head(V1[0], 37.5, 'd', HOT)
    # V2 top port: evaporative cooling (pump, V6, tower coil)
    pipe([(V2[0], V2[1]+s), (V2[0], V6[1]), (V6[0]-s, V6[1])], HOT); pump(45.0, V6[1], 'r'); head(88.0, V6[1], 'r', HOT)
    # V2 bottom port: chiller (pump, evaporator)
    pipe([(V2[0], V2[1]-s), (V2[0], 10.0), (18.0, 10.0)], HOT); pump(15.2, 10.0, 'r', r=1.1); head(V2[0], 20.0, 'd', HOT)
    # condenser outlet (pump) into the right port of V6; V6 side port up to the tower coil
    pipe([(64.0, 10.0), (100.2, 10.0), (100.2, V6[1]), (V6[0]+s, V6[1])], HOT); pump(69.5, 10.0, 'r'); head(100.2, 22.0, 'u', HOT); head(99.2, V6[1], 'l', HOT, 0.6)
    pipe([(V6[0], V6[1]+s), (V6[0], 53.0), (93.0, 53.0)], HOT); head(V6[0], 48.0, 'u', HOT); head(95.0, 53.0, 'l', HOT)
    # cooled water leaving the tower coil into the side port of V5
    pipe([(67.0, 57.0), (V5[0], V5[1]+s)], COLD); head(V5[0], 46.0, 'd', COLD)
    # V5 left port: evaporative cooling, to the top port of V4
    pipe([(V5[0]-s, V5[1]), (V4[0], V5[1]), (V4[0], V4[1]+s)], COLD); head(45.0, V5[1], 'l', COLD)
    # V5 right port: chiller, to the condenser inlet
    pipe([(V5[0]+s, V5[1]), (76.0, V5[1]), (76.0, 20.0), (64.0, 20.0)], COLD); head(76.0, 30.0, 'd', COLD); head(67.5, 20.0, 'l', COLD)
    # evaporator outlet (side port of its frame) to the bottom port of V4
    pipe([(18.0, 20.0), (16.0, 20.0), (16.0, 24.5), (V4[0], 24.5), (V4[0], V4[1]-s)], COLD); head(V4[0], 31.0, 'u', COLD)
    # V4 side port to the bottom port of V3
    pipe([(V4[0]-s, V4[1]), (V3[0], V4[1]), (V3[0], V3[1]-s)], COLD); head(V3[0], 44.5, 'u', COLD)
    # free-cooler outlet (top of the header) to the top port of V3
    pipe([(32.4, 66.8), (32.4, 71.0), (26.5, 71.0), (26.5, 54.3), (V3[0], 54.3), (V3[0], V3[1]+s)], COLD); head(26.5, 62.0, 'd', COLD); head(24.2, 54.3, 'l', COLD, 0.7)
    # V3 side port: 18 C water back to the ITE heat exchanger
    pipe([(V3[0]-s, V3[1]), (19.0, V3[1]), (19.0, 59.0)], COLD); head(19.0, 56.2, 'u', COLD)
    valve3(*V1, 'v', 'l'); valve3(*V2, 'v', 'l'); valve3(*V3, 'v', 'l'); valve3(*V4, 'v', 'l')
    valve3(*V5, 'h', 'u'); valve3(*V6, 'h', 'u')

    fig.savefig(out, dpi=600, bbox_inches='tight', pad_inches=0.03)
    os.makedirs('figs/eps', exist_ok=True)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        fig.savefig('figs/eps/fig02_schematic.eps', format='eps', bbox_inches='tight', pad_inches=0.03)
    assert not [x for x in w if 'transparen' in str(x.message).lower()]
    fig.savefig('figs/eps/fig02_schematic.pdf', format='pdf', bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)


if __name__ == '__main__':
    fig2()
    print('figure 2 done')
