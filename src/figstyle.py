"""One style for every figure of the article. The rules are written out in docs/FIGURES.md.

    import figstyle as fs; fs.use()

Size        180 mm (two columns) or 88 mm (one column) wide; PDF in vector with embedded TrueType, PNG at 300 dpi
Type        one sans family, 7 pt for axes and ticks, 6 to 6.5 pt inside the panels, panel letters 8 pt bold capitals
Lines       axes 0.5 pt, data 1.0 to 1.2 pt, reference lines 0.5 to 0.7 pt dashed; no top or right spine, no title, no legend frame
Text        always ink or muted ink, never the colour of a series; variables in italics, exponents as superscripts
Colour      one meaning per colour in the whole article (PORTAL, SINUSOID, HEPATIC, EVENT below);
            an ordered parameter (b, h, baseline gradient) takes one hue from light to dark (ramp)
"""
import logging, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

MM = 1 / 25.4
DOUBLE, SINGLE = 180 * MM, 88 * MM                       # figure widths

INK, MUTED, GRID, BAND, NEUTRAL = '#0b0b0b', '#52514e', '#e6e5e0', '#f0efec', '#d9d8d2'
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, RED, PURPLE = '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#e34948', '#b33fb0'
CATEGORICAL = (BLUE, ORANGE, AQUA, YELLOW, MAGENTA)      # identity of unordered series, always in this order
PORTAL, SINUSOID, HEPATIC = PURPLE, ORANGE, BLUE         # as in the atlases: portal vein purple, hepatic veins blue; green is left to bile, red to the artery
EVENT = RED                                              # the one recipient or group singled out for an adverse outcome
MODEL_TINT = dict(color=BLUE, alpha=0.16, lw=0)          # the range the network model reaches
_RAMP = ('#b7d3f6', '#86b6ef', '#5598e7', '#2a78d6', '#184f95', '#0d366b')

FS, FS_IN, FS_SMALL, FS_LETTER = 7, 6.5, 6, 8            # axes and ticks; text inside a panel; second lines; panel letters
LW, LW_REF, LW_AXIS = 1.2, 0.7, 0.5
DASH = (0, (4, 2))
DOT = dict(mec='white', mew=0.6)                         # a surface ring around every marker


def ramp(n):
    """n steps of the sequential hue, light to dark, for an ordered parameter."""
    pick = {1: (3,), 2: (1, 4), 3: (1, 3, 5), 4: (1, 2, 3, 5), 5: (0, 1, 2, 3, 5), 6: (0, 1, 2, 3, 4, 5)}[n]
    return [_RAMP[k] for k in pick]


def family():
    have = {f.name for f in font_manager.fontManager.ttflist}
    return next((f for f in ('Liberation Sans', 'Arial', 'Helvetica', 'DejaVu Sans') if f in have), 'sans-serif')


def use():
    f = family()
    plt.rcParams.update({
        'font.family': [f, 'DejaVu Sans'], 'font.size': FS, 'axes.labelsize': FS, 'xtick.labelsize': FS, 'ytick.labelsize': FS, 'legend.fontsize': FS_IN,
        'mathtext.fontset': 'custom', 'mathtext.rm': f, 'mathtext.it': f + ':italic', 'mathtext.bf': f + ':bold', 'mathtext.default': 'it',
        'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK,
        'axes.linewidth': LW_AXIS, 'xtick.major.width': LW_AXIS, 'ytick.major.width': LW_AXIS, 'xtick.minor.width': 0.35, 'ytick.minor.width': 0.35,
        'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.minor.size': 1.5, 'ytick.minor.size': 1.5, 'axes.labelpad': 2.5,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.unicode_minus': True,
        'lines.linewidth': LW, 'lines.solid_capstyle': 'round', 'lines.markersize': 4.5,
        'legend.frameon': False, 'legend.handlelength': 1.6, 'legend.labelspacing': 0.3, 'legend.borderaxespad': 0.2, 'legend.handletextpad': 0.5,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'savefig.dpi': 300, 'figure.dpi': 150})


def figure(height_mm, width=None):
    return plt.figure(figsize=(width or DOUBLE, height_mm * MM))


def axes(fig, left, bottom, width, height, **kw):
    """Axes placed in millimetres from the lower left corner of the figure."""
    W, H = fig.get_size_inches() / MM
    return fig.add_axes((left / W, bottom / H, width / W, height / H), **kw)


def letter(fig, ax, s, dx=9.0, dy=3.5):
    """Panel letter dx mm to the left of the axes and dy mm above them."""
    W, H = fig.get_size_inches() / MM; b = ax.get_position()
    fig.text(b.x0 - dx / W, b.y1 + dy / H, s, fontsize=FS_LETTER, fontweight='bold', va='center', ha='left', color=INK)


def bare(ax, keep='bottom'):
    """An axis that carries rows of labels instead of a scale on its left side."""
    ax.set_yticks([]); ax.spines['left'].set_visible(False)


def save(fig, name, out):
    os.makedirs(out, exist_ok=True)
    fig.savefig(f'{out}/{name}.pdf'); fig.savefig(f'{out}/{name}.png', dpi=300); plt.close(fig); print('saved', name)
