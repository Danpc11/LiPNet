"""The figures of the article are drawn in one style (src/figstyle.py, src/paper_figures.py)."""
import os, sys
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
import figstyle as fs
import paper_figures as pf


def test_one_meaning_per_colour():
    assert len({fs.PORTAL, fs.SINUSOID, fs.HEPATIC, fs.EVENT}) == 4
    assert fs.AQUA not in (fs.PORTAL, fs.SINUSOID, fs.HEPATIC) and fs.RED not in (fs.PORTAL, fs.SINUSOID, fs.HEPATIC), 'green is bile and red the artery in every atlas'
    for n in range(1, 7):
        r = fs.ramp(n); light = [sum(int(c[k:k + 2], 16) for k in (1, 3, 5)) for c in r]
        assert len(set(r)) == n and light == sorted(light, reverse=True), 'an ordered parameter runs from light to dark'


def test_every_figure_of_the_article_has_a_function():
    assert list(pf.FIGURES) == ['1', '2', '3', '4', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7']


def test_a_figure_is_drawn_at_the_width_of_the_page(tmp_path):
    import matplotlib.pyplot as plt
    pf.OUT = str(tmp_path); pf.figure_s2(); pf.figure_s4()                     # these two need no cached result
    for name, width in (('Figure_S2', fs.DOUBLE), ('Figure_S4', fs.SINGLE)):
        assert os.path.getsize(tmp_path / f'{name}.pdf') > 5000 and os.path.exists(tmp_path / f'{name}.png')
        w = plt.imread(tmp_path / f'{name}.png').shape[1]; assert abs(w - width * 300) <= 2
    assert plt.rcParams['pdf.fonttype'] == 42 and plt.rcParams['font.size'] == fs.FS
