"""Merge the per-size outputs of graph_scaling.py and network3d.py into the tables used by plots.py.
    python src/model/merge_campaigns.py   (reads results/graph_scaling_R*.tsv, results/graph3d_R*.tsv -> data/)"""
import glob, pandas as pd
for pat, out in (('results/graph_scaling_R*.tsv', 'data/graph_scaling_all.tsv'), ('results/graph3d_R*.tsv', 'data/graph3d_all.tsv')):
    f = sorted(glob.glob(pat))
    if f: pd.concat([pd.read_csv(x, sep='\t') for x in f]).sort_values(['b', 'N']).to_csv(out, sep='\t', index=False); print('wrote', out, 'from', len(f), 'files')
