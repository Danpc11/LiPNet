#!/usr/bin/env bash
# Rebuilds every result, every plot (one file each) and the calculator. Run from the repository root.
set -e
export PYTHONPATH=src
python src/indices.py          # indices + pooled logistic fit -> results/
python src/plots.py all        # 14 individual plots -> results/<name>.pdf/.png (network_2d/3d recompute the model, ~1 min)
python src/build_app.py        # -> sfss_calculator.html
# Long campaigns behind the cached tables in data/ (optional, 1-2 h):
#   for R in 3.0 4.5 6.5 9.0 12.0; do python src/model/graph_scaling.py $R; done      # -> results/graph_scaling_R*.tsv
#   for R in 3.0 4.0 5.0 6.5 8.5;  do python src/model/network3d.py $R results; done  # -> results/graph3d_R*.tsv
#   python src/model/certify.py                                                          # -> results/certificate_small.tsv
