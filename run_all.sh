#!/usr/bin/env bash
# Rebuilds every result, every plot (one file each) and the calculator. Run from the repository root.
set -e
export PYTHONPATH=src
python src/indices.py          # indices, stratified and hierarchical fits, validation
python src/bayes.py            # primary hierarchical model, priors, overlap sets, validation (NUTS, a few minutes)
python src/predictions.py      # the four predictions of the model and the reference sweep -> results/predictions.json
python src/plots.py all        # 19 stand-alone plots -> results/<name>.pdf/.png (network_2d/3d recompute the model, ~1 min)
python src/build_app.py        # -> sfss_calculator.html
python src/ipd/run_ipd.py all --quiet     # ecological-bias correction -> results/ipd_*.json (~1 min)
python src/fem/run_fem.py all --quiet   # lobule finite elements -> results/fem_*.json and 4 figures (~1 min)
# Long model campaigns behind the cached tables in data/ (optional, 1-2 h):
#   for R in 3.0 4.5 6.5 9.0 12.0; do python src/model/graph_scaling.py $R; done
#   for R in 3.0 4.0 5.0 6.5 8.5;  do python src/model/network3d.py $R results; done
#   python src/model/merge_campaigns.py
#   python src/model/certify.py && cp results/certificate_small.tsv data/
