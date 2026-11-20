#!/usr/bin/env bash
# Rebuilds every result, every plot (one file each) and the calculator. Run from the repository root.
set -e
export PYTHONPATH=src
# --clean: delete every generated result, ignore every cache and recompute the cached tables of the model. The run to do before quoting a number (about 20 minutes).
if [ "$1" = "--clean" ]; then
  export LIPNET_REBUILD=1
  find results -maxdepth 1 -type f ! -name '.gitkeep' -delete
  echo "results/ emptied; every fit will be recomputed"
fi
python src/indices.py          # indices, stratified and hierarchical fits, validation
if [ -n "$LIPNET_REBUILD" ]; then python src/bayes.py; unset LIPNET_REBUILD; export LIPNET_REBUILT=1; fi   # primary hierarchical model (NUTS, several minutes)
python src/predictions.py      # the predictions of the model; refits bayes.py if the cached fits are not those of this table and code
python src/plots.py all        # 19 stand-alone plots -> results/<name>.pdf/.png (network_2d/3d recompute the model, ~1 min)
python src/build_app.py        # -> sfss_calculator.html
python src/ipd/run_ipd.py all --quiet     # ecological-bias correction -> results/ipd_*.json (~1 min)
python src/fem/run_fem.py all --quiet   # lobule finite elements -> results/fem_*.json and 4 figures (~1 min)
python src/model/liver.py > /dev/null    # the model applied to a human liver and a partial graft -> results/model_liver.json
python src/model/human.py > /dev/null    # the model in physical units, growth and regrowth -> results/model_human.json
python src/model/distension.py > /dev/null   # what a passively distensible bed can account for, as bounds -> results/model_distension.json
python src/audit.py                     # what the table and the fits support -> results/audit.json
if [ -n "$LIPNET_REBUILT" ] || [ ! -f results/campaign_check.json ]; then python src/model/verify_campaigns.py; fi   # the cached tables of the model against the code that writes them -> results/campaign_check.json
python src/paper_figures.py             # the eleven figures of the article in one style -> results/Figure_*.pdf/.png (docs/FIGURES.md)
python src/manifest.py                  # hashes of inputs and results, versions -> results/MANIFEST.json (python src/manifest.py check: numbers, figures, tables of the model)
# The tables of the model in data/ (graph_scaling_all, graph3d_all, certificate_small) are cached, not rebuilt above.
# verify_campaigns.py reruns the code behind them and compares every number (3 to 4 minutes); to write them again:
#   for R in 3.0 4.5 6.5 9.0; do python src/model/graph_scaling.py $R; done
#   for R in 3.0 4.0 5.0 6.5 8.5;  do python src/model/network3d.py $R results; done
#   python src/model/merge_campaigns.py
#   python src/model/certify.py && cp results/certificate_small.tsv data/
