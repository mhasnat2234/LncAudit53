#!/usr/bin/env bash
# Reproduces every table, figure and the PDF report from the raw inputs in data/raw/.
set -euo pipefail
cd "$(dirname "$0")/results/intermediate"
python3 ../../src/01_digitize_heatmap.py
python3 ../../src/02_core_analysis.py
python3 ../../src/03_sensitivity_and_simulation.py
python3 ../../src/04_make_figures.py
python3 ../../src/05_build_report.py
echo "Done. Report: report/GSE159185_Final_Report.pdf"
