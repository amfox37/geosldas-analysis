# matlab2python

Staging area for migrating legacy MATLAB workflows into Python. The `python_calc_plot_ObsFcstAna` subfolder already contains Python equivalents of common MATLAB notebooks, while `matlab_postprocess/` and `obs_scaling_params/` retain reference MATLAB assets used during translation.

## Key assets
- `python_calc_plot_ObsFcstAna/OFA_figure_maker_041125.ipynb` – Python port of the OFA figure generator.
- `python_calc_plot_ObsFcstAna/ObsFcstAna_notebook.ipynb` – general ObsFcstAna reader/plotter translated from MATLAB.
- `scripts/OMF_monthly_040125.py` – example of moving monthly OmF calculations from MATLAB into Python.

Use this directory when cleaning up MATLAB code or experimenting with new Python replacements before integrating them into `utils/` or project-specific folders.
