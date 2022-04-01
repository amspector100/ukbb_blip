# An application of BLiP to UK Biobank Data

This repository contains the code used in [Spector and Janson (2022)](https://arxiv.org/abs/2203.17208) to perform a genetic fine-mapping analysis on UK Biobank data. In particular, we applied BLiP on top of SuSiE models fit by [Weissbrod et al. (2019)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7710571/).

To run this analysis, you will need to install [pyblip](https://github.com/amspector100/pyblip).

# Replicating the main analysis

To rerun the analysis, do as follows:

1. Download the requesite data by running ``python3.9 getdata.py``. This may take a while as some files are quite large.
2. Run BLiP on top of SuSiE by running ``python3.9 blip_wrapper.py``. This should take roughly ten minutes. However, the vast majority of the time is spent loading the data into memory or finding the rejection sets from SuSiE. Running BLiP on top of SuSiE should take less than one minute per trait.
3. To process the outputs, run ``python3.9 postprocessing.py``. This will be slow due to the need to merge large files.
4. The plots can then be generated using the ``ukbb_blip_final_plots.ipynb`` notebook.

Final plots will be located in ``output/final/plots/`` and the final set of rejections will be located in ``output/final/rejections``.

The data in final/rejections has already been pushed to GitHub, so it should be possible to run ``ukbb_blip_final_plots.ipynb`` without running any of the prior analysis.

# Rerunning the simulations

Simulations can be run using the ``blip_sims.py`` script. The specific parameters of our simulations are contained in ``sims.sh``, which can be run using the command ``bash sims.sh``.


