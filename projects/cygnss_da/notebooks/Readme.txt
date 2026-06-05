README file

Created by Andrew Fox, Rolf Reichle and Qing Liu, Feb 2026

===============================================================================

This directory contains data used in the following paper:

Fox, A. M., R. H. Reichle and Q. Liu, 2026: Assimilation of CYGNSS Soil Moisture
Observations into the NASA GEOS Land Data Assimilation System,
Journal of Hydrometeorology, submitted.

As a condition of using the data provided here, please reference this paper.

===============================================================================

Top-level directory structure

  figure_data/     Derived datasets and analysis outputs used to generate the
                   figures in the paper
  model_config/   Experiment-specific configuration information for the
                   GEOS-LDAS runs
  model_output/   Subsetted GEOS-LDAS model output and observation-space
                   assimilation diagnostics

===============================================================================

Experiments and data coverage

Data are provided for the following experiments, which are described in detail
in Fox et al. (2026):

  CNTL     --> Model-only control simulation
  CYG_DA   --> CYGNSS soil moisture data assimilation
  SSA_DA   --> Multi-sensor DA excluding CYGNSS
               (SMAP Tb + SMOS Tb + ASCAT SM)
  ALL_DA   --> Multi-sensor data assimilation
               (CYGNSS SM + SMAP Tb + SMOS Tb + ASCAT SM)

Data are provided for the period from August 2018 to June 2024.

===============================================================================
===============================================================================

model_output/

This directory contains subsetted ensemble-average GEOS Land Data Assimilation System output
and ensemble-average observation-space diagnostics.

Directory organization:

  model_output/<Expt>/Y<year>/M<month>/D<day>/

Two output file types are provided:

  (1) <Expt>.ens_avg.ldas_ObsFcstAna.<timestamp>.bin
  (2) <Expt>.tavg24_1d_lnd_Nt.<timestamp>.nc4

-------------------------------------------------------------------------------

(1) <Expt>.ens_avg.ldas_ObsFcstAna.<timestamp>.bin

These binary files contain ensemble-average observation-space diagnostics from
the GEOS Land Data Assimilation System. Eight files, corresponding to the
3-hourly assimilation windows, are available each day.

The file format and contents are identical across all experiments and
consistent with previous GEOS-LDAS releases.

A a set of Python scripts are provided in https://github.com/GEOS-ESM/GEOSldas_GridComp/tree/develop/GEOSldas_App/util/postproc/ObsFcstAna_stats for reading these files.
The variables stored in each file are:

  (i)    Total number of observations in file
  (ii)   Time stamp entry
  (iii)  Observation assimilation flag
  (iv)   Observation species
  (v)    Tile number
  (vi)   Longitude
  (vii)  Latitude
  (viii) Observation value
  (ix)   Observation variance
  (x)    Observation-space model forecast value
  (xi)   Observation-space model forecast variance
  (xii)  Observation-space model analysis value
  (xiii) Observation-space model analysis variance

-------------------------------------------------------------------------------

(2) <Expt>.tavg24_1d_lnd_Nt.<timestamp>.nc4

These files are provided in self-describing NetCDF-4 (nc4) format and contain
daily time-averaged land surface diagnostics on the native GEOS land tile grid.

Dimensions:
  tile = 70773
  time = 1

Coordinate variables:
  lon(tile)    Longitude of land tile centers [degrees east]
  lat(tile)    Latitude of land tile centers  [degrees north]

Science variables (ensemble mean):

  SFMC(time,tile)
    Surface soil moisture content
    Units: m3 m-3

  RZMC(time,tile)
    Root-zone soil moisture content
    Units: m3 m-3

  TSURFLAND(time,tile)
    Surface temperature over land, including snow
    Units: K

  TSOIL1(time,tile)
    Soil temperature of the uppermost soil layer
    Units: K

  PRECTOTCORRLAND(time,tile)
    Total precipitation rate over land
    Units: kg m-2 s-1

All variables represent 24-hour time-averaged values ending at the file
timestamp.

Missing values are indicated by _FillValue = 1.e+15.

===============================================================================
===============================================================================

figure_data/

This directory contains derived datasets and analysis outputs used to generate
the figures and tables in Fox et al. (2026).

The subdirectory

  cygnss_da_combined_inputs/

contains processed model output, evaluation data, and
derived statistics used in the evaluation of the experiments. These datasets
are derived from the native model_output files and include, as applicable:

  • Spatially aggregated statistics (e.g., global and regional maps)
  • Temporally aggregated statistics (e.g., daily, monthly, or climatological)
  • Triple-collocation and instrumental-variable diagnostics
  • Validation metrics against Cal/Val sites, SCAN, and USCRN
  • Ancillary fields used for stratification (e.g., aridity indices)

Files are provided in a mixture of NetCDF (.nc4), MATLAB (.mat), Python pickle
(.pkl), and binary formats, reflecting the analysis workflows used to generate
the paper figures. A manifest file is included to document the contents of this
directory.

The Jupyter notebooks

  cygnss_da_combined_paper_figures.ipynb

  cygnss_da_combined_revised_paper_figures.ipynb

reproduce the figures in the manuscript using the datasets in
cygnss_da_combined_inputs/. The first notebook corresponds to the originally
submitted paper; the second produces figures for the revised manuscript.
Figure numbering and content correspond directly to the respective versions
of the paper.

===============================================================================
===============================================================================

model_config/

This directory contains experiment-specific configuration information for the
GEOS Land Data Assimilation System runs used in Fox et al. (2026).

Subdirectories are provided for each experiment:

  CNTL
  CYG_DA
  SSA_DA
  ALL_DA

Each experiment directory contains the configuration files that define the
assimilation setup, including (as applicable):

  • Experiment metadata and descriptions
  • Namelist files controlling ensemble propagation and update settings
  • User-level configuration script to use with the ObsFcstAna Python processing scripts

These files document the  configuration used for each experiment and are
provided to support transparency and reproducibility, rather than to serve as a
turnkey re-run package.

===============================================================================

The total data volume is approximately 85 GB.

=============================== EOF ==========================================