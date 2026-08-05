### Datasets

This directory contains the datasets accompanying the following study:

Ji, Q., Dey, I., & Dunham, E. (2026). Turbulent seismoacoustic imprints during a hurricane landfall. Science. https://doi.org/10.1126/science.adt7323

Dataset published at: https://doi.org/10.25740/zk709vr3334

The studied event is Hurricane Isaac (2012). Place this `Data/` folder at the top level of the code archive (next to `Notebooks/` and `Codes/`); the Jupyter Notebooks locate it automatically (see Notebooks/README.txt). Total size is about 4.5 GB.

- Isaac_HQ.png: Satellite basemap image for the maps (Maps.ipynb).

### ASOS (Automated Surface Observing Systems)

Surface wind records from ASOS weather stations, each paired with a nearby TA seismic station.

- {TA station}_{ASOS station}.csv (e.g. 645A_HUM.csv): ASOS data closest to the given TA station.
- 645A_HUM_Jan.csv: ASOS data for January 2012 at station HUM (closest to TA.645A).
- LA_meta.csv: Metadata of the Louisiana ASOS network.

### TA (Transportable Array)

Seismic and infrasound data and derived spectral products (~2.2 GB).

- TA_list.csv: Station list (network, station name, longitude, latitude).

- mseeds/: Instrument-corrected time series (downloaded from EarthScope/IRIS in the notebooks).
  - TA_645A.mseed, TA_544A.mseed: Pressure (LDF) and seismic displacement (LHZ/LHN/LHE) channels at 1 Hz, 2012-08-24 to 2012-09-03.
  - BDF_645A.mseed: Infrasound (BDF) channel at 40 Hz for station 645A.
  - BDF data of the other stations are not included here; they can be re-downloaded with calc_turb.ipynb.

- spectrogram/: Wavelet (fCWT) spectral products from calc_spectrogram.ipynb - PSD, cross-spectra and coherence between all channel pairs (P, Z, N, E).
  - TA_645A.npz, TA_544A.npz: During the hurricane passage.
  - TA_645A_Jan.npz: Reference quiet period (January 2012).

- turbulence/: Turbulent pressure products from calc_turb.ipynb for 19 TA stations.
  - Cpp_{station}.npz: Pressure PSD and RMS amplitude in 15-minute windows.
  - Dpp_{station}.npz: Second-order pressure structure functions in 15-minute windows.

### fcmp_tower (Florida Coastal Monitoring Program portable towers)

- raw_data/: Empty placeholder for the raw 32-Hz wind tower records, which are available upon request (see Acknowledgments of the paper).
- turbulence/: Wind spectra (Suu_{height}.npz) and structure functions (Duu_{height}.npz) at heights of 5, 7.5, 10, 12.5 and 15 m, computed from the raw records by calc_obs_spec.ipynb.

### Isaac (hurricane observations)

Atmospheric observations during the landfall of Hurricane Isaac (~620 MB, mostly ERA5).

- obs_spectra.npz: Observed spectra during the LES modeling hour - wind tower velocity/temperature spectra at several heights, and infrasound/seismic PSD with P-Z coherence at TA.645A (from calc_obs_spec.ipynb).
- dropsonde/: NOAA aircraft dropsonde profiles for missions 20120828H2 and 20120829A1 (raw .avp and quality-controlled *QC.nc files).
- ERA5/: ERA5 reanalysis around the landfall - single-level files (data_stream-oper*.nc, data_stream-wave.nc) and pressure-level file (data_pres_level.nc).
- NEXRAD/: Wind profiles from NEXRAD radar station KLIX.
  - KLIX_VAD.nc: Gathered velocity-azimuth display (VAD) wind profiles.
  - KLIX_NVW/: Raw NVW wind profile files.
- radar/: Radar wind analyses for missions 20120828H1 and 20120828H2 in GrADS format (.ctl/.dat, read via xgrads in CM1_inputs.ipynb).
- radiosonde/: Radiosonde soundings of station 72233 (Slidell, LA) and the IGRA2 station list.

### TC_track (hurricane track database)

- IBTrACS.ALL.v04r00.nc: IBTrACS best-track database (v04r00), used by TC_analysis/tc_process.py to obtain the hurricane track. A newer version (v04r01) is also provided, which will slightly change the track in the zoom-in plot Figure 1B.

### CM1 (large-eddy simulations)

CM1 hurricane boundary layer simulations and elastic modeling results (~1.9 GB). Simulation directories are named V{wind}_Cd{drag x 1000} (e.g. V42_Cd22: input gradient wind 42 m/s, drag coefficient 0.022) or R{radius} (radius to hurricane center in km).

- inputs/: CM1 configuration - input_sounding, namelist.input (main run), namelist.input_dense (1 Hz output run), and cm1_config.txt (CM1 version r21.0 and configuration summary). Same as in Codes/CM1_run.

- outputs/: Reference simulation (V42_Cd22) outputs and elastic modeling results.
  - cm1out_stats.nc, base_state.nc: Domain statistics and base state.
  - evolution.mat, diag_profile_6.mat, analysis_profile.mat: Temporal evolution and time/domain-averaged vertical profiles of diagnostic variables (from Codes/CM1_analysis/matlab).
  - cm1out_prs_15m.nc: Pressure field extracted at 15 m height from the dense output run (1 Hz); input to the elastic modeling.
  - LES_ref.nc, LES_snapshot.npz: Modeled surface pressure and vertical displacement fields (full time series and a single snapshot; Modeling.ipynb).
  - spectra_ref.npz: Modeled spectra of the reference case (see spectra/ below).

- spectra/: Spectra of LES fields and elastic modeling results (Modeling.ipynb, Transfer.ipynb).
  - spectra_{prs,th,wind}_all.npz: LES pressure, temperature and wind spectra at height levels (gathered by Codes/CM1_analysis/python/gather_spec.py).
  - spectra_{prs,th,wind}_all_dx10m.npz: Same for the higher-resolution run with 10-m grid spacing.
  - spectra_ref.npz: Modeled pressure/displacement PSD and P-Z transfer function for the fitted velocity model.
  - spectra_site.npz: Same for the initial velocity model from the site survey.
  - spectra_hs_645A.npz, spectra_hs_544A.npz: Same for halfspace (Vs30) models at each station.
  - spectra_kernel_{1-7}.npz: Results for layer-perturbed velocity models (layer sensitivity; from Codes/elastic_modeling/model_kernel.m).

- corr_xt/: Space-time correlation functions of the surface pressure field along different directions (corr_{angle}deg.mat, from Codes/CM1_analysis/matlab/convec_vel.m), used to estimate the convection velocity.

- param_search/: Grid search over input wind V (38-46 m/s) and drag coefficient Cd (0.018-0.026).
  - V{V}_Cd{Cd}/: Outputs of each run (same file types as outputs/ above, without the dense-run products).
  - param_search.mat: Database gathered from all runs (from Codes/CM1_analysis/matlab/create_db.m).
  - psd_30-200s.csv: Mean pressure PSD in the 30-200 s band for each run, used for the misfit (Modeling.ipynb).

- param_search_radius/: Simulations at different radii to the hurricane center (R80-R300 km).
  - R{radius}/: Outputs of each run.
  - param_search_radius.mat: Database gathered from all runs.
  - convec_summary.txt: Convection velocity for each radius.

- prs_movie.mp4: Movie of the simulated surface pressure fluctuations.

### vel_model (elastic velocity models)

- vel_model_fit.csv: Velocity model for the main text (used by Codes/elastic_modeling and Modeling.ipynb).
- vel_model_site.csv: Initial velocity model from the site survey.
- vel_model_site.png: Figure of the site velocity model.
