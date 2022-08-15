This is the code for Voelcker et al., 2022, Nature Communications.

This employs 'processed' data; for raw data, please contact Peron lab.

1. Getting started: you will first need to modify get_two_whisker_settings.m to 
point to appropriate locations for all the files.  Ensure that the "src"
directory is in your MATLAB path.

2. Generating specific figures

Figure 1: figure_task_kinematics.m
Figure 2: figure_imaging_encoding.m
Figure 3: figure_rf_broadening.m 
Figure 5: figure_decode.m 
Figure 6: figure_correlations_basic.m 
Figure 7: figure_ensemble_basic.m 

3. Helper functions, may be of interest:

assign_types.m : this assigns touch type to neurons; see here if you want to
                 alter or understand how neurons are given specific touch type
get_two_whisker_data.m : core function for pulling up data.
get_two_whisker_settings.m : contains settings; you will have to alter this no
                             matter what.

4. Helper functions, unlikely to be of interest:

colormap_human.m : plotting utility
extern_plotSpread.m : plotting utility
get_best_subvol_list.m : provides list of best subvolumes for each animal for
                         analyses that require concurrent data (e.g.,
                         correlation)
get_ensemble_ids.m : returns ensemble ID for each ensemble; used in Fig. 7
get_full_corrmat.m : returns a full correlation matrix ;for Fig. 6
get_preferred_scores.m : returns touch score for *preferred* touch type
get_subgroups.m : to facilitate plotting (Fig. 3), returns cells in particular
                  groups.
get_two_whisker_decode_data.m : specifically gets decoding data (Fig. 5)
get_two_whisker_ensemble_data.m : pulls ensemble data (Fig. 7)
get_two_whisker_per_trial_data.m : pulls per-trial data (for Fig. 4)

get_two_whisker_touch_ca_ev_data.m : pulls event-aligned data (Fig. 4)
plot_3d_two_whiskers.m : generates a 3D plot (Figs. 2, 3)
plot_fancy_hist.m : plotter of histograms for Fig. 6
sort_ensembles.m : returns ensembles sorted in specific manner (Fig. 7)
struct_nan_blanks.m : utility function for generating vectors from
                      structure-embedded variables. 

