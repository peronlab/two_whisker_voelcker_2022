%
% This is where ensemble membership is determined; if you want different methods of doing so, put them here
%
function ens_dat = get_two_whisker_ensemble_data(force_redo)
    if (nargin < 1) ; force_redo = 0 ; end

    params.exclude_matrix_name = 'excluded_for_correlations';
    params.force_redo = 0;
    dat = get_two_whisker_data(params);
    pt_dat = get_two_whisker_per_trial_data;
    ens_dat = identify_ensembles(dat,pt_dat, force_redo);

function ens_dat = identify_ensembles(dat, pt_dat, force_redo)
    global ens_dat;
    if (isempty(ens_dat) | force_redo)
        ens_dat_filename = [dat.settings.summary_data_root filesep 'two_whisker_ensemble_summary_data.mat'];
        if (exist(ens_dat_filename, 'file') & ~force_redo)
            disp (['Loading ' ens_dat_filename]);
            load(ens_dat_filename);
        else
            disp (['Generating ' ens_dat_filename]);

            % first grab corr matrix
            corr_dat_filename = [dat.settings.summary_data_root filesep 'two_whisker_corr_summary_data.mat'];
            global corr_dat;
            if (isempty(corr_dat))
                % if you get an error here, its because you must run figure_correlation_basic to aggregate correlation data
                load(corr_dat_filename);
            end

            % loop away
            corr_type = 'allTime';
%corr_type = 'whiskerNoTouchTrials';
            for ani=1:length(pt_dat)
                disp(['Processing ' dat.settings.animals(ani).name ' ... ']);
                excluded_ids = dat.anim_data(ani).types_by_id.excluded_for_correlations;
                for svi=1:length(pt_dat(ani).per_subvol)
                    disp(['    subvolume ' num2str(svi)]);

                    % pull correlation matrix and order by id ; fix subvolume
                    sv_ids = dat.anim_data(ani).ids(find(dat.anim_data(ani).subvol_idx == svi));
                    sv_ids = setdiff(sv_ids, excluded_ids);
                    fname = [dat.settings.data_path_root filesep dat.settings.animals(ani).name filesep 'session_neuropilone_final_merged' ...
                             filesep 'pairwiseCorrelations' filesep dat.settings.animals(ani).name sprintf('_vol_%02d_sess__%s.mat', svi + dat.settings.animals(ani).subvol_offset, corr_type)];
                    fdat = load(fname);
                    fdat = fdat.dat;

                    vali = find(ismember(fdat.cellIds, sv_ids));
                    M = fdat.corrMat(vali,vali);
%ens_dat.anim(ani).subvol(svi).all_corr.M_original_notouch = M;
%continue;

                    % settings for how to do ensemble identification
                    gen_ens_algorithm = 'ward_then_brute';
                    gen_ens_params = [];
                    gen_ens_params.min_ens_size = 3;
                    gen_ens_params.min_ens_overlap_frac_to_merge = 0.75;

                    % we set thresholds based on quantiles ...
                    [fM nM] = get_full_corrmat(M);
                    nM(find(isnan(nM))) = 0;
                    gen_ens_params.min_within_ens_overall_corr = 2*quantile(nM(:), .995);
                    gen_ens_params.min_individual_pairwise_corr = quantile(nM(:), .995)/2;
                    gen_ens_params.min_seed_individual_pairwise_corr = 2*quantile(nM(:), .995);
                    
%                    gen_ens_params.min_within_ens_overall_corr = 0.4;
%                    gen_ens_params.min_individual_pairwise_corr = 0.2;
%                    gen_ens_params.min_seed_individual_pairwise_corr = 0.4;
gen_ens_params                    
                    [ens_member_matrix ens_sort] = generate_ensembles(M, gen_ens_algorithm, gen_ens_params);
    
                    % the stuff we actually use:
                    ens_dat.anim(ani).subvol(svi).all_corr.M_original = M;
                    ens_dat.anim(ani).subvol(svi).all_corr.ids = fdat.cellIds(vali);
                    ens_dat.anim(ani).subvol(svi).all_corr.ens_matrix = ens_member_matrix;
                    ens_dat.anim(ani).subvol(svi).all_corr.ens_sort = ens_sort;
                end
            end

            save(ens_dat_filename,'ens_dat','-v7.3');
        end
    end


% --- THIS DOES NOTHING NOW ---
% a very early clustering attempt by yours truly
function [M_sorted grp_idx_original_sorting] = get_clusters_in_matrix(M, settings)
    % settings
    n_repeats = 1000; % how many times to run k-means
    km_seed_cells_per_group = 20;
    km_n_groups = max(2,floor(size(M,1)/km_seed_cells_per_group)); % 20 cells per group approx, at least 2
    min_frac_group = 0.5; % must be in a given group this proportion of runs to be assigned to said group

    if (nargin >= 2)
        fn = fieldnames(settings);
        for f=1:length(fn)
            eval([fn{f} '=settings.' fn{f} ';']);
        end
    end

    if (km_n_groups > size(M,1))
        M_sorted = [];
        grp_idx_original_sorting = [];
        return;
    end

    % restore to full mat
    for i1=1:size(M,1)    
        M(i1,i1) = 1;
        for i2=i1+1:size(M,2)    
            M(i2,i1) = M(i1,i2);
        end
    end

    % n_repeats of k-means 
    cluster_mat = nan*zeros(n_repeats, length(M));
    for n=1:n_repeats
        % do k-means, sorting 
        [irr irr nki] = sort_kmeans(M, km_n_groups);
        cluster_mat(n,:) = nki;
    end

    % now find the neurons that are consistently in a particular cluster
    grp_idx_original_sorting = nan*zeros(size(M,1),1);
    for g=1:km_n_groups
        count_mat = 0*cluster_mat;
        count_mat(find(cluster_mat == g)) = 1;
        frac_group_g = sum(count_mat)/length(cluster_mat);
        vali = find(frac_group_g > min_frac_group);

        unassigned_i = find(isnan(grp_idx_original_sorting));
        vali = intersect(vali, unassigned_i);
    
        grp_idx_original_sorting(vali) = g;
    end

    grp_idx_original_sorting(find(isnan(grp_idx_original_sorting))) = km_n_groups+1; % residual, effectively
    [irr n_sorti] = sort(grp_idx_original_sorting, 'ascend');

    M_sorted = M(n_sorti,n_sorti);

