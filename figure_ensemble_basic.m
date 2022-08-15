%
% a) corrmat sorted by ensemble ; include ens membership, score, etc ; float nonmembers
% b) dff from these neurons, as well ; float nonmembers
% c) example plane
% d) ensemble size: fraction of cells in ensemble
% e) ensemble overlap
% f) frac mw/sw cells in ensembles, by lyr
%% f) what frac of touch evoked ca events are from ensemble members?
%% g) ensemble decoding compared to others
%
function figure_ensemble_basic
    individual_plot('an280201','L2', 1);

function individual_plot (corr_mat_ex_animal, corr_mat_ex_layer, plot_raw)
    % --- data gather 
    params.exclude_matrix_name = 'excluded_for_correlations';
    params.force_redo = 0;

    global dat;
    global ens_dat;
    global pt_dat;
    if (isempty(ens_dat))
        dat = get_two_whisker_data(params);
        pt_dat = get_two_whisker_per_trial_data;
        ens_dat = get_two_whisker_ensemble_data;
    end
    [subvol_idx l_tag] = get_best_subvol_list (dat);

    % --- setup figure
    fsize = 10;
    max_xy = 700;
    fh_main = figure ('Position',[0 0 1800 900], 'Name', [corr_mat_ex_animal '-' corr_mat_ex_layer]);

    top_y = 0.55;

    example_ensemble_corr_mat_ax = axes(fh_main,'Position', [.025 top_y .15 .3]);
    example_ensemble_corr_mat_cbar_ax = axes(fh_main,'Position', [.18 top_y .025/3 .05]);

    example_ensemble_membership_ax = axes(fh_main,'Position', [.19 top_y .025 .3]);
    example_pref_whisker_ax = axes(fh_main,'Position', [.95 .1 2*.025/3 .3]);
    example_touch_type_ax = axes(fh_main,'Position', [.22 top_y .025 .3]);
    example_ensemble_overall_score_ax = axes(fh_main,'Position', [.265 top_y .025/3 .3]);
    example_ensemble_w1_score_ax = axes(fh_main,'Position', [.28 top_y .025/3 .3]);
    example_ensemble_w1p_score_ax = axes(fh_main,'Position', [.29 top_y .025/3 .3]);
    example_ensemble_w1r_score_ax = axes(fh_main,'Position', [.30 top_y .025/3 .3]);
    example_ensemble_w2_score_ax = axes(fh_main,'Position', [.315 top_y .025/3 .3]);
    example_ensemble_w2p_score_ax = axes(fh_main,'Position', [.325 top_y .025/3 .3]);
    example_ensemble_w2r_score_ax = axes(fh_main,'Position', [.335 top_y .025/3 .3]);
    
    example_ensemble_overall_score_cbar_ax = axes(fh_main,'Position', [.265 top_y+.32 .025/3 .05]);
    example_ensemble_w1_score_cbar_ax = axes(fh_main,'Position', [.28 top_y+.32 .025/3 .05]);
    example_ensemble_w1p_score_cbar_ax = axes(fh_main,'Position', [.29 top_y+.32 .025/3 .05]);
    example_ensemble_w1r_score_cbar_ax = axes(fh_main,'Position', [.30 top_y+.32 .025/3 .05]);
    example_ensemble_w2_score_cbar_ax = axes(fh_main,'Position', [.315 top_y+.32 .025/3 .05]);

    example_ensemble_dff_ax = axes(fh_main,'Position', [.35 top_y .3 .3]);
    example_ensemble_dkappa_w2_ax = axes(fh_main,'Position', [.35 .86 .3 .03]);
    example_ensemble_dkappa_w1_ax = axes(fh_main,'Position', [.35 .895 .3 .03]);
    example_ensemble_theta_w1_ax = axes(fh_main,'Position', [.35 .93 .3 .06]);
   
    summ_frac_in_ens_ax = axes(fh_main,'Position', [.025 .4 .1 .1]);
    summ_ens_size_ax = axes(fh_main,'Position', [.025 .25 .1 .1]);
    summ_ens_count_ax = axes(fh_main,'Position', [.025 .1 .1 .1]);

    summ_frac_in_ens_by_type_ax(1) = axes(fh_main,'Position', [.15 .4 .1 .1]);
    summ_frac_in_ens_by_type_ax(2) = axes(fh_main,'Position', [.15 .25 .1 .1]);
    summ_frac_in_ens_by_type_ax(3) = axes(fh_main,'Position', [.15 .1 .1 .1]);

    summ_hist_frac_in_t_ens_by_type_ax(1) = axes(fh_main,'Position', [.275 .4 .1 .1]);
    summ_hist_frac_in_t_ens_by_type_ax(2) = axes(fh_main,'Position', [.275 .25 .1 .1]);
    summ_hist_frac_in_t_ens_by_type_ax(3) = axes(fh_main,'Position', [.275 .1 .1 .1]);

    summ_frac_in_t_ens_ax = axes(fh_main,'Position', [.4 .4 .1 .1]);
    summ_t_ens_size_ax = axes(fh_main,'Position', [.4 .25 .1 .1]);
    summ_t_ens_count_ax = axes(fh_main,'Position', [.4 .1 .1 .1]);

    summ_event_frac_in_t_ens_ax = axes(fh_main,'Position', [.525 .4 .1 .1]);
    summ_prob_resp_on_touch_t_ens_ax = axes(fh_main,'Position', [.525 .25 .1 .1]);
    summ_prob_resp_on_touch_non_t_ens_ax = axes(fh_main,'Position', [.525 .1 .1 .1]);

    summ_overlap_non_t_ens_ax = axes(fh_main,'Position', [.65 .4 .1 .1]);
    summ_prob_resp_at_least_one_on_touch_t_ens_ax = axes(fh_main,'Position', [.65 .25 .1 .1]);
    summ_prob_resp_at_least_one_on_touch_non_t_ens_ax = axes(fh_main,'Position', [.65 .1 .1 .1]);

    summ_overlap_t_ens_ax = axes(fh_main,'Position', [.775 .4 .1 .1]);
    summ_pairwise_dist_in_t_ens_ax = axes(fh_main,'Position', [.775 .25 .1 .1]);
    summ_pairwise_dist_out_t_ens_ax = axes(fh_main,'Position', [.775 .1 .1 .1]);

    example_ensemble_spatial_ax = axes(fh_main, 'Position', [.8 top_y .15 .3]);

    % --- visualize correlation matrices
    ani = find(strcmp(corr_mat_ex_animal, {dat.settings.animals.name}));
    li = find(strcmp(corr_mat_ex_layer, l_tag));
    svi = subvol_idx(ani, li);
    disp(['Plotting ' corr_mat_ex_animal ' ' corr_mat_ex_layer ', svi: ' num2str(svi)]);

    % apply restriction to ens_dat - laminar but this also enforces no duplicate cells/ROIs
    l_ids = {dat.anim_data(ani).types_by_id.L2, dat.anim_data(ani).types_by_id.L3, dat.anim_data(ani).types_by_id.L4};
    vali = find(ismember(ens_dat.anim(ani).subvol(svi).all_corr.ids, l_ids{li}));
    if (length(vali) ~= length(ens_dat.anim(ani).subvol(svi).all_corr.ids)) % trimmage needed?
        ens_dat.anim(ani).subvol(svi).all_corr.ids = ens_dat.anim(ani).subvol(svi).all_corr.ids(vali);
        ens_dat.anim(ani).subvol(svi).all_corr.M_original = ens_dat.anim(ani).subvol(svi).all_corr.M_original(vali,vali);
        ens_dat.anim(ani).subvol(svi).all_corr.ens_matrix = ens_dat.anim(ani).subvol(svi).all_corr.ens_matrix(:,vali);

        [eT sort_eT] = sort_ensembles(ens_dat.anim(ani).subvol(svi).all_corr.ens_matrix, ens_dat.anim(ani).subvol(svi).all_corr.M_original, 'sumcorr-first');
        ens_dat.anim(ani).subvol(svi).all_corr.ens_sort = sort_eT;
    end

    % now pull everything
    M = get_full_corrmat(ens_dat.anim(ani).subvol(svi).all_corr.M_original);
    ens_member_mat = ens_dat.anim(ani).subvol(svi).all_corr.ens_matrix;
    ens_sort = ens_dat.anim(ani).subvol(svi).all_corr.ens_sort;
    memberi = find(nansum(ens_member_mat)>0);
    sort_memberi = find(ismember(ens_sort,memberi));
    ens_sort = ens_sort(sort_memberi);

    M_sorted = M(ens_sort,ens_sort);

    axes(example_ensemble_corr_mat_ax);
    imagesc([M_sorted], [0 0.75]);
    cm = colormap_human(256);
    cm (end,:) = [0 0 0];

    colormap(example_ensemble_corr_mat_ax, cm);

    set(example_ensemble_corr_mat_ax, 'TickDir', 'out', 'FontSize', fsize, 'XTick',[]);
    ylabel('Cell');

    axes(example_ensemble_corr_mat_cbar_ax);
    imagesc((0:.01:1)', [0 1]);
    colormap(example_ensemble_corr_mat_cbar_ax, cm);
    set(example_ensemble_corr_mat_cbar_ax, 'TickDir', 'out', 'YTick', [], 'XTick',[]);



    % ensemble membership "matrix"
    axes(example_ensemble_membership_ax);
    imagesc(ens_member_mat(:,ens_sort)', [0 1]);
    colormap(example_ensemble_membership_ax, flipud(gray));
    set(example_ensemble_membership_ax, 'TickDir', 'out', 'FontSize', fsize, 'YTick',[]);
    xlabel(example_ensemble_membership_ax,'ensemble');
 
    % whisker typeZ
    w1i = find(ismember(ens_dat.anim(ani).subvol(svi).all_corr.ids, dat.anim_data(ani).types_by_id.w1ExclusiveTouchTrials));
    w2i = find(ismember(ens_dat.anim(ani).subvol(svi).all_corr.ids, dat.anim_data(ani).types_by_id.w2ExclusiveTouchTrials));
    mwi = find(ismember(ens_dat.anim(ani).subvol(svi).all_corr.ids, dat.anim_data(ani).types_by_id.mwCellIds));
    w1i = setdiff(w1i,mwi);
    w2i = setdiff(w2i,mwi);
    unidirwi = find(ismember(ens_dat.anim(ani).subvol(svi).all_corr.ids, intersect(dat.anim_data(ani).types_by_id.unidirCellIds, union(dat.anim_data(ani).types_by_id.w1ExclusiveTouchTrials, dat.anim_data(ani).types_by_id.w2ExclusiveTouchTrials))));
    bidirwi = find(ismember(ens_dat.anim(ani).subvol(svi).all_corr.ids, dat.anim_data(ani).types_by_id.bidirCellIds));

    type_mat = ones(size(M,1),5,3);

    % column one - whisker etc
    dat.settings.colors
    for i=1:3
        type_mat(bidirwi,2,i) = dat.settings.colors.bidiColor(i);
        type_mat(unidirwi,1,i) = dat.settings.colors.uniColor(i);
        type_mat(mwi,3,i) = dat.settings.colors.mwColor(i);

        type_mat(w1i,4,i) = dat.settings.colors.w1Color(i);
        type_mat(mwi,4,i) = dat.settings.colors.w1Color(i);

        type_mat(w2i,5,i) = dat.settings.colors.w2Color(i);
        type_mat(mwi,5,i) = dat.settings.colors.w2Color(i);

    end

    axes(example_pref_whisker_ax);
    imagesc(type_mat(ens_sort,4:5,:),[0 1]);
    set(example_pref_whisker_ax, 'TickDir', 'out', 'FontSize', fsize, 'YTick',[],'XTick',[]);
    xlabel(example_pref_whisker_ax,'PW');

    axes(example_touch_type_ax);
    imagesc(type_mat(ens_sort,1:3,:),[0 1]);
    set(example_touch_type_ax, 'TickDir', 'out', 'FontSize', fsize, 'YTick',[],'XTick',[]);
    xlabel(example_touch_type_ax,'Type');

    % whisker encoding scores
    celli = find(ismember(dat.anim_data(ani).ids, ens_dat.anim(ani).subvol(svi).all_corr.ids));
    overall_score = [dat.anim_data(ani).touch_scores(celli).allTrials];
    w1_score = [dat.anim_data(ani).touch_scores(celli).w1ExclusiveTouchTrials];
    w1p_score = [dat.anim_data(ani).touch_scores(celli).w1ProExclusiveTouchTrialsBalanced]; % balance b/c rare 
    w1r_score = [dat.anim_data(ani).touch_scores(celli).w1RetExclusiveTouchTrialsBalanced];
    w2_score = [dat.anim_data(ani).touch_scores(celli).w2ExclusiveTouchTrials];
    w2p_score = [dat.anim_data(ani).touch_scores(celli).w2ProExclusiveTouchTrialsBalanced];
    w2r_score = [dat.anim_data(ani).touch_scores(celli).w2RetExclusiveTouchTrialsBalanced];

    axes(example_ensemble_overall_score_ax);
    imagesc(overall_score(ens_sort)', [0 0.75]);
    colormap(example_ensemble_overall_score_ax, [linspace(0,1,100)' linspace(0,1,100)' linspace(0,1,100)']);
%    colormap(example_ensemble_overall_score_ax, jet(256));
    xlabel('all');
    axes(example_ensemble_overall_score_cbar_ax);
    imagesc((0:.01:1)', [0 1]);
    colormap(example_ensemble_overall_score_cbar_ax, [linspace(0,1,100)' linspace(0,1,100)' linspace(0,1,100)']);
    set(example_ensemble_overall_score_cbar_ax, 'TickDir', 'out', 'YTick', [], 'XTick',[]);

    axes(example_ensemble_w1_score_ax);
    imagesc(w1_score(ens_sort)', [0 0.75]);
    colormap(example_ensemble_w1_score_ax, [linspace(0,1,100)' linspace(0,0,100)' linspace(0,0,100)']);
    xlabel('w1');
    axes(example_ensemble_w1_score_cbar_ax);
    imagesc((0:.01:1)', [0 1]);
    colormap(example_ensemble_w1_score_cbar_ax, [linspace(0,1,100)' linspace(0,0,100)' linspace(0,0,100)']);
    set(example_ensemble_w1_score_cbar_ax, 'TickDir', 'out', 'YTick', [], 'XTick',[]);

    axes(example_ensemble_w1p_score_ax);
    imagesc(w1p_score(ens_sort)', [0 0.5]);
    colormap(example_ensemble_w1p_score_ax, [linspace(0,0,100)' linspace(0,1,100)' linspace(0,1,100)']);
    xlabel('Pro');
    axes(example_ensemble_w1p_score_cbar_ax);
    imagesc((0:.01:1)', [0 1]);
    colormap(example_ensemble_w1p_score_cbar_ax, [linspace(0,0,100)' linspace(0,1,100)' linspace(0,1,100)']);
    set(example_ensemble_w1p_score_cbar_ax, 'TickDir', 'out', 'YTick', [], 'XTick',[]);

    axes(example_ensemble_w1r_score_ax);
    imagesc(w1r_score(ens_sort)', [0 1]);
    colormap(example_ensemble_w1r_score_ax, [linspace(0,1,100)' linspace(0,0.5,100)' linspace(0,0,100)']);
    xlabel('Ret');
    axes(example_ensemble_w1r_score_cbar_ax);
    imagesc((0:.01:1)', [0 1]);
    colormap(example_ensemble_w1r_score_cbar_ax, [linspace(0,1,100)' linspace(0,0.5,100)' linspace(0,0,100)']);
    set(example_ensemble_w1r_score_cbar_ax, 'TickDir', 'out', 'YTick', [], 'XTick',[]);

 
    axes(example_ensemble_w2_score_ax);
    imagesc(w2_score(ens_sort)', [0 0.75]);
    colormap(example_ensemble_w2_score_ax, [linspace(0,0,100)' linspace(0,0.5,100)' linspace(0,1,100)']);
    xlabel('w2');
    axes(example_ensemble_w2_score_cbar_ax);
    imagesc((0:.01:1)', [0 1]);
    colormap(example_ensemble_w2_score_cbar_ax, [linspace(0,0,100)' linspace(0,0.5,100)' linspace(0,1,100)']);
    set(example_ensemble_w2_score_cbar_ax, 'TickDir', 'out', 'YTick', [], 'XTick',[]);

    axes(example_ensemble_w2p_score_ax);
    imagesc(w2p_score(ens_sort)', [0 1]);
    colormap(example_ensemble_w2p_score_ax, [linspace(0,0,100)' linspace(0,1,100)' linspace(0,1,100)']);
    xlabel('Pro');

    axes(example_ensemble_w2r_score_ax);
    imagesc(w2r_score(ens_sort)', [0 .5]);
    colormap(example_ensemble_w2r_score_ax, [linspace(0,1,100)' linspace(0,0.5,100)' linspace(0,0,100)']);
    xlabel('Ret');

    set(example_ensemble_overall_score_ax, 'TickDir', 'out', 'FontSize', fsize, 'XTick',[],'YTick',[]);
    set(example_ensemble_w1_score_ax, 'TickDir', 'out', 'FontSize', fsize, 'XTick',[],'YTick',[]);
    set(example_ensemble_w1p_score_ax, 'TickDir', 'out', 'FontSize', fsize, 'XTick',[],'YTick',[]);
    set(example_ensemble_w1r_score_ax, 'TickDir', 'out', 'FontSize', fsize, 'XTick',[],'YTick',[]);
    set(example_ensemble_w2_score_ax, 'TickDir', 'out', 'FontSize', fsize, 'XTick',[],'YTick',[]);
    set(example_ensemble_w2p_score_ax, 'TickDir', 'out', 'FontSize', fsize, 'XTick',[],'YTick',[]);
    set(example_ensemble_w2r_score_ax, 'TickDir', 'out', 'FontSize', fsize, 'XTick',[],'YTick',[]);
    
    % --- ensemble statistics
    ens_count = zeros(length(ens_dat.anim), 3);
    frac_cells_in_any_ens = zeros(length(ens_dat.anim), 3);
    ens_overlap_fraction = zeros(length(ens_dat.anim), 3);
    ens_size_mu = nan*zeros(length(ens_dat.anim), 3);
    for ani=1:length(ens_dat.anim)
        l_ids = {dat.anim_data(ani).types_by_id.L2, dat.anim_data(ani).types_by_id.L3, dat.anim_data(ani).types_by_id.L4};
    
        for li=1:length(l_tag)
            svi = subvol_idx(ani, li);
            
            ens_member_mat = ens_dat.anim(ani).subvol(svi).all_corr.ens_matrix;
            idx_layer = find(ismember(ens_dat.anim(ani).subvol(svi).all_corr.ids, l_ids{li}));
            ens_member_mat = ens_member_mat(:,idx_layer);
            
            ens_count(ani,li) = size(ens_member_mat,1);

            n_cells_subvol = length(find(dat.anim_data(ani).subvol_idx == svi & ismember(dat.anim_data(ani).ids, l_ids{li})));
            n_cells_in_any_ens = length(find(sum(ens_member_mat) > 0));
            frac_cells_in_any_ens(ani,li) = n_cells_in_any_ens / n_cells_subvol;

            % ensemble size
%            ens_size_mu(ani, li) = nanmean(sum(ens_dat.anim(ani).subvol(svi).all_corr.ens_matrix'));
            ens_size_mu(ani, li) = nanmean(sum(ens_member_mat'));

            % loop thru ensembles to compute overlap
            n_ens = size(ens_member_mat,1);
            overlap = [];
            if (n_ens > 1)
                for e1=1:n_ens-1
                    e1i = find(ens_member_mat(e1,:));
                    for e2=e1+1:n_ens
                        e2i = find(ens_member_mat(e2,:));
                        t_olp = length(intersect(e1i,e2i))/length(union(e1i,e2i)); % e1 = e2 yields 1 ; no overlap yields 0
                        overlap = [overlap t_olp];
                    end
                end
            end
            ens_overlap_fraction(ani,li) = nanmean(overlap);
        end
    end

    ens_size = nan*zeros(length(ens_dat.anim),3,max(ens_count(:))); % animal X layer X ensemble ; nan means this ensemble not present in this an/layer
    for ani=1:length(ens_dat.anim)
        for li=1:length(l_tag)
            svi = subvol_idx(ani, li);
            ens_member_mat = ens_dat.anim(ani).subvol(svi).all_corr.ens_matrix;
            idx_layer = find(ismember(ens_dat.anim(ani).subvol(svi).all_corr.ids, l_ids{li}));
            ens_member_mat = ens_member_mat(:,idx_layer); 
            if (~isempty(ens_member_mat)) 
                for u=1:size(ens_member_mat,1)-1
                    ens_size(ani,li,u) = length(find(ens_member_mat(u,:)));
                end
            end
        end
    end

    % whisker type ens membership
    t_ens_membership_frac_by_type = zeros(length(ens_dat.anim), 3, 3); % anim X layer X type "for each cell type, what fraction of those cells are in T ens?"
    t_ens_composition_frac_by_type = zeros(length(ens_dat.anim), 3, 4); % anim X layer X type (fourth type: non-touch cells)
    ens_ids = get_ensemble_ids;
    for ani=1:length(ens_dat.anim)
        l_ids = {dat.anim_data(ani).types_by_id.L2, dat.anim_data(ani).types_by_id.L3, dat.anim_data(ani).types_by_id.L4};
        for li=1:length(l_tag)
            svi = subvol_idx(ani, li);

            ens_member_mat = ens_dat.anim(ani).subvol(svi).all_corr.ens_matrix;
            ensmi = find(sum(ens_member_mat) > 0);
            ens_id = ens_ids(ani).touch_ensemble_ids;
            ens_id = intersect(l_ids{li}, ens_id);

            mw_id =  dat.anim_data(ani).types_by_id.mwCellIds;
            uni_id =  dat.anim_data(ani).types_by_id.unidirCellIds;
            bi_id =  dat.anim_data(ani).types_by_id.bidirCellIds;
            non_t_id = setdiff(dat.anim_data(ani).ids, union(bi_id, union(mw_id, uni_id)));

            fn = {'mw_id','uni_id','bi_id','non_t_id'};
            for f=1:length(fn)
                eval([fn{f} ' = intersect(' fn{f} ', l_ids{li});']);
            end

            t_ens_membership_frac_by_type(ani, li, 1) = length(intersect(ens_id, mw_id))/length(mw_id);
            t_ens_membership_frac_by_type(ani, li, 2) = length(intersect(ens_id, uni_id))/length(uni_id);
            t_ens_membership_frac_by_type(ani, li, 3) = length(intersect(ens_id, bi_id))/length(bi_id);

            t_ens_composition_frac_by_type(ani, li, 1) = length(intersect(ens_id, mw_id))/length(ens_id);
            t_ens_composition_frac_by_type(ani, li, 2) = length(intersect(ens_id, uni_id))/length(ens_id);
            t_ens_composition_frac_by_type(ani, li, 3) = length(intersect(ens_id, bi_id))/length(ens_id);
            t_ens_composition_frac_by_type(ani, li, 4) = length(intersect(ens_id, non_t_id))/length(ens_id);
        end
    end

    % summary bar plots:
    marker_col = [.5 .5 .5];
    bar_col = [1 1 1]*0.8;    

    bar_col = dat.settings.colors.mwColor;     
    bar_3l_plot(summ_frac_in_ens_by_type_ax(1), dat, squeeze(t_ens_composition_frac_by_type(:,:,1)), 'Frac. T ens. MW',bar_col, marker_col);
    bar_col = dat.settings.colors.bidiColor
    bar_3l_plot(summ_frac_in_ens_by_type_ax(2), dat, squeeze(t_ens_composition_frac_by_type(:,:,3)), 'Frac. T ens. bSW',bar_col, marker_col);
    bar_col =  dat.settings.colors.uniColor;
    bar_3l_plot(summ_frac_in_ens_by_type_ax(3), dat, squeeze(t_ens_composition_frac_by_type(:,:,2)), 'Frac. T ens. uSW',bar_col, marker_col);
    for l=1:3
        type_str = {'MW','uSW','bSW','nonT'};
        plotted_val = squeeze(t_ens_composition_frac_by_type(:,l,:));
        [pval_anov zz anov_stats] = anova1(plotted_val, [], 'off');
        disp(['ANOVA pval for ens composition ' l_tag{l} ' : ' num2str(pval_anov)]);
        if (pval_anov < 0.05)
            % multiple comparison post-hoc test
            c = multcompare(anov_stats, 'CType', 'hsd','Display','off');
            for v=1:size(c,1);
                pval_sr = signrank(plotted_val(:,c(v,1)), plotted_val(:,c(v,2)));
                [h_tt pval_tt] = ttest(plotted_val(:,c(v,1)), plotted_val(:,c(v,2)));
                if(c(v,6) < 0.05)
                    disp(sprintf(' ** %s vs. %s p: %0.3f ; RAW paired ttest pvalue: %0.3f', type_str{c(v,1)}, type_str{c(v,2)}, c(v,6), pval_tt));
                else
                    disp(sprintf(' %s vs. %s p: %0.3f ; RAW paired ttest pvalue %0.3f', type_str{c(v,1)}, type_str{c(v,2)}, c(v,6), pval_tt));
                end
            end
        end
    end

    xmax = 0;
    for t=1:3
        aa = axis(summ_frac_in_ens_by_type_ax(t));
        xmax = max(xmax, aa(2));
    end
    for t=1:3
        aa = axis(summ_frac_in_ens_by_type_ax(t));
        axis(summ_frac_in_ens_by_type_ax(t),[aa(1) xmax aa(3) aa(4)]);
    end

    % within-layer comparison (ANOVA)

    % --- touch ensemble statistics 
    min_touch_cell_frac_touch_ensemble = 0.5; % @ least 50%
        
    % who are the tuch ensembles? for each ensemble, frac T cells
    frac_t_cells_per_layer = {[],[],[]}; % for aggregate histograms across animals 
    frac_t_cells_per_layer_ani_idx = {[],[],[]}; % for aggregate histograms across animals 

    t_ens_member_mat_all = {}; % ani X li
    t_ens_count = nan*zeros(length(ens_dat.anim), 3);
    t_ens_size_mu = nan*zeros(length(ens_dat.anim), 3);
    t_ens_overlap_fraction = zeros(length(ens_dat.anim), 3);
    frac_cells_in_any_t_ens = nan*zeros(length(ens_dat.anim), 3);

    non_t_ens_member_mat_all = {}; % ani X li
    non_t_ens_count = nan*zeros(length(ens_dat.anim), 3);
    non_t_ens_size_mu = nan*zeros(length(ens_dat.anim), 3);
    non_t_ens_overlap_fraction = zeros(length(ens_dat.anim), 3);
    frac_cells_in_any_non_t_ens = nan*zeros(length(ens_dat.anim), 3);

    for ani=1:length(ens_dat.anim)
        l_ids = {dat.anim_data(ani).types_by_id.L2, dat.anim_data(ani).types_by_id.L3, dat.anim_data(ani).types_by_id.L4};
        t_ids = union(dat.anim_data(ani).types_by_id.mwCellIds, union(dat.anim_data(ani).types_by_id.unidirCellIds,dat.anim_data(ani).types_by_id.bidirCellIds));
    
        for li=1:length(l_tag)
            svi = subvol_idx(ani, li);
            
            ens_member_mat = ens_dat.anim(ani).subvol(svi).all_corr.ens_matrix;
            idx_touch = find(ismember(ens_dat.anim(ani).subvol(svi).all_corr.ids, t_ids));

            touch_ensi = [];
            for e=1:size(ens_member_mat,1)
                idx_ens = find(ens_member_mat(e,:));
                frac_t_cells_per_layer{li}(end+1) = length(intersect(idx_touch,idx_ens))/length(idx_ens);
                frac_t_cells_per_layer_ani_idx{li}(end+1) = ani;

                if (frac_t_cells_per_layer{li}(end) >= min_touch_cell_frac_touch_ensemble)
                    touch_ensi = [touch_ensi e];
                end
            end

            % store for now n later
            t_ens_member_mat_all{ani,li} = ens_member_mat(touch_ensi,:);
            non_touch_ensi = setdiff(1:size(ens_member_mat,1), touch_ensi);
            non_t_ens_member_mat_all{ani,li} = ens_member_mat(non_touch_ensi,:);

            % now get stats
            idx_layer = find(ismember(ens_dat.anim(ani).subvol(svi).all_corr.ids, l_ids{li}));
            t_ens_member_mat = t_ens_member_mat_all{ani,li}(:,idx_layer);
            non_t_ens_member_mat = non_t_ens_member_mat_all{ani,li}(:,idx_layer);
            
            t_ens_count(ani,li) = size(t_ens_member_mat,1);

            n_cells_subvol = length(find(dat.anim_data(ani).subvol_idx == svi & ismember(dat.anim_data(ani).ids, l_ids{li})));
            n_cells_in_any_t_ens = length(find(sum(t_ens_member_mat) > 0));
            frac_cells_in_any_t_ens(ani,li) = n_cells_in_any_t_ens / n_cells_subvol;

            t_ens_size_mu(ani, li) = nanmean(sum(t_ens_member_mat'));
            
            non_t_ens_count(ani,li) = size(non_t_ens_member_mat,1);

            n_cells_in_any_non_t_ens = length(find(sum(non_t_ens_member_mat) > 0));
            frac_cells_in_any_non_t_ens(ani,li) = n_cells_in_any_non_t_ens / n_cells_subvol;

            non_t_ens_size_mu(ani, li) = nanmean(sum(non_t_ens_member_mat'));      


            % loop thru ensembles to compute overlap
            n_ens = size(t_ens_member_mat,1);
            overlap = [];
            if (n_ens > 1)
                for e1=1:n_ens-1
                    e1i = find(t_ens_member_mat(e1,:));
                    for e2=e1+1:n_ens
                        e2i = find(t_ens_member_mat(e2,:));
                        t_olp = length(intersect(e1i,e2i))/length(union(e1i,e2i)); % e1 = e2 yields 1 ; no overlap yields 0
                        overlap = [overlap t_olp];
                    end
                end
            end
            t_ens_overlap_fraction(ani,li) = nanmean(overlap);
           
            n_ens = size(non_t_ens_member_mat,1);
            overlap = [];
            if (n_ens > 1)
                for e1=1:n_ens-1
                    e1i = find(non_t_ens_member_mat(e1,:));
                    for e2=e1+1:n_ens
                        e2i = find(non_t_ens_member_mat(e2,:));
                        t_olp = length(intersect(e1i,e2i))/length(union(e1i,e2i)); % e1 = e2 yields 1 ; no overlap yields 0
                        overlap = [overlap t_olp];
                    end
                end
            end
            non_t_ens_overlap_fraction(ani,li) = nanmean(overlap);
 
        end
    end

    t_ens_size_mu(find(isnan(t_ens_size_mu))) = 0;
    non_t_ens_size_mu(find(isnan(non_t_ens_size_mu))) = 0;
    frac_cells_in_any_t_ens(find(isnan(frac_cells_in_any_t_ens))) = 0;
    frac_cells_in_any_non_t_ens(find(isnan(frac_cells_in_any_non_t_ens))) = 0;
    t_ens_overlap_fraction(find(isnan(t_ens_overlap_fraction))) = 0;
    non_t_ens_overlap_fraction(find(isnan(non_t_ens_overlap_fraction))) = 0;

    bar_3l_plot(summ_frac_in_t_ens_ax, dat, frac_cells_in_any_t_ens, 'Fraction of cells in touch ensemble', bar_col, marker_col);
    bar_3l_plot(summ_t_ens_size_ax, dat, t_ens_size_mu, 'Mean touch ensemble size',bar_col, marker_col);
    bar_3l_plot(summ_t_ens_count_ax, dat, t_ens_count, 'Mean touch ensemble count',bar_col, marker_col);
    bar_3l_plot(summ_overlap_t_ens_ax, dat, t_ens_overlap_fraction, 'Mean touch ensemble overlap',bar_col, marker_col);

    bar_3l_plot(summ_frac_in_ens_ax, dat, frac_cells_in_any_non_t_ens, 'Fraction of cells in nonT ensemble', bar_col, marker_col);
    bar_3l_plot(summ_ens_size_ax, dat, non_t_ens_size_mu, 'Mean nonT ensemble size',bar_col, marker_col);
    bar_3l_plot(summ_ens_count_ax, dat, non_t_ens_count, 'Mean nonT ensemble count',bar_col, marker_col);
    bar_3l_plot(summ_overlap_non_t_ens_ax, dat, non_t_ens_overlap_fraction, 'Mean non touch ensemble overlap',bar_col, marker_col);
    for l=1:3 
        [h pval] = ttest(frac_cells_in_any_t_ens(:,l), frac_cells_in_any_non_t_ens(:,l));
        disp(sprintf('%s T vs. not T frac in ens p-value ttest: %0.3f', l_tag{l}, pval));
    end         
    for l=1:3 
        [h pval] = ttest(t_ens_count(:,l), non_t_ens_count(:,l));
        disp(sprintf('%s T vs. not T ens count p-value ttest: %0.3f', l_tag{l}, pval));
    end         
    for l=1:3 
        [h pval] = ttest(t_ens_size_mu(:,l), non_t_ens_size_mu(:,l));
        disp(sprintf('%s T vs. not T ens cell count p-value ttest: %0.3f', l_tag{l}, pval));
    end        
    for l=1:3 
        [h pval] = ttest(t_ens_overlap_fraction(:,l), non_t_ens_overlap_fraction(:,l));
        disp(sprintf('%s T vs. not T ens overlap p-value ttest: %0.3f', l_tag{l}, pval));
    end   

    % histograms
    w = 0.08;
    for li=1:length(l_tag)
        axes(summ_hist_frac_in_t_ens_by_type_ax(li));
        [n x] = hist(frac_t_cells_per_layer{li},0:.1:1);;
        for ii=1:length(x)
            rectangle('Position',[x(ii)-w/2 0 w n(ii)], 'EdgeColor','none', 'FaceColor', [1 1 1]*0.8);
        end
        if (li == 3)
            set(summ_hist_frac_in_t_ens_by_type_ax(li),'TickDir','out','XTick',0:.25:1, 'FontSize', fsize);
            xlabel(summ_hist_frac_in_t_ens_by_type_ax(li), 'Fraction touch neurons');
        else
            set(summ_hist_frac_in_t_ens_by_type_ax(li),'TickDir','out','XTick',[], 'FontSize', fsize);
        end   
        ylabel(summ_hist_frac_in_t_ens_by_type_ax(li), 'Ensemble count');
    end


    % --- frac of ca events related to touch that comes from ensemble ; prob that ens cell will respond vs. non-ens

    % gather the data - in/out ensemble touches for all cell types, layers
    te_dat = get_two_whisker_touch_ca_ev_data;
    event_frac_in_t_ens = zeros(length(ens_dat.anim), length(l_tag)); % animal X layer
    prob_resp_t_ens_member = zeros(length(ens_dat.anim), length(l_tag)); % animal X layer
    prob_resp_t_ens_nonmember = zeros(length(ens_dat.anim), length(l_tag)); % animal X layer
    n_cells_multi = 10;
    prob_resp_frac_trial_thresh_multi = 0.1;
    n_repetitions = 1000; % how many times we pick n_cells_multi from the touch ensembles or from non touch ensemble pool
    prob_resp_multi_t_ens_member = nan*zeros(length(ens_dat.anim), length(l_tag)); % animal X layer
    prob_resp_multi_t_ens_nonmember = nan*zeros(length(ens_dat.anim), length(l_tag)); % animal X layer
    pairwise_nt_corr_t_cell_t_ens_member = nan*zeros(length(ens_dat.anim), length(l_tag)); % animal X layer
    pairwise_nt_corr_t_cell_t_ens_nonmember = nan*zeros(length(ens_dat.anim), length(l_tag)); % animal X layer
    for ani=1:length(ens_dat.anim)
        for li=1:length(l_tag)
            l_ids = {dat.anim_data(ani).types_by_id.L2, dat.anim_data(ani).types_by_id.L3, dat.anim_data(ani).types_by_id.L4};
            svi = subvol_idx(ani, li);

            t_ens_member_mat = t_ens_member_mat_all{ani,li};

            t_ensmi = find(sum(t_ens_member_mat) > 0);
            t_ens_id = ens_dat.anim(ani).subvol(svi).all_corr.ids(t_ensmi);
            t_ens_id = intersect(l_ids{li}, t_ens_id);
            non_t_ens_id = setdiff(intersect(l_ids{li}, dat.anim_data(ani).ids(find(dat.anim_data(ani).subvol_idx == svi))), t_ens_id);
            all_t_id = union(dat.anim_data(ani).types_by_id.mwCellIds, dat.anim_data(ani).types_by_id.swCellIds);
            non_t_ens_touch_cell_id = intersect(non_t_ens_id, all_t_id);

            touch_ev_dat = te_dat(ani).sv_dat{svi};
            binary_response_by_trial = 0*touch_ev_dat.ev_per_trial;
            binary_response_by_trial(find(touch_ev_dat.ev_per_trial > 0)) = 1;

            all_touchi = find(touch_ev_dat.trial_type_vec);

            ini = find(ismember(touch_ev_dat.cell_ids, t_ens_id));
            outi = find(ismember(touch_ev_dat.cell_ids, non_t_ens_touch_cell_id));

            % pairwise correlations during no touch epoch
            cini = find(ismember(ens_dat.anim(ani).subvol(svi).all_corr.ids, t_ens_id));
            couti = find(ismember(ens_dat.anim(ani).subvol(svi).all_corr.ids, non_t_ens_touch_cell_id));
            pairwise_nt_corr_t_cell_t_ens_member(ani,li) = nanmean(reshape(ens_dat.anim(ani).subvol(svi).all_corr.M_original_notouch(cini,cini),[],1));
            pairwise_nt_corr_t_cell_t_ens_nonmember(ani,li) = nanmean(reshape(ens_dat.anim(ani).subvol(svi).all_corr.M_original_notouch(couti,couti),[],1));

            % fraction of events carried by specific cell types:
            touch_ev_count_per_cell = nansum(binary_response_by_trial(:,all_touchi)');
            ept = binary_response_by_trial(:,all_touchi);
            ept(find(ept > 1)) = 1;
            touch_ev_prob_per_cell = nansum(ept')/size(ept,2);
            touch_ev_count_population = sum(touch_ev_count_per_cell);

            event_frac_in_t_ens(ani,li) = sum(touch_ev_count_per_cell(ini))/touch_ev_count_population;

            if (length(ini) > 0)
                prob_resp_t_ens_member(ani, li) = nanmean(touch_ev_prob_per_cell(ini));
            end
            if (length(outi) > 0)
                prob_resp_t_ens_nonmember(ani, li) = nanmean(touch_ev_prob_per_cell(outi));
            end

            % population
            if (length(ini) > n_cells_multi & length(outi) > n_cells_multi)
                n_trials = length(all_touchi);
                prob_resp_this_subvol_t_ens = nan*zeros(1,n_repetitions);
                prob_resp_this_subvol_non_t_ens = nan*zeros(1,n_repetitions);
                for n=1:n_repetitions
                    i_t_ens = ini(randperm(length(ini), n_cells_multi));
                    resp_frac_per_trial = sum(binary_response_by_trial (i_t_ens, all_touchi))/length(i_t_ens);
                    prob_resp_this_subvol_t_ens(n) = length(find(resp_frac_per_trial >= prob_resp_frac_trial_thresh_multi))/n_trials;

                    o_t_ens = outi(randperm(length(outi), n_cells_multi));
                    resp_frac_per_trial = sum(binary_response_by_trial (o_t_ens, all_touchi))/length(o_t_ens);
                    prob_resp_this_subvol_non_t_ens(n) = length(find(resp_frac_per_trial > prob_resp_frac_trial_thresh_multi))/n_trials;                    
                end

                prob_resp_multi_t_ens_member (ani,li) = nanmedian(prob_resp_this_subvol_t_ens);
                prob_resp_multi_t_ens_nonmember (ani,li) = nanmedian(prob_resp_this_subvol_non_t_ens);
            end
        end
    end

    prob_resp_t_ens_member(find(isnan(prob_resp_t_ens_member))) = 0; % does not make sense, no ensemble means ZERO, not skip
    prob_resp_t_ens_nonmember(find(isnan(prob_resp_t_ens_nonmember))) = 0;
    prob_resp_multi_t_ens_member(find(isnan(prob_resp_multi_t_ens_member))) = 0; % does not make sense, no ensemble means ZERO, not skip
    prob_resp_multi_t_ens_nonmember(find(isnan(prob_resp_multi_t_ens_nonmember))) = 0;
 
    for l=1:3 
        [h pval] = ttest(prob_resp_t_ens_member(:,l), prob_resp_t_ens_nonmember(:,l));
        disp(sprintf('%s T ens vs. out T ens T cell prob resp 1 cell p-value ttest: %0.3f', l_tag{l}, pval));
    end  

    for l=1:3 
        [h pval] = ttest(prob_resp_multi_t_ens_member(:,l), prob_resp_multi_t_ens_nonmember(:,l));
        disp(sprintf('%s T ens vs. out T ens T cell prob resp 10 cells p-value ttest: %0.3f', l_tag{l}, pval));
    end  

    marker_col = [.5 .5 .5];
    bar_col = [1 1 1]*0.8;    
    bar_3l_plot(summ_event_frac_in_t_ens_ax, dat, event_frac_in_t_ens, 'Fraction of Ca events in touch ens', bar_col, marker_col);
    bar_3l_plot(summ_prob_resp_on_touch_t_ens_ax, dat, prob_resp_t_ens_member, 'P response, in touch ens', bar_col, marker_col);
    bar_3l_plot(summ_prob_resp_on_touch_non_t_ens_ax, dat, prob_resp_t_ens_nonmember, 'P response, out T ens T cell', bar_col, marker_col);
  
    axis(summ_prob_resp_on_touch_t_ens_ax, [0    0.3000   -0.5000    2.5000]);
    axis(summ_prob_resp_on_touch_non_t_ens_ax, [0    0.3000   -0.5000    2.5000]);

    marker_col = [.5 .5 .5];
    bar_col = [1 1 1]*0.8;    
    bar_3l_plot(summ_prob_resp_at_least_one_on_touch_t_ens_ax, dat, prob_resp_multi_t_ens_member, {'P resp, any one of 10','in touch ens'}, bar_col, marker_col);
    bar_3l_plot(summ_prob_resp_at_least_one_on_touch_non_t_ens_ax, dat, prob_resp_multi_t_ens_nonmember, {'P resp, any one of 10','out touch ens'}, bar_col, marker_col);
    axis(summ_prob_resp_at_least_one_on_touch_t_ens_ax, [0   1   -0.5000    2.5000]);
    axis(summ_prob_resp_at_least_one_on_touch_non_t_ens_ax, [0    1 -0.5000    2.5000]);

    % --- visualize a subvolume projected in 2D and with connecting lines for correlation, limited to high corrs ; give spatial stats
    if (~isempty (ens_dat.anim(ani).subvol(svi).all_corr.ens_matrix))
        ani = find(strcmp(corr_mat_ex_animal, {dat.settings.animals.name}));
        li = find(strcmp(corr_mat_ex_layer, l_tag));
        svi = subvol_idx(ani,li);

        % get the guys that are layer/subvolume correct
        l_ids = {dat.anim_data(ani).types_by_id.L2, dat.anim_data(ani).types_by_id.L3, dat.anim_data(ani).types_by_id.L4};
        ids_idx = find(ens_dat.anim(ani).subvol(svi).all_corr.ens_matrix(1,:));
        ids_plotted = ens_dat.anim(ani).subvol(svi).all_corr.ids(ids_idx); 
        vali = find(ismember(ids_plotted, l_ids{li}));
        ids_idx = ids_idx(vali);
        ids_plotted = ids_plotted(vali);

        % neurons will be sized/shaded based on correlation
        M = make_symm_nan_diag (ens_dat.anim(ani).subvol(svi).all_corr.M_original(ids_idx,ids_idx));
        score_vec = nanmean(M');
        cmap = [linspace(0.8,0.2,100)' linspace(.8,.2,100)' linspace(.8,.2,100)'  ];
        plot_others=0;
        plot_example_subvol_projected(example_ensemble_spatial_ax, dat, ani, svi, l_tag{li}, ids_plotted, score_vec, cmap, max_xy, plot_others);

        if (0) % set to 1 to plot additional ensembles 2 and 3
            plot_others=0;
            cmap = [linspace(0.5,0,100)' linspace(.5,0,100)' linspace(0.5,0,100)'  ];
            ids_idx = find(ens_dat.anim(ani).subvol(svi).all_corr.ens_matrix(2,:));
            ids_plotted = ens_dat.anim(ani).subvol(svi).all_corr.ids(ids_idx); 
            M = make_symm_nan_diag (ens_dat.anim(ani).subvol(svi).all_corr.M_original(ids_idx,ids_idx));
            score_vec = nanmean(M');
            plot_example_subvol_projected(example_ensemble_spatial_ax, dat, ani, svi, l_tag{li}, ids_plotted, score_vec, cmap, max_xy, plot_others);
            cmap = [linspace(.75,0,100)' .8*linspace(1,1,100)' .8*linspace(1,1,100)'  ];
            ids_idx = find(ens_dat.anim(ani).subvol(svi).all_corr.ens_matrix(3,:));
            ids_plotted = ens_dat.anim(ani).subvol(svi).all_corr.ids(ids_idx); 
            M = make_symm_nan_diag (ens_dat.anim(ani).subvol(svi).all_corr.M_original(ids_idx,ids_idx));
            score_vec = nanmean(M');
            plot_example_subvol_projected(example_ensemble_spatial_ax, dat, ani, svi, l_tag{li}, ids_plotted, score_vec, cmap, max_xy, plot_others);
        end
    end

    ens_dist_fname = [dat.settings.summary_data_root filesep 'ens_distance_data.mat'];
    if (~exist(ens_dist_fname))
        median_dist_in_ens = zeros(length(ens_dat.anim), length(l_tag)); % animal X layer
        median_dist_out_ens = zeros(length(ens_dat.anim), length(l_tag)); % animal X layer
        for ani=1:length(ens_dat.anim)
            l_ids = {dat.anim_data(ani).types_by_id.L2, dat.anim_data(ani).types_by_id.L3, dat.anim_data(ani).types_by_id.L4};

            xy = dat.anim_data(ani).xy_micron;
            z = dat.anim_data(ani).z_micron_corrected;
            xyz = [xy ; z]';

            d_mat = nan*zeros(length(dat.anim_data(ani).ids));
            d_mat_ids = dat.anim_data(ani).ids;
            for d1=1:length(d_mat_ids)
                for d2=d1+1:length(d_mat_ids)
                    d_mat(d1,d2) = sqrt( sum((xyz(d1,:) - xyz(d2,:)).^2));
                    d_mat(d2,d1) = d_mat(d1,d2);
                end
            end

            % gather spatial stats
            for li = 1:length(l_tag)
                svi = subvol_idx(ani, li);

                t_ens_member_mat = t_ens_member_mat_all{ani,li};

                t_ensmi = find(sum(t_ens_member_mat) > 0);
                t_ens_id = ens_dat.anim(ani).subvol(svi).all_corr.ids(t_ensmi);
                t_ens_id = intersect(l_ids{li}, t_ens_id);
                non_t_ens_id = setdiff(intersect(l_ids{li}, dat.anim_data(ani).ids(find(dat.anim_data(ani).subvol_idx == svi))), t_ens_id);

                ensi = find(ismember(d_mat_ids, t_ens_id));
                non_ensi = find(ismember(d_mat_ids, non_t_ens_id));

                disp(sprintf('a: %d s: %d ensi: %d non-ensi: %d ', ani, li, length(ensi), length(non_ensi)));

                median_dist_in_ens(ani,li) = nanmedian(reshape(d_mat(ensi,ensi), [], 1));
                median_dist_out_ens(ani,li) = nanmedian(reshape(d_mat(non_ensi,non_ensi), [], 1));
            end
        end

        save(ens_dist_fname, 'median_dist_in_ens', 'median_dist_out_ens');
    else
        load(ens_dist_fname);
    end
    

    marker_col = [.5 .5 .5];
    bar_col = [1 1 1]*0.8;    
    bar_3l_plot(summ_pairwise_dist_in_t_ens_ax, dat, median_dist_in_ens, {'Meidan pairwise distance','in touch ens'}, bar_col, marker_col);
    bar_3l_plot(summ_pairwise_dist_out_t_ens_ax, dat, median_dist_out_ens, {'Median pairwise distance','out touch ens'}, bar_col, marker_col);
    axis(summ_pairwise_dist_in_t_ens_ax, [0  400   -0.5000    2.5000]);
    axis(summ_pairwise_dist_out_t_ens_ax, [0  400 -0.5000    2.5000]);
 
    for l=1:3 
        [h pval] = ttest(median_dist_in_ens(:,l), median_dist_out_ens(:,l));
        disp(sprintf('%s T vs. out T ens distance p-value ttest: %0.3f', l_tag{l}, pval));
    end  

function bar_3l_plot(ax, dat, plotted_val, x_lab_str, bar_col, marker_col)
    symbols = {'^','o','v','d','s','p','h'};
    h = 0.8;
    fsize = 10;
    offs = linspace(-1*h/4, h/4, length(dat.anim_data));
    msize = 5;
    pval_thresh = 0.05;
    
    plotted_val(find(isinf(plotted_val))) = nan;

    axes(ax);
    hold on;
    bar_type = 'horizontal';
    N = size(plotted_val,1);
    if (strcmp(bar_type, 'horizontal'))
        rectangle('Position',[0 2-h/2 nanmean(plotted_val(:,1)) h], 'EdgeColor','none', 'FaceColor', bar_col);
        plot(nanmean(plotted_val(:,1))+[1 -1]*nanstd(plotted_val(:,1))/sqrt(N), [2 2],  'k-', 'LineWidth', 2);
        rectangle('Position',[0 1-h/2 nanmean(plotted_val(:,2)) h], 'EdgeColor','none', 'FaceColor', bar_col);
        plot( nanmean(plotted_val(:,2))+[1 -1]*nanstd(plotted_val(:,2))/sqrt(N), [1 1], 'k-', 'LineWidth', 2);
        rectangle('Position',[0 0-h/2 nanmean(plotted_val(:,3)) h], 'EdgeColor','none', 'FaceColor', bar_col); 
        plot( nanmean(plotted_val(:,3))+[1 -1]*nanstd(plotted_val(:,3))/sqrt(N), [0 0], 'k-', 'LineWidth', 2);
        % per animal markers to look at individual data pts
        for a=1:length(dat.anim_data)
            plot (plotted_val(a,1), 2+offs(a), symbols{a}, 'MarkerFaceColor',marker_col , ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);
            plot (plotted_val(a,2), 1+offs(a), symbols{a}, 'MarkerFaceColor',marker_col , ...
                        'Color',[1 1 1]*0.75, 'MarkerSize', msize);
            plot (plotted_val(a,3), 0+offs(a), symbols{a}, 'MarkerFaceColor',marker_col , ...
                        'Color',[1 1 1]*0.75, 'MarkerSize', msize);        
        end
        aa = axis(ax);
        axis(ax, [0 aa(2) -0.5 2.5]);
        set(ax,'TickDir','out', 'YTick',[0 1 2], 'YTickLabel',{'L4','L3','L2'}, 'FontSize', fsize);
        xlabel(x_lab_str);
    end

    % STATS - probably going to use normal assumption, but lets just put the nonparametric in just in case we change mind
    l_tag = {'L2','L3','L4'};
    [pval_anov zz anov_stats] = anova1(plotted_val, [], 'off');
    if (iscell(x_lab_str)) ; x_lab_str = [x_lab_str{1} ' ' x_lab_str{2}] ; end
    disp(['ANOVA pval for ' x_lab_str  ': ' num2str(pval_anov)]);
    if (pval_anov < pval_thresh) 
        % multiple comparison post-hoc test
        c = multcompare(anov_stats, 'CType', 'hsd','Display','off');
        for v=1:size(c,1);
            pval_sr = signrank(plotted_val(:,c(v,1)), plotted_val(:,c(v,2)));
            [h_tt pval_tt] = ttest(plotted_val(:,c(v,1)), plotted_val(:,c(v,2)));
            if(c(v,6) < pval_thresh)
                disp(sprintf(' ** %s vs. %s p: %0.3f ; RAW paired ttest pvalue: %0.3f', l_tag{c(v,1)}, l_tag{c(v,2)}, c(v,6), pval_tt));
            else
                disp(sprintf(' %s vs. %s p: %0.3f ; RAW paired ttest pvalue %0.3f', l_tag{c(v,1)}, l_tag{c(v,2)}, c(v,6), pval_tt));
            end
        end
    end

    [hz pval_23] = ttest(plotted_val(:,1), plotted_val(:,2));
    [hz pval_34] = ttest(plotted_val(:,3), plotted_val(:,2));
    disp(sprintf(' L4: %0.3f / %0.3f  L3: %0.3f / %0.3f  L2: %0.3f / %0.3f ', nanmean(plotted_val(:,3)), nanstd(plotted_val(:,3)), nanmean(plotted_val(:,2)), nanstd(plotted_val(:,2)), nanmean(plotted_val(:,1)), nanstd(plotted_val(:,1))));
    disp(sprintf(' L4 vs. L3 p: %0.3f ', pval_34));
    disp(sprintf(' L3 vs. L2 p: %0.3f ', pval_23));



function M = make_symm_nan_diag (M)
    for i1=1:size(M,1)
        M(i1,i1) = nan;
        for i2=i1+1:size(M,1)
            M(i2,i1) = M(i1,i2);
        end
    end

function plot_example_subvol_projected(ax, dat, ani, svi, l_str, ids_plotted, score_vec, cmap, max_xy, plot_others)
% this makes the plots where all cells in a layer are projected into a single plane

    scaling_factor = 1000; % increase in ball size
    color_other = [1 1 1]*0.5;
    size_other = 5;

    % stuff to use in loopz
    xy = dat.anim_data(ani).xy_micron;
    hold(ax, 'on');

    % indices
    idx_layer = find(ismember(dat.anim_data(ani).ids, dat.anim_data(ani).types_by_id.(l_str)));
    idx_plotted = find(ismember(dat.anim_data(ani).ids, ids_plotted));
    idx_other = setdiff(intersect(idx_layer, find(dat.anim_data(ani).subvol_idx == svi)), idx_plotted);

    % color assignment
    M = max(score_vec);
    m = min(score_vec);
    M = 0.5;
    m = 0.15;
    cols = zeros(3,length(idx_plotted));
    for i=1:length(idx_plotted)
        f = (score_vec(i)-m)/(M-m);
        f = max(1,round(f*size(cmap,1)));
        cols(:,i) = cmap(f,:);
    end
    
    % plotting
    if (plot_others)
        scatter(ax,xy(1,idx_other), xy(2,idx_other), size_other, color_other, 'filled');
    end
    scatter(ax,xy(1,idx_plotted), xy(2,idx_plotted), score_vec*scaling_factor, cols', 'filled');
    scatter(ax, [20 100],[700 700], [m M]*scaling_factor, cmap([1 end],:), 'filled');

    axis(ax, [0 max_xy 0 max_xy]);
    set(ax,'XTick',[],'YTick',[]);

function plot_raw_example (ani, svi, dat, ens_dat, ens_sess, tr, ens_sort, example_ensemble_dff_ax, example_ensemble_dkappa_w1_ax, example_ensemble_dkappa_w2_ax, example_ensemble_theta_w1_ax, fsize)
    ids = ens_dat.anim(ani).subvol(svi).all_corr.ids;
    s_ids = ens_sess.caTSA.ids;
    vm_idx = nan*ens_sort;
    for i=1:length(ens_sort)
        vm_idx(i) = find(s_ids == ids(ens_sort(i)));
    end

    dt_sec = mode(diff(ens_sess.caTSA.time))/1000;
    x_idx = round(tr(1)/dt_sec):round(tr(2)/dt_sec);

    axes(example_ensemble_dff_ax);
    submat = ens_sess.caTSA.dffTimeSeriesArray.valueMatrix(vm_idx,x_idx);
    imagesc(submat, [0 quantile(submat(:),.975)]);
    cm = parula(256);
    cm = flipud(gray(256));

    colormap(example_ensemble_dff_ax, cm);

    set(example_ensemble_dff_ax, 'TickDir', 'out', 'FontSize', fsize, 'YTick',[]);
    xlabel('time(s)');
    set(example_ensemble_dff_ax, 'XTick',[1 range(x_idx)/2 range(x_idx)], 'XTickLabel', {'0', sprintf('%0.3d', range(tr)/2), sprintf('%0.3d', tr(2)-tr(1))});

    % plot max dKappa
    lw = 1;
    cmap = [ linspace(1,0,100) linspace(0,0,100); linspace(0,0,100) linspace(0,1,100); linspace(1,0,100) linspace(0,1,100)]';
 
    axes(example_ensemble_dkappa_w1_ax);
    w1_dk = ens_sess.derivedDataTSA.getTimeSeriesById(20112).value(x_idx);
    plot(x_idx,w1_dk, '-', 'LineWidth', lw, 'Color', dat.settings.colors.w1Color);
    aa = axis(example_ensemble_dkappa_w1_ax);
    qv = quantile(abs(w1_dk(find(abs(w1_dk) > 0))), .9);
    axis(example_ensemble_dkappa_w1_ax, [x_idx(1) x_idx(end) -1*qv qv]);

    axes(example_ensemble_dkappa_w2_ax);
    w2_dk = ens_sess.derivedDataTSA.getTimeSeriesById(20113).value(x_idx);
    plot(x_idx,w2_dk, '-', 'LineWidth', lw, 'Color', dat.settings.colors.w2Color);
    aa = axis(example_ensemble_dkappa_w2_ax);
    qv = quantile(abs(w2_dk(find(abs(w2_dk) > 0))), .9);
    axis(example_ensemble_dkappa_w2_ax, [x_idx(1) x_idx(end) -1*qv qv]);

    axes(example_ensemble_theta_w1_ax);
    theta = ens_sess.derivedDataTSA.getTimeSeriesById(10000).value(x_idx);
    plot(x_idx,theta, '-', 'LineWidth', lw, 'Color', [0 0.8 0]);
    aa = axis(example_ensemble_dkappa_w2_ax);
    qv = quantile(abs(theta(find(abs(theta) > 0))), [0.01 .99]);
    qvd = diff(qv)/4;
    qv = qv + [-qvd qvd]; 
    axis(example_ensemble_theta_w1_ax, [x_idx(1) x_idx(end) qv(1) qv(2)]);

    set(example_ensemble_dkappa_w1_ax,'XTick',[],'YTick',[]); 
    ylabel(example_ensemble_dkappa_w1_ax, 'w1');
    set(example_ensemble_dkappa_w2_ax,'XTick',[],'YTick',[]); 
    ylabel(example_ensemble_dkappa_w2_ax, 'w2');
    set(example_ensemble_theta_w1_ax,'XTick',[],'YTick',[]); 
    ylabel(example_ensemble_theta_w1_ax, 'setpoint');

