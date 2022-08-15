%
% a) example corr mats for one animal, subvolume
% b) distro comparison for all auto-corrs, one animal
% c) all 6 comparisons for all mice
% d) all 6 comparisons, all mice non-stim epoch
% e) distance dependence for 3 self-types
% f) distance at key position
% g) lambdas for distance-v-corr
%
function figure_correlations_basic (force_redo)
    if (nargin == 0) ; force_redo = 0; end % force_redo will regather data for two_whisker_corr_summary_data.mat

    params.force_redo = 0;
    params.exclude_matrix_name = 'excluded_for_correlations';
    dat = get_two_whisker_data(params);
    pt_dat = get_two_whisker_per_trial_data(force_redo);
    corr_dat_filename = [dat.settings.summary_data_root filesep 'two_whisker_corr_summary_data.mat'];

    % general settings
    max_xy = 700; % size of plane

    % w1/w2 or pro/ret
    if (1) % w1/w2
        mat_list = {'mw_mw','usw_usw','bsw_bsw'};
        ids_list = {'mw','usw','bsw'};
        mat_col= {dat.settings.colors.mwColor dat.settings.colors.uniColor dat.settings.colors.bidiColor};
        mat_pcol = {[.8 .6 .8], [.8 .8 .8], [.5 .75 1]};
        corr_labels = {'MW','uSW','bSW','uSW-MW','bSW-MW','uSW-bSW'};
    end

    %  --- setup figures
    fsize = 10;

    fh_main = figure ('Position',[0 0 1800 900]);
    for i=1:3
        example_trace_corr_mat_ax(i) = axes(fh_main,'Position', [.05 .95-0.125*i .05 .1]);
        trace_corr_hist_ax(i) = axes(fh_main, 'Position', [.125 .95-0.125*i .1 .1]);
        trace_corr_touch_stats_ax(i) = axes(fh_main, 'Position', [.25 .95-0.125*i .1 .1]);
        trace_corr_quiet_stats_ax(i) = axes(fh_main, 'Position', [.375 .95-0.125*i .1 .1]);

        trace_corr_dist_touch_ax(i) = axes(fh_main, 'Position', [.5 .95-0.125*i .1 .1]);
        trace_corr_norm_dist_touch_ax(i) = axes(fh_main, 'Position', [.65 .95-0.125*i .1 .1]);
        trace_corr_summ_bar_ax(i) = axes(fh_main, 'Position', [.8 .95-0.125*i .1 .1]);
    end
    example_corr_mat_cb_ax = axes(fh_main, 'Position', [0.05 .05 .05 .1]);

    % --- data gathering - first time around will take a bit
    % loop thru subvolumes and, by cell type and epoch, compute / gather correlations

    if (~exist(corr_dat_filename,'file') | force_redo)
        corr_dat = [];

        for a=1:length(dat.anim_data)
            xyz = [dat.anim_data(a).xy_micron; dat.anim_data(a).z_micron_corrected];

            % multi v single
            mw_ids = dat.anim_data(a).types_by_id.mwCellIds;
            bsw_ids = {dat.anim_data(a).types_by_id.w1BidirCellIds, dat.anim_data(a).types_by_id.w2BidirCellIds};
            usw_ids = {dat.anim_data(a).types_by_id.w1pCellIds, dat.anim_data(a).types_by_id.w1rCellIds, dat.anim_data(a).types_by_id.w2pCellIds, dat.anim_data(a).types_by_id.w2rCellIds};

            for s=1:length(pt_dat(a).per_subvol)
                % trace correlations
                fname_wc = [dat.settings.data_path_root filesep dat.settings.animals(a).name filesep 'session_neuropilone_final_merged' ...
                            filesep 'pairwiseCorrelations' filesep dat.settings.animals(a).name sprintf('_vol_%02d_sess__*.mat', s + dat.settings.animals(a).subvol_offset)];
                trace_corr_dat = get_trace_corrs(fname_wc, dat, pt_dat, a, s, mw_ids, usw_ids, bsw_ids);

                corr_dat.anim(a).subvol(s).trace_corrs_mw_mw = trace_corr_dat.corrs_mw_mw;
                corr_dat.anim(a).subvol(s).trace_corrs_usw_usw = trace_corr_dat.corrs_usw_usw;
                corr_dat.anim(a).subvol(s).trace_corrs_bsw_bsw = trace_corr_dat.corrs_bsw_bsw;

                corr_dat.anim(a).subvol(s).trace_corrs_usw_mw = trace_corr_dat.corrs_usw_mw;
                corr_dat.anim(a).subvol(s).trace_corrs_bsw_mw = trace_corr_dat.corrs_bsw_mw;
                corr_dat.anim(a).subvol(s).trace_corrs_usw_bsw = trace_corr_dat.corrs_usw_bsw;

                corr_dat.anim(a).subvol(s).trace_corrs_mw_ids = trace_corr_dat.mw_ids;
                corr_dat.anim(a).subvol(s).trace_corrs_usw_ids = trace_corr_dat.usw_ids;
                corr_dat.anim(a).subvol(s).trace_corrs_bsw_ids = trace_corr_dat.bsw_ids;
            end
        end

        save(corr_dat_filename, 'corr_dat');
    else % load it and globaliZE
        global corr_dat;
        if (isempty(corr_dat))
            load(corr_dat_filename);
        end
    end

    % additional processing - for each animal-layer combo, figure out which subvolume has most cells for that layer 
    [subvol_idx l_tag] = get_best_subvol_list (dat)
 
    % --- correlation histogram
    corr_mat_ex_animal = 'an280201';
    
    lw = 3;
    ani = find(strcmp(corr_mat_ex_animal,{dat.settings.animals.name}));
    epoch_used = 'allTime';
    l_ids = {dat.anim_data(ani).types_by_id.L2, dat.anim_data(ani).types_by_id.L3, dat.anim_data(ani).types_by_id.L4};
 
    corr_bins = -0.2:0.05:1;
    for l=1:3
        svi = subvol_idx(ani,l);  
        hold (trace_corr_hist_ax(l), 'on');
        touch_ids = [];

        % predone                     
        for m=1:3
            mat_name = ['trace_corrs_' mat_list{m}];
            ids_name = ['trace_corrs_' ids_list{m} '_ids'];
            tmat = corr_dat.anim(ani).subvol(svi).(mat_name).(epoch_used);
            ids = corr_dat.anim(ani).subvol(svi).(ids_name);
            touch_ids = [touch_ids ids];
            if (~iscell(ids)) ; ids = {ids} ; tmat={tmat}; end
            corr_vec = [];
            for i=1:length(ids)
                vali = find(ismember(ids{i}, l_ids{l}));
                stmat = tmat{i}(vali,vali); 
                for i=1:size(stmat,1) 
                    stmat(i,i) = nan;
                end

                corr_vec = [corr_vec stmat(:)'];
            end
            plot_fancy_hist (trace_corr_hist_ax(l), corr_vec, corr_bins, mat_col{m}, lw)
        end
    end

    % --- trace correlation during touch/no-touch stats
    plot_relative_corrs(trace_corr_touch_stats_ax, dat, corr_dat, 'allTime', subvol_idx, 'trace', mat_list, ids_list, corr_labels);
    plot_relative_corrs(trace_corr_quiet_stats_ax, dat, corr_dat, 'whiskerNoTouchTrials', subvol_idx, 'trace', mat_list, ids_list, corr_labels);
    title(trace_corr_touch_stats_ax(1), 'Trace correlations, all time');
    title(trace_corr_quiet_stats_ax(1), 'Trace correlations, no-touch epoch');
    for i=1:3
        ylabel(trace_corr_touch_stats_ax(i), l_tag{i});
        axis(trace_corr_touch_stats_ax(i), [0 7 0 0.325]);
        axis(trace_corr_quiet_stats_ax(i), [0 7 0 0.325]);
    end
 
function plot_relative_corrs(ax, dat, corr_dat, corrs_used, subvol_idx, mat_type, mat_list, ids_list, corr_labels)
    ms = 7;
    layer_data_mat = nan*zeros(3, length(dat.anim_data), 4);

    for l=1:3
        % aggregate across animals
        mu_mw = [];
        mu_usw = [];
        mu_bsw = [];
        mu_usw_mw = [];
        mu_bsw_mw = [];
        mu_usw_bsw = []; 
        for ani=1:length(dat.anim_data)
            l_ids = {dat.anim_data(ani).types_by_id.L2, dat.anim_data(ani).types_by_id.L3, dat.anim_data(ani).types_by_id.L4};
            svi = subvol_idx(ani,l);

            ids_mw = [mat_type '_corrs_' ids_list{1} '_ids'];
            vali_mw = find(ismember(corr_dat.anim(ani).subvol(svi).(ids_mw), l_ids{l}));
            M = corr_dat.anim(ani).subvol(svi).(sprintf('%s_corrs_%s', mat_type, mat_list{1})).(corrs_used);
            mu_mw(ani) = nanmean( reshape(M(vali_mw,vali_mw),[],1));


            ids_usw = [mat_type '_corrs_' ids_list{2} '_ids'];
            used_ids_all_usw = [];
            ids_all_usw = [];
            all_corrs_usw = [];
            for i=1:length(corr_dat.anim(ani).subvol(svi).(ids_usw))
                vali_usw = find(ismember(corr_dat.anim(ani).subvol(svi).(ids_usw){i}, l_ids{l}));
                used_ids_all_usw = [used_ids_all_usw corr_dat.anim(ani).subvol(svi).(ids_usw){i}(vali_usw)];
                ids_all_usw = [ids_all_usw corr_dat.anim(ani).subvol(svi).(ids_usw){i}];
                M = corr_dat.anim(ani).subvol(svi).(sprintf('%s_corrs_%s', mat_type, mat_list{2})).(corrs_used){i};
                all_corrs_usw = [all_corrs_usw M(:)'];
            end
            mu_usw(ani) = nanmean(all_corrs_usw);

            ids_bsw = [mat_type '_corrs_' ids_list{3} '_ids'];
            used_ids_all_bsw = [];
            ids_all_bsw = [];
            all_corrs_bsw = [];
            for i=1:length(corr_dat.anim(ani).subvol(svi).(ids_bsw))
                vali_bsw = find(ismember(corr_dat.anim(ani).subvol(svi).(ids_bsw){i}, l_ids{l}));
                used_ids_all_bsw = [used_ids_all_bsw corr_dat.anim(ani).subvol(svi).(ids_bsw){i}(vali_bsw)];
                ids_all_bsw = [ids_all_bsw corr_dat.anim(ani).subvol(svi).(ids_bsw){i}];
                M = corr_dat.anim(ani).subvol(svi).(sprintf('%s_corrs_%s', mat_type, mat_list{3})).(corrs_used){i};
                all_corrs_bsw = [all_corrs_bsw M(:)'];
            end
            mu_bsw(ani) = nanmean(all_corrs_bsw);

            vali_all_usw = find(ismember(ids_all_usw, used_ids_all_usw));
            vali_all_bsw = find(ismember(ids_all_bsw, used_ids_all_bsw));

            M = corr_dat.anim(ani).subvol(svi).(sprintf('%s_corrs_%s_%s', mat_type, ids_list{2}, ids_list{1})).(corrs_used);
            mu_usw_mw(ani) = nanmean( reshape(M(vali_all_usw,vali_mw),[],1)); 
            M = corr_dat.anim(ani).subvol(svi).(sprintf('%s_corrs_%s_%s', mat_type, ids_list{3}, ids_list{1})).(corrs_used);
            mu_bsw_mw(ani) = nanmean( reshape(M(vali_all_bsw,vali_mw),[],1));
            M = corr_dat.anim(ani).subvol(svi).(sprintf('%s_corrs_%s_%s', mat_type, ids_list{2}, ids_list{3})).(corrs_used);
            mu_usw_bsw(ani) = nanmean( reshape(M(vali_all_usw,vali_all_bsw),[],1)); 
        end

        hold(ax(l), 'on')
        symbols = {'^','o','v','d','s','p','h'};
        w = 0.6;
        offs = linspace(-w/3, w/3, length(dat.anim_data));
        y_base = 0;

        Nsq = sqrt(length(dat.anim_data));
        rectangle(ax(l), 'Position',[4-(w/2) y_base w nanmean(mu_mw)], 'EdgeColor','none', 'FaceColor', dat.settings.colors.mwColor);
        plot(ax(l), [1 1]*4, nanmean(mu_mw)+nanstd(mu_mw)*[1 -1]/Nsq, 'k-','LineWidth',1);
        rectangle(ax(l), 'Position',[2-(w/2) y_base w nanmean(mu_usw)], 'EdgeColor','none', 'FaceColor',  dat.settings.colors.uniColor);
        plot(ax(l), [1 1]*2, nanmean(mu_usw)+nanstd(mu_usw)*[1 -1]/Nsq, 'k-','LineWidth',1);
        rectangle(ax(l), 'Position',[3-(w/2) y_base w nanmean(mu_bsw)], 'EdgeColor','none', 'FaceColor',  dat.settings.colors.bidiColor);
        plot(ax(l), [1 1]*3, nanmean(mu_bsw)+nanstd(mu_bsw)*[1 -1]/Nsq, 'k-','LineWidth',1);

        for a=1:length(dat.anim_data)
            plot(ax(l), 4+offs(a), mu_mw(a), symbols{a}, 'MarkerFaceColor', [0 0 0 ], 'Color', [1 1 1]*.5, 'MarkerSize', ms); 
            plot(ax(l), 2+offs(a), mu_usw(a), symbols{a}, 'MarkerFaceColor', [1 1 1]*0, 'Color', [1 1 1]*.5, 'MarkerSize', ms); 
            plot(ax(l), 3+offs(a), mu_bsw(a), symbols{a}, 'MarkerFaceColor', [1 1 1]*0, 'Color',[1 1 1]*.5, 'MarkerSize', ms); 
        end

        aa = axis(ax(l));
        axis(ax(l) ,[0 7 aa(3) aa(4)])
        set(ax(l),'TickDir', 'out', 'XTick',[]);
        if (l == 3) 
            set(ax(l),'XTick', [1 2 3 4], 'XTickLabel', {'non T','uSW', 'bSW', 'MW'}, 'XTickLabelRotation',45);
        end

        %% STaTZ
        data_mat = nan*zeros(length(dat.anim_data), 4);
        data_mat(:,2) = mu_usw;
        data_mat(:,3) = mu_bsw;
        data_mat(:,4) = mu_mw;
        layer_data_mat(l,:,:) = data_mat;
        type_str = {'non T','uSW', 'bSW', 'MW'};
        type_used = 2:4;
        if (l == 3) ; type_used = [2 3]; end
        type_str = type_str(type_used);
        data_mat = data_mat(:,type_used);
        [pval_kw zz kw_stats] = kruskalwallis(data_mat, [], 'off');
        [pval_anov zz anov_stats] = anova1(data_mat, [], 'off');
        disp(['KW pval for ' corrs_used ' for L' num2str(l+1) ': ' num2str(pval_kw)  ' ANOVA: ' num2str(pval_anov)]);
        pval_thresh = 0.05;
        if (pval_anov < pval_thresh) 

            % multiple comparison test..
            c = multcompare(anov_stats, 'CType', 'hsd','Display','off');
            for v=1:size(c,1);
                pval_sr = signrank(data_mat(:,c(v,1)), data_mat(:,c(v,2)));
                [h_tt pval_tt] = ttest(squeeze(data_mat(:,c(v,1))), squeeze(data_mat(:,c(v,2))));
                if(c(v,6) < pval_thresh)
                    disp(sprintf(' ** %s vs. %s p: %0.3f ; RAW paired ttest pvalue: %0.3f', type_str{c(v,1)}, type_str{c(v,2)}, c(v,6), pval_tt));
                else
                    disp(sprintf(' %s vs. %s p: %0.3f ; RAW paired ttest pvalue %0.3f', type_str{c(v,1)}, type_str{c(v,2)}, c(v,6), pval_tt));
                end
            end
        end
    end
 
    type_str = {'non T','uSW', 'bSW', 'MW'};
    for m=1:4
        [h p43] = ttest(squeeze(layer_data_mat(3,:,m)),squeeze(layer_data_mat(2,:,m)));
        [h p32] = ttest(squeeze(layer_data_mat(1,:,m)),squeeze(layer_data_mat(2,:,m)));
        disp(sprintf(' %s L4 v. L3 p: %0.3f', type_str{m}, p43));
        disp(sprintf(' %s L3 v. L2 p: %0.3f', type_str{m}, p32));
    end

% sorts into n_groups, with a further sort by average value of group
function [M ski] = sort_kmeans(M, n_groups)
    ki = kmeans(M, n_groups);
    tM = M;
    for i=1:size(tM,1) ; tM(i,i) = nan; end
    for n=1:n_groups
        mu_val(n) = nanmean(reshape(tM(find(ki == n), find(ki == n)),[],1));
    end
    [irr vki] = sort(mu_val,'descend');

    nki = ki*nan;
    for n=1:n_groups
        ni = find(ki == vki(n));
        nki(ni) = n;
    end
    [irr ski] = sort (nki);
    M = M(ski,ski);

function trace_corr_dat = get_trace_corrs(fname_wc, dat, pt_dat, a, s, mw_ids, usw_ids, bsw_ids)
    % for referencing within p_dat per-trial dff matrices
    si = find(pt_dat(a).subvol == s);
    subvol_ids = pt_dat(a).ids(si);

    mw_idx = find(ismember(subvol_ids, mw_ids));
    for i=1:length(usw_ids)
        usw_idx{i} = find(ismember(subvol_ids, usw_ids{i}));
    end
    for i=1:length(bsw_ids)
        bsw_idx{i} = find(ismember(subvol_ids, bsw_ids{i}));
    end

    % file loop - this will be for ONE subvolume, but handles (3 usually) possible correlation types
    fl = dir(fname_wc);
    base_dir = fileparts_crossplatform(fname_wc);
    for f=1:length(fl)
        disp(['Processing ' base_dir filesep fl(f).name]);
        C = load([base_dir filesep fl(f).name]);
        i_tag = strfind(fl(f).name,'__');
        tag = fl(f).name(i_tag+2:end-4);
       
        mw_idx = find(ismember(C.dat.cellIds, mw_ids));
        trace_corr_dat.corrs_mw_mw.(tag) = C.dat.corrMat(mw_idx,mw_idx);
        trace_corr_dat.mw_ids = C.dat.cellIds(mw_idx);

        all_usw_idx = [];
        for i=1:length(usw_ids)
            usw_idx = find(ismember(C.dat.cellIds, usw_ids{i}));
            all_usw_idx = [all_usw_idx usw_idx];
            trace_corr_dat.corrs_usw_usw.(tag){i}  = C.dat.corrMat(usw_idx,usw_idx);
            trace_corr_dat.usw_ids{i} = C.dat.cellIds(usw_idx);
        end

        all_bsw_idx = [];
        for i=1:length(bsw_ids)
            bsw_idx = find(ismember(C.dat.cellIds, bsw_ids{i}));
            all_bsw_idx = [all_bsw_idx bsw_idx];
            trace_corr_dat.corrs_bsw_bsw.(tag){i}  = C.dat.corrMat(bsw_idx,bsw_idx);
            trace_corr_dat.bsw_ids{i} = C.dat.cellIds(bsw_idx);
        end

        trace_corr_dat.corrs_usw_mw.(tag)  = C.dat.corrMat(all_usw_idx,mw_idx);
        trace_corr_dat.corrs_bsw_mw.(tag)  = C.dat.corrMat(all_bsw_idx,mw_idx);
        trace_corr_dat.corrs_usw_bsw.(tag)  = C.dat.corrMat(all_usw_idx,all_bsw_idx);
    end

