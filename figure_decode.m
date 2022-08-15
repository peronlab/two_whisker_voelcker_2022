%% DECODING figure
% 1) measure d' in dff space for weak v. strong touches across neural types, depths
% 3) look at population decoding of force - pool dff across neurons of specific type and estimate d'
% 2) look at decoding of: whisker identity (w1 v w2 v mw v sw), multi v. single (mw, sw, and direction discrim; by cell type 
% 3) look at population decoding of binary variables above

function figure_decode
    %% --- get the data   
    dat = get_two_whisker_data(0);
    force_redo = 0;
    dec_dat = get_two_whisker_decode_data (force_redo);
    [subvol_idx l_tag] = get_best_subvol_list (dat);
    ens_ids = get_ensemble_ids;

    ens_touch_color = [1 1 1]*0.75;
    ens_non_touch_color = [1 1 1]*0.75;

    %% --- plot setup
    fsize = 10;
    symbols = {'^','o','v','d','s','p','h'};
    offs = linspace(-0.2,0.2,length(dat.anim_data));

    fh_main = figure ('Position',[0 0 1500 750]);
    y_1 = 0.75;
    example_dff_v_dkappa_ax = axes(fh_main,'Position', [.05 y_1 .1 .2]);
    example_histogram_ax(1) = axes(fh_main, 'Position', [.2 y_1 0.1 0.2]);
    example_histogram_ax(2) = axes(fh_main, 'Position', [.35 y_1 0.1 0.2]);
    single_mouse_dprime_ax = axes(fh_main,'Position', [.5 y_1 .2 .2]);
    stats_dprime_ax = axes(fh_main,'Position', [.75 y_1 .2 .2]);

    y_2 = y_1 - 0.3;
    stats_auc_pooled_ax = axes(fh_main,'Position', [.75 y_2 .2 .2]);
    roc_example_ax = axes(fh_main,'Position', [.55 y_2 .1 .2]);
    stats_auc_w1mw_ax = axes(fh_main,'Position', [.05 y_2 .2 .2]);
    stats_auc_w2mw_ax = axes(fh_main,'Position', [.3 y_2 .2 .2]);

    y_3 = y_2 - 0.3;
    stats_auc_w1_w2_ax = axes(fh_main,'Position', [.05 y_3 .2 .2]);
    stats_auc_w1pr_ax = axes(fh_main,'Position', [.3 y_3 .2 .2]);
    stats_auc_w2pr_ax = axes(fh_main,'Position', [.55 y_3 .2 .2]);

    %% 2) d' for 2 unitary touch types for MW, SWb, SWu cells ; a) beeswarm plots single mouse by type/layer b) cross-animal by type/layer
   
    % a1) beeswarm for example mouse, single touch
    ex_ani = find(strcmp('an280201', {dat.settings.animals.name}));
    ex_field_idx = 1;
    field_idx = [1 1 3 1 3 3 1]; %  most numerous touch type per animal 1=W1P 2=W1R 3 4 same W2
    prim_scores = [];
    ids = [];

    for svi=1:length(dec_dat(ex_ani).subvol)
        prim_scores = [prim_scores dec_dat(ex_ani).subvol(svi).scores(:,field_idx(ex_ani))'];
        ids = [ids dec_dat(ex_ani).subvol(svi).ids];
    end
    aggregate_scores = prim_scores;
    aggregate_auc = normcdf(aggregate_scores/sqrt(2));

    axes(single_mouse_dprime_ax);
    for li=1:length(l_tag)
        mw_l_i = find(ismember(ids, intersect(dat.anim_data(ex_ani).types_by_id.mwCellIds, dat.anim_data(ex_ani).types_by_id.(l_tag{li}))));
        swb_l_i = find(ismember(ids, intersect(dat.anim_data(ex_ani).types_by_id.bidirCellIds, dat.anim_data(ex_ani).types_by_id.(l_tag{li}))));
        swu_l_i = find(ismember(ids, intersect(dat.anim_data(ex_ani).types_by_id.unidirCellIds, dat.anim_data(ex_ani).types_by_id.(l_tag{li}))));
        ens_t_l_i = find(ismember(ids, intersect(ens_ids(ex_ani).touch_ensemble_ids, dat.anim_data(ex_ani).types_by_id.(l_tag{li}))));
        ens_non_t_l_i = find(ismember(ids, intersect(ens_ids(ex_ani).non_touch_ensemble_ids, dat.anim_data(ex_ani).types_by_id.(l_tag{li}))));

        switch ex_field_idx
            case 1
                prim_vali =  find(ismember(ids, dat.anim_data(ex_ani).types_by_id.w1ProExclusiveTouchTrials));
            case 2
                prim_vali =  find(ismember(ids, dat.anim_data(ex_ani).types_by_id.w1RetExclusiveTouchTrials));
            case 3
                prim_vali =  find(ismember(ids, dat.anim_data(ex_ani).types_by_id.w2ProExclusiveTouchTrials));
            case 4
                prim_vali =  find(ismember(ids, dat.anim_data(ex_ani).types_by_id.w2RetExclusiveTouchTrials));
        end

        mw_l_i = intersect(mw_l_i,prim_vali);
        swb_l_i = intersect(swb_l_i,prim_vali);
        swu_l_i = intersect(swu_l_i,prim_vali);
        ens_t_l_i = intersect(ens_t_l_i,prim_vali);
        ens_non_t_l_i = intersect(ens_non_t_l_i,prim_vali);

        extern_plotSpread({aggregate_auc(mw_l_i)', aggregate_auc(swb_l_i)', aggregate_auc(swu_l_i)', aggregate_auc(ens_t_l_i)', aggregate_auc(ens_non_t_l_i)'},  ...
                         'distributionColors',{dat.settings.colors.mwColor, dat.settings.colors.bidiColor, dat.settings.colors.uniColor, ens_touch_color, ens_non_touch_color}, ...
                         'distributionMarkers','o', 'showMM', 3, 'xValues', [1 2 3 4 5] + 6*(li-1), 'markerSize', 3); % include median
    end
    set(single_mouse_dprime_ax, 'TickDir','out','FontSize',fsize);
    set(single_mouse_dprime_ax, 'YLabel', text('String','AUC_h_i_-_l_o'), 'XTick', [1 2 3 4 5 7 8 9 10 11 13 14 15 16 17],  ...
        'XTickLabel', {'L2 MW','SWb','SWu', 'TEns', 'NTEns','L3 MW','SWb','SWu','TEns', 'NTEns','L4 MW','SWb','SWu','TEns', 'NTEns'}, 'XTickLabelRotation', 45);
    title(single_mouse_dprime_ax, 'Example mouse w1 pro; must be sig w1p cell');
    aa = axis(single_mouse_dprime_ax);
    axis(single_mouse_dprime_ax,[aa(1) aa(2) 0.5 1]);

    % b1) data across all mice - median d' for each type/layer combo for most numerous touch of that mouse
    median_dprime = nan*zeros(length(dat.anim_data), length(l_tag), 6); % animal,layer,cell type
    thresh_dprime = 1;
    for ani=1:length(dat.anim_data)
        prim_scores = [];
        ids = [];
        for svi=1:length(dec_dat(ani).subvol)
            prim_scores = [prim_scores dec_dat(ani).subvol(svi).scores(:,field_idx(ani))'];
            ids = [ids dec_dat(ani).subvol(svi).ids];
        end
        prim_auc = normcdf(prim_scores/sqrt(2));

        switch field_idx(ani)
            case 1
                prim_vali =  find(ismember(ids, dat.anim_data(ani).types_by_id.w1ProExclusiveTouchTrials));
            case 2
                prim_vali =  find(ismember(ids, dat.anim_data(ani).types_by_id.w1RetExclusiveTouchTrials));
            case 3
                prim_vali =  find(ismember(ids, dat.anim_data(ani).types_by_id.w2ProExclusiveTouchTrials));
            case 4
                prim_vali =  find(ismember(ids, dat.anim_data(ani).types_by_id.w2RetExclusiveTouchTrials));
        end

        for li=1:length(l_tag)
            mw_l_i = find(ismember(ids, intersect(dat.anim_data(ani).types_by_id.mwCellIds, dat.anim_data(ani).types_by_id.(l_tag{li}))));
            swb_l_i = find(ismember(ids, intersect(dat.anim_data(ani).types_by_id.bidirCellIds, dat.anim_data(ani).types_by_id.(l_tag{li}))));
            swu_l_i = find(ismember(ids, intersect(dat.anim_data(ani).types_by_id.unidirCellIds, dat.anim_data(ani).types_by_id.(l_tag{li}))));
            ens_t_l_i = find(ismember(ids, intersect(ens_ids(ani).touch_ensemble_ids, dat.anim_data(ani).types_by_id.(l_tag{li}))));
            ens_non_t_l_i = find(ismember(ids, intersect(ens_ids(ani).non_touch_ensemble_ids, dat.anim_data(ani).types_by_id.(l_tag{li}))));

            all_ids = unique([dat.anim_data(ani).types_by_id.mwCellIds dat.anim_data(ani).types_by_id.bidirCellIds dat.anim_data(ani).types_by_id.unidirCellIds]);
            all_l_i = find(ismember(ids, intersect(all_ids, dat.anim_data(ani).types_by_id.(l_tag{li}))));

            mw_l_i = intersect(mw_l_i,prim_vali);
            swb_l_i = intersect(swb_l_i,prim_vali);
            swu_l_i = intersect(swu_l_i,prim_vali);
            ens_t_l_i = intersect(ens_t_l_i,prim_vali);
            ens_non_t_l_i = intersect(ens_non_t_l_i,prim_vali);
            all_l_i = intersect(all_l_i,prim_vali);

            qq = 0.5;

            median_auc(ani, li, 1) = quantile(prim_auc(mw_l_i),qq);
            median_auc(ani, li, 2) = quantile(prim_auc(swb_l_i),qq);
            median_auc(ani, li, 3) = quantile(prim_auc(swu_l_i),qq);
            median_auc(ani, li, 4) = quantile(prim_auc(ens_t_l_i),qq);
            median_auc(ani, li, 5) = quantile(prim_auc(ens_non_t_l_i),qq);
            median_auc(ani, li, 6) = quantile(prim_auc(all_l_i),qq);
        end
    end

    % now plot
    cell_types = {'mw','bidir','unidir','ensT','ensNT', 'allT'};
    cell_type_colors = {dat.settings.colors.mwColor, dat.settings.colors.bidiColor, dat.settings.colors.uniColor, ens_touch_color, ens_non_touch_color, [1 1 1]*0.5};
    cell_type_str = {'MW','SWb','SWu','ensT','ensNT', 'allT'};

    cti = 1:4;
    bar_plot_layer_type(stats_dprime_ax, median_auc(:,:,cti), 0.5, cell_type_colors(cti), cell_type_str(cti), fsize, 'all, for mouse most numerous uniwh type; must be sig cell')
    aa=axis(stats_dprime_ax);
    axis(stats_dprime_ax, [aa(1) aa(2) 0.5 1]);
    set(stats_dprime_ax, 'YLabel', text('String','AUC_h_i_-_l_o'));

    % b2) data across all mice - sum of d' for every cell, taken across each type/layer combo
    median_auc_pooled = zeros(length(dat.anim_data), length(l_tag), 6); % animal,layer,cell type
    median_auc_insig_pooled = zeros(length(dat.anim_data), length(l_tag), 6); % animal,layer,cell type
    for ani=1:length(dat.anim_data)
        ids = [];
        for svi=1:length(dec_dat(ani).subvol)
            ids = [ids dec_dat(ani).subvol(svi).ids];
        end
        auc_pooled = 0*ids;
        auc_pooled_zero_insig = 0*ids;

        for li=1:length(l_tag)
            mw_l_io = find(ismember(ids, intersect(dat.anim_data(ani).types_by_id.mwCellIds, dat.anim_data(ani).types_by_id.(l_tag{li}))));
            swb_l_io = find(ismember(ids, intersect(dat.anim_data(ani).types_by_id.bidirCellIds, dat.anim_data(ani).types_by_id.(l_tag{li}))));
            swu_l_io = find(ismember(ids, intersect(dat.anim_data(ani).types_by_id.unidirCellIds, dat.anim_data(ani).types_by_id.(l_tag{li}))));
            ens_t_l_io = find(ismember(ids, intersect(ens_ids(ani).touch_ensemble_ids, dat.anim_data(ani).types_by_id.(l_tag{li}))));
            ens_non_t_l_io = find(ismember(ids, intersect(ens_ids(ani).non_touch_ensemble_ids, dat.anim_data(ani).types_by_id.(l_tag{li}))));

            all_ids = unique([dat.anim_data(ani).types_by_id.mwCellIds dat.anim_data(ani).types_by_id.bidirCellIds dat.anim_data(ani).types_by_id.unidirCellIds]);
            all_l_io = find(ismember(ids, intersect(all_ids, dat.anim_data(ani).types_by_id.(l_tag{li}))));

            single_touch_ids = union(union(dat.anim_data(ani).types_by_id.w2ProExclusiveTouchTrials, dat.anim_data(ani).types_by_id.w2RetExclusiveTouchTrials), ...
                                     union(dat.anim_data(ani).types_by_id.w1ProExclusiveTouchTrials, dat.anim_data(ani).types_by_id.w1RetExclusiveTouchTrials));
            auc_pooled(find(isnan(auc_pooled))) = 0;
            auc_pooled_zero_insig(find(isnan(auc_pooled_zero_insig))) = 0;

            for field_idx = 1:4
                prim_scores = [];
                for svi=1:length(dec_dat(ani).subvol)
                    prim_scores = [prim_scores dec_dat(ani).subvol(svi).scores(:,field_idx)'];
                end

                % convert to AUC
                prim_auc = normcdf(prim_scores/sqrt(2));
                prim_auc(find(prim_auc < 0.5)) = 1-prim_auc(find(prim_auc < 0.5));
                prim_auc = prim_auc - 0.5; % baseline to 0 

                prim_vali =  find(ismember(ids, single_touch_ids));
                mw_l_i = intersect(mw_l_io,prim_vali);
                swb_l_i = intersect(swb_l_io,prim_vali);
                swu_l_i = intersect(swu_l_io,prim_vali);
                ens_t_l_i = intersect(ens_t_l_io,prim_vali);
                ens_non_t_l_i = intersect(ens_non_t_l_io,prim_vali);
                all_l_i = intersect(all_l_io,prim_vali);

                auc_pooled(mw_l_i) = auc_pooled(mw_l_i) + prim_auc(mw_l_i);
                auc_pooled(swb_l_i) = auc_pooled(swb_l_i) + prim_auc(swb_l_i);
                auc_pooled(swu_l_i) = auc_pooled(swu_l_i) + prim_auc(swu_l_i);
                auc_pooled(ens_t_l_i) = auc_pooled(ens_t_l_i) + prim_auc(ens_t_l_i);
                auc_pooled(ens_non_t_l_i) = auc_pooled(ens_non_t_l_i) + prim_auc(ens_non_t_l_i);
                auc_pooled(all_l_i) = auc_pooled(all_l_i) + prim_auc(all_l_i);

                switch field_idx
                    case 1
                        prim_vali =  find(ismember(ids, dat.anim_data(ani).types_by_id.w1ProExclusiveTouchTrials));
                    case 2
                        prim_vali =  find(ismember(ids, dat.anim_data(ani).types_by_id.w1RetExclusiveTouchTrials));
                    case 3
                        prim_vali =  find(ismember(ids, dat.anim_data(ani).types_by_id.w2ProExclusiveTouchTrials));
                    case 4
                        prim_vali =  find(ismember(ids, dat.anim_data(ani).types_by_id.w2RetExclusiveTouchTrials));
                end

                mw_l_i = intersect(mw_l_io,prim_vali);
                swb_l_i = intersect(swb_l_io,prim_vali);
                swu_l_i = intersect(swu_l_io,prim_vali);
                ens_t_l_i = intersect(ens_t_l_io,prim_vali);
                ens_non_t_l_i = intersect(ens_non_t_l_io,prim_vali);
                all_l_i = intersect(all_l_i,prim_vali);

                auc_pooled_zero_insig(mw_l_i) = auc_pooled_zero_insig(mw_l_i) + prim_auc(mw_l_i);
                auc_pooled_zero_insig(swb_l_i) = auc_pooled_zero_insig(swb_l_i) + prim_auc(swb_l_i);
                auc_pooled_zero_insig(swu_l_i) = auc_pooled_zero_insig(swu_l_i) + prim_auc(swu_l_i);
                auc_pooled_zero_insig(ens_t_l_i) = auc_pooled_zero_insig(ens_t_l_i) + prim_auc(ens_t_l_i);
                auc_pooled_zero_insig(ens_non_t_l_i) = auc_pooled_zero_insig(ens_non_t_l_i) + prim_auc(ens_non_t_l_i);
                auc_pooled_zero_insig(all_l_i) = auc_pooled_zero_insig(all_l_i) + prim_auc(all_l_i);
            end

            auc_pooled(find(auc_pooled == 0)) = nan;
            auc_pooled_zero_insig(find(auc_pooled_zero_insig == 0)) = nan;
            
            median_auc_pooled(ani, li, 1) = nanmedian(auc_pooled(mw_l_io));
            median_auc_pooled(ani, li, 2) = nanmedian(auc_pooled(swb_l_io));
            median_auc_pooled(ani, li, 3) = nanmedian(auc_pooled(swu_l_io)); 
            median_auc_pooled(ani, li, 4) = nanmedian(auc_pooled(ens_t_l_io)); 
            median_auc_pooled(ani, li, 5) = nanmedian(auc_pooled(ens_non_t_l_io)); 
            median_auc_pooled(ani, li, 6) = nanmedian(auc_pooled(all_l_io)); 

            median_auc_insig_pooled(ani, li, 1) = nanmedian(auc_pooled_zero_insig(mw_l_io));
            median_auc_insig_pooled(ani, li, 2) = nanmedian(auc_pooled_zero_insig(swb_l_io));
            median_auc_insig_pooled(ani, li, 3) = nanmedian(auc_pooled_zero_insig(swu_l_io)); 
            median_auc_insig_pooled(ani, li, 4) = nanmedian(auc_pooled_zero_insig(ens_t_l_io)); 
            median_auc_insig_pooled(ani, li, 5) = nanmedian(auc_pooled_zero_insig(ens_non_t_l_io)); 
            median_auc_insig_pooled(ani, li, 6) = nanmedian(auc_pooled_zero_insig(all_l_io)); 
        end
    end

    mdp = 0.5+(median_auc_insig_pooled/4);
    for l=1:3
        [h pval_tens_mw] = ttest(squeeze(mdp(:,l,4)), squeeze(mdp(:,l,1)));
        disp(sprintf('For layer %s, MW vs. t ens discrim hilo p-val: %0.3f ', l_tag{l}, pval_tens_mw));
    end

    % now plot
    %cti = 1:3;
    bar_plot_layer_type(stats_auc_pooled_ax, 0.5+(median_auc_insig_pooled(:,:,cti)/4), 0, cell_type_colors(cti), cell_type_str(cti), fsize, 'all, dprime avg across all 4 primitives')
    aa=axis(stats_auc_pooled_ax);
    axis(stats_auc_pooled_ax, [aa(1) aa(2) 0.5 1]);
    set(stats_auc_pooled_ax, 'YLabel', text('String','Mean AUC_h_i_-_l_o'));

    %% 3) decode for whisker identity, multi v. single
    cell_types = {'mw','bidir','unidir','w1Unidir','w2Unidir','ensT','ensNT','allT'};
    cell_type_colors = {dat.settings.colors.mwColor, dat.settings.colors.bidiColor, dat.settings.colors.uniColor, dat.settings.colors.w1Color, dat.settings.colors.w2Color, ens_touch_color, ens_non_touch_color, [1 1 1]*0.5};
    cell_type_str = {'MW','SWb','SWu','SWuW1','SWuW2','ensT','ensNT','allT'};
 
    cti = [1:3 6];
    median_auc_wh_identity = get_median_auc ('W1vW2DiscrimAUC', dat, dec_dat, l_tag, ens_ids)
    bar_plot_layer_type(stats_auc_w1_w2_ax, median_auc_wh_identity(:,:,cti), 0.5, cell_type_colors(cti), cell_type_str(cti), fsize, 'W1 v W2');
    aa=axis(stats_auc_w1_w2_ax);
    axis(stats_auc_pooled_ax, [aa(1) aa(2) 0.5 0.8]);

    for l=1:3
        [h pval_tens_bw] = ttest(squeeze(mdp(:,l,4)), squeeze(mdp(:,l,2)));
        disp(sprintf('For layer %s, bSW vs. t ens discrim t type p-val: %0.3f ', l_tag{l}, pval_tens_bw));
    end

    cti = [1 2 4 5];
    median_auc_wh1pr = get_median_auc ('W1ProW1RetDiscrimAUC', dat, dec_dat, l_tag, ens_ids);
    bar_plot_layer_type(stats_auc_w1pr_ax, median_auc_wh1pr(:,:,cti), 0.5, cell_type_colors(cti), cell_type_str(cti), fsize, 'W1p v W1r');
    aa=axis(stats_auc_w1_w2_ax);
    axis(stats_auc_pooled_ax, [aa(1) aa(2) 0.5 0.8]); 

    median_auc_wh2pr = get_median_auc ('W2ProW2RetDiscrimAUC', dat, dec_dat, l_tag, ens_ids);
    bar_plot_layer_type(stats_auc_w2pr_ax, median_auc_wh2pr(:,:,cti), 0.5, cell_type_colors(cti), cell_type_str(cti), fsize, 'W2p v W2r');
    aa=axis(stats_auc_w1_w2_ax);
    axis(stats_auc_pooled_ax, [aa(1) aa(2) 0.5 0.8]);

function bar_plot_layer_type_pop_data(ax, data_mat, y_base, type_colors, type_str, fsize, tstr)
    offs = linspace(-0.2,0.2,size(data_mat,1));
    symbols = {'^','o','v','d','s','p','h'};
    l_tag = {'L2','L3','L4'};

    pval_thresh = 0.05;
    
    x = 1;
    w = 0.8;
    ms = 5;
    hold(ax,'on');
    X = [];
    x_label = {};
    for li=1:size(data_mat,3)
        [pval_anov zz anov_stats] = anova1(squeeze(data_mat(:,:,li)), [], 'off');
        disp(['ANOVA pval for n=10 ' tstr ' for ' l_tag{li} ': ' num2str(pval_anov)]);
    
        for m=1:size(data_mat,2) % cell type
            vvec = squeeze(data_mat(:,m,li));
            Nsq = sqrt(length(find(~isnan(vvec))));
            h = nanmean(vvec);
            if (~isnan(h))
                rectangle(ax, 'Position',[x-(w/2) y_base w h-y_base], 'EdgeColor','none', 'FaceColor', type_colors{m});
                plot(ax, [1 1]*x, (h+(nanstd(vvec-0.5)*[1 -1]/Nsq)), 'k-','LineWidth',1);
                X(end+1) = x;
                x_label{end+1} = [l_tag{li} ' ' type_str{m}];
                for ani=1:size(data_mat,1)
                    plot(ax, x+offs(ani), data_mat(ani,m, li), symbols{ani},'MarkerFaceColor', [0 0 0], 'Color', [1 1 1]*.5, 'MarkerSize', ms); 
                end
            end
            x = x+1;
        end
        x = x+1;

        if (pval_anov < pval_thresh) 

            c = multcompare(anov_stats, 'CType', 'hsd','Display','off');
            for v=1:size(c,1);
                [h_tt pval_tt] = ttest(squeeze(data_mat(:,c(v,1),li)), squeeze(data_mat(:,c(v,2), li)));
                if(c(v,6) < pval_thresh)
                    disp(sprintf(' ** %s vs. %s p: %0.3f ; RAW paired ttest pvalue: %0.3f', type_str{c(v,1)}, type_str{c(v,2)}, c(v,6), pval_tt));
                else
                    disp(sprintf(' %s vs. %s p: %0.3f ; RAW paired ttest pvalue %0.3f', type_str{c(v,1)}, type_str{c(v,2)}, c(v,6), pval_tt));
                end
            end

            x_label{end} = ['**' x_label{end}]; 
        end        
    end

    for m=1:size(data_mat,3)
        [hyp pval_l24] = ttest(squeeze(data_mat(:,m,1)), squeeze(data_mat(:,m,3)));
        disp(sprintf(' %s L2 vs.L4 p: %0.3f ', type_str{m}, pval_l24));
        [hyp pval_l23] = ttest(squeeze(data_mat(:,m,1)), squeeze(data_mat(:,m,2)));
        disp(sprintf(' %s L2 vs.L3 p: %0.3f ', type_str{m}, pval_l23)); 
        [hyp pval_l34] = ttest(squeeze(data_mat(:,m,2)), squeeze(data_mat(:,m,3)));
        disp(sprintf(' %s L3 vs.L4 p: %0.3f ', type_str{m}, pval_l34));

    end

    set(ax, 'TickDir','out','FontSize',fsize);
    set(ax, 'YLabel', text('String','AUC (unitless)'), 'XTick', X, 'XTickLabel', x_label, 'XTickLabelRotation', 45);
    title(ax, tstr);
    aa = axis(ax);
    axis(ax, [aa(1) aa(2) 0.5 1]);

function bar_plot_layer_type(ax, data_mat, y_base, type_colors, type_str,  fsize, tstr)
    offs = linspace(-0.2,0.2,size(data_mat,1));
    symbols = {'^','o','v','d','s','p','h'};
    l_tag = {'L2','L3','L4'};
    
    pval_thresh = 0.05;
     
    disp(' ');
    x = 1;
    w = 0.8;
    ms = 5;
    hold(ax,'on');
    X = [];
    x_label = {};    
    for li=1:size(data_mat,2)
        for m=1:size(data_mat,3)
            h = nanmean(squeeze(data_mat(:,li,m)));
            rectangle(ax, 'Position',[x-(w/2) y_base w h-y_base], 'EdgeColor','none', 'FaceColor', type_colors{m});
            se_offs = nanstd(squeeze(data_mat(:,li,m)))/sqrt(size(data_mat,1));
            plot(ax, [x x], [h-se_offs h+se_offs], 'k-', 'LineWidth', 2);

            X(end+1) = x;
            x_label{end+1} = [l_tag{li} ' ' type_str{m}];
            
            for ani=1:size(data_mat,1)
               plot(ax, x+offs(ani), data_mat(ani,li,m), symbols{ani},'MarkerFaceColor', [0 0 0], 'Color', [1 1 1]*.5, 'MarkerSize', ms); 
            end
            x = x+1;
        end
        x = x+1;

        [pval_kw zz kw_stats] = kruskalwallis(squeeze(data_mat(:,li,:)), [], 'off');
        [pval_anov zz anov_stats] = anova1(squeeze(data_mat(:,li,:)), [], 'off');
        disp(['KW pval for ' tstr ' for ' l_tag{li} ': ' num2str(pval_kw)  ' ANOVA: ' num2str(pval_anov)]);
        if (pval_anov < pval_thresh) 

            % multiple comparison test.. (Dunn's test, for some reason Dunn-Sidak in MATLAB with KW works same as Dunn test)
            c = multcompare(anov_stats, 'CType', 'hsd','Display','off');
%            c = multcompare(kw_stats, 'CType', 'dunn-sidak','Display','off');
            for v=1:size(c,1);
                pval_sr = signrank(data_mat(:,li,c(v,1)), data_mat(:,li,c(v,2)));
                [h_tt pval_tt] = ttest(squeeze(data_mat(:,li,c(v,1))), squeeze(data_mat(:,li,c(v,2))));
                if(c(v,6) < pval_thresh)
                    disp(sprintf(' ** %s vs. %s p: %0.3f ; RAW paired ttest pvalue: %0.3f', type_str{c(v,1)}, type_str{c(v,2)}, c(v,6), pval_tt));
                else
                    disp(sprintf(' %s vs. %s p: %0.3f ; RAW paired ttest pvalue %0.3f', type_str{c(v,1)}, type_str{c(v,2)}, c(v,6), pval_tt));
                end
            end

            x_label{end} = ['**' x_label{end}]; 
        end
    end

    for m=1:size(data_mat,3)
        [hyp pval_l24] = ttest(squeeze(data_mat(:,1,m)), squeeze(data_mat(:,3,m)));
        disp(sprintf(' %s L2 vs.L4 p: %0.3f ', type_str{m}, pval_l24));
        [hyp pval_l23] = ttest(squeeze(data_mat(:,1,m)), squeeze(data_mat(:,2,m)));
        disp(sprintf(' %s L2 vs.L3 p: %0.3f ', type_str{m}, pval_l23)); 
        [hyp pval_l34] = ttest(squeeze(data_mat(:,2,m)), squeeze(data_mat(:,3,m)));
        disp(sprintf(' %s L3 vs.L4 p: %0.3f ', type_str{m}, pval_l34));
    end
    
    set(ax, 'TickDir','out','FontSize',fsize);
    set(ax, 'YLabel', text('String','AUC (unitless)'), 'XTick', X, 'XTickLabel', x_label, 'XTickLabelRotation', 45);
    title(ax, tstr)
    aa = axis(ax);
    axis(ax, [aa(1) aa(2) 0.5 aa(4)]);

function median_auc = get_median_auc (field_used, dat, dec_dat, l_tag, ens_ids)
    median_auc = nan*zeros(length(dat.anim_data), length(l_tag), 8); % animal,layer,cell type
    for ani=1:length(dat.anim_data)
        auc = [];
        ids = [];
        for svi=1:length(dec_dat(ani).subvol)
            if (isempty(dec_dat(ani).subvol(svi).(field_used)))
                auc = [auc nan*dec_dat(ani).subvol(svi).ids];
            else
                auc = [auc dec_dat(ani).subvol(svi).(field_used)];
            end
            ids = [ids dec_dat(ani).subvol(svi).ids];
        end

        % apply 'inversion'
        auc(find(auc < 0.5)) = 1-auc(find(auc < 0.5));

        for li=1:length(l_tag)
            mw_l_i = find(ismember(ids, intersect(dat.anim_data(ani).types_by_id.mwCellIds, dat.anim_data(ani).types_by_id.(l_tag{li}))));
            swb_l_i = find(ismember(ids, intersect(dat.anim_data(ani).types_by_id.bidirCellIds, dat.anim_data(ani).types_by_id.(l_tag{li}))));
            swu_l_i = find(ismember(ids, intersect(dat.anim_data(ani).types_by_id.unidirCellIds, dat.anim_data(ani).types_by_id.(l_tag{li}))));
            uw1_l_i = find(ismember(ids, intersect(dat.anim_data(ani).types_by_id.w1UnidirCellIds, dat.anim_data(ani).types_by_id.(l_tag{li}))));
            uw2_l_i = find(ismember(ids, intersect(dat.anim_data(ani).types_by_id.w2UnidirCellIds, dat.anim_data(ani).types_by_id.(l_tag{li}))));
            ens_t_l_i = find(ismember(ids, intersect(ens_ids(ani).touch_ensemble_ids, dat.anim_data(ani).types_by_id.(l_tag{li}))));
            ens_non_t_l_i = find(ismember(ids, intersect(ens_ids(ani).non_touch_ensemble_ids, dat.anim_data(ani).types_by_id.(l_tag{li}))));
            
            all_ids = unique([dat.anim_data(ani).types_by_id.mwCellIds dat.anim_data(ani).types_by_id.bidirCellIds dat.anim_data(ani).types_by_id.unidirCellIds]);
            all_l_i = find(ismember(ids, intersect(all_ids, dat.anim_data(ani).types_by_id.(l_tag{li}))));

            median_auc(ani, li, 1) = nanmedian(auc(mw_l_i));
            median_auc(ani, li, 2) = nanmedian(auc(swb_l_i));
            median_auc(ani, li, 3) = nanmedian(auc(swu_l_i));
            median_auc(ani, li, 4) = nanmedian(auc(uw1_l_i));
            median_auc(ani, li, 5) = nanmedian(auc(uw2_l_i));
            median_auc(ani, li, 6) = nanmedian(auc(ens_t_l_i));
            median_auc(ani, li, 7) = nanmedian(auc(ens_non_t_l_i));
            median_auc(ani, li, 8) = nanmedian(auc(all_l_i));
        end
    end


function [DVa DVb] = plot_example_cell(example_cell, dff_v_dkappa_ax, hist_ax, dat, dec_dat, fsize)
% example_cell: { anim_str, cell_id, whisker to show}
% *_ax: axes variables
% dat, dec_dat: data
% fsize: font size
    ani = find(strcmp(example_cell{1}, {dat.settings.animals.name}));
    ci = find(ismember(dat.anim_data(ani).ids, example_cell{2}));
    svi = dat.anim_data(ani).subvol_idx(ci);
    ci = find(ismember(dec_dat(ani).subvol(svi).ids, example_cell{2}));

    dffVDKappa = get_dff_v_dkappa (dat, ani, svi, ci);

    if (strcmp(example_cell{3},'W1'))
        dffVDKappa = dffVDKappa(1:2);
        scores = dec_dat(ani).subvol(svi).scores(ci,1:2);
    else
        dffVDKappa = dffVDKappa(3:4);
        scores = dec_dat(ani).subvol(svi).scores(ci,3:4);
    end

    % DFF v dkappa
    ms = 8; % marker size
    borcol = [1 1 1] *0.9;
    hold(dff_v_dkappa_ax,'on');
    for di=1:2
        dKappaVec = get_log_adjusted_dkappas (dffVDKappa{di}(:,1));
        goodi = find(~isnan(dKappaVec));
        dKappaVec = dKappaVec(goodi);
        dffVec = dffVDKappa{di}(goodi,2);

        if (di == 1) % do the magicv
            [irr sorti] = sort(abs(dKappaVec), 'ascend');
            g1 = sorti(1:round(length(sorti)/3));
            g2 = sorti((1+(round(length(sorti)/3)):((round(2*length(sorti)/3)-1))));
            g3 = sorti(round(2*length(sorti)/3):end);

            plot(dff_v_dkappa_ax, dKappaVec(g1), dffVec(g1), 'o','Color', borcol, 'MarkerFaceColor', [1 0.75 .75], 'MarkerSize', ms);
            plot(dff_v_dkappa_ax, dKappaVec(g2), dffVec(g2), 'o','Color', borcol,'MarkerFaceColor', [1 0.5 .5], 'MarkerSize', ms);
            plot(dff_v_dkappa_ax, dKappaVec(g3), dffVec(g3), 'o','Color', borcol,'MarkerFaceColor', [1 0.1 0.1], 'MarkerSize', ms);

            % Hi-lo histograms for PROTRACTION
            hist_range = 0:0.2:ceil(max(dffVec));
            if (ceil(max(dffVec)) == 1) ; hist_range = 0:0.05:1; end
            histogram(hist_ax(1), dffVec(g1), hist_range, 'EdgeColor','none','FaceColor',[1 0.75 .75]);
            histogram(hist_ax(2), dffVec(g3), hist_range, 'EdgeColor','none','FaceColor',[1 0.1 0.1]);

            DVa = dffVec(g1);
            DVb = dffVec(g3);
        else
            plot(dff_v_dkappa_ax, dKappaVec, dffVec, 'o','Color', borcol,'MarkerFaceColor', [1 1 1]*0, 'MarkerSize', ms);
        end
    end

    axi = axis(dff_v_dkappa_ax);
    axis( dff_v_dkappa_ax, [-5 5 -0.5 axi(4)]);
    fixAdjustedLogXLabels(dff_v_dkappa_ax);
    line('Parent', dff_v_dkappa_ax, 'XData', [0 0], 'YData', [-0.5 axi(4)], 'Color', [0 0 0], 'LineWidth', 1);
%    title(dff_v_dkappa_ax, sprintf('%s d-prime hi/lo pro: %0.2g ret: %0.2g', example_cell{3}, scores(1), scores(2)));
    auc = normcdf(scores/sqrt(2));
    title(dff_v_dkappa_ax, sprintf('%s AUC hi/lo pro: %0.2g ret: %0.2g', example_cell{3}, auc(1), auc(2)));
       
% log-adjusted axis plotting
function fixAdjustedLogXLabels(ax);
    minLogValue = -7; % in logarithmic space this is minimal
    zti = find(ax.XTick == 0);
    ax.XTickLabel{zti} = num2str(minLogValue);
    nti = find(ax.XTick < 0);
    for n=1:length(nti)
        ax.XTickLabel{nti(n)} = num2str(-1*ax.XTick(nti(n)) + minLogValue);
    end      
    pti = find(ax.XTick > 0);
    for p=1:length(pti)
        ax.XTickLabel{pti(p)} = num2str(ax.XTick(pti(p)) + minLogValue);
    end           

function dffVDKappa = get_dff_v_dkappa (dat, ani, svi, ci)
    oDir = pwd;
    global dec_example_sess;
    if (~isobject(dec_example_sess))

        datapath = [dat.settings.data_path_root filesep dat.settings.animals(ani).name filesep 'session_neuropilone_final_merged'];
        cd (datapath);

        fl = dir('an*sess.mat');
        load(fl(svi).name);
        dec_example_sess = s;
    end
    obj=dec_example_sess;

    [trials touchES] = obj.getWhiskerTouchTrialTypes;
    usedES = {touchES.w1ProExclusiveES , touchES.w1RetExclusiveES, touchES.w2ProExclusiveES, touchES.w2RetExclusiveES};

    timeWindow = [-2 5];
    excludeTimeWindow = [-1 1];
    timeUnit = pldo.timeSeries.second;
    allowOverlap = 0;

    w1Name = obj.whiskerTag{1};
    w2Name = obj.whiskerTag{2};

    dKappaTS{1} = obj.derivedDataTSA.getTimeSeriesById(20112);
    dKappaTS{2} = obj.derivedDataTSA.getTimeSeriesById(20113); 

    % prep returned variables
    scores = nan*ones(length(obj.caTSA),length(usedES));
    dffVDKappa = {};
    ids = obj.caTSA.ids;

    % grab cell
    ts = obj.caTSA.dffTimeSeriesArray.getTimeSeriesByIdx(ci);
    fitField = sprintf('dffBased_%sAbsMaxKappaZeroNotouch-%sAbsMaxKappaZeroNotouchtwoWhiskerModelAllDataFit', w1Name, w2Name);
    twmFit = obj.caTSA.cellFeatures.get(fitField);
    modelFit = twmFit(ci);

    % dff v force
    dffPredicted = ts.predictDffFromTwoWhiskerModelFit(obj.caTSA.trialIndices, trials, dKappaTS{1}, dKappaTS{2}, modelFit);

    useWhiskerForDffVDKappa = [1 1 2 2 ];
    proRetOnly = 'prpr'; % b = both; p = pro ; r = ret
    
    for e=1:length(usedES)
        [dataMat timeMat idxMat plotTimeVec ieIdxVec] = ts.getValuesAroundEvents(usedES{e}, timeWindow, timeUnit, allowOverlap, [], excludeTimeWindow);

        % pull kappa TS's
        vali = find(idxMat > 0);
        kappaDataMat = nan*dataMat;
        kappaDataMat(vali) = dKappaTS{useWhiskerForDffVDKappa(e)}.value(idxMat(vali));

        % pull model dff
        modelDataMat = nan*dataMat;
        modelDataMat(vali) = dffPredicted(idxMat(vali));

        % compute dff v. kappa - post touch
        timeWindowForTouchAndKappaAggregation = [0 2];
        tidxdff = find(plotTimeVec >= timeWindowForTouchAndKappaAggregation(1) & plotTimeVec < timeWindowForTouchAndKappaAggregation(2));

        dffVDKappaMat = nan*zeros(size(dataMat,1),2); % first value: dKappa ; second: dff
        modelDffVDKappaMat = nan*zeros(size(dataMat,1),2); % first value: dKappa ; second: dff
        for t=1:size(dataMat,1) % trial loop
            dffVDKappaMat(t,:) = [nanmean(kappaDataMat(t,tidxdff)) nanmean(dataMat(t,tidxdff)) ];
            modelDffVDKappaMat(t,:) = [nanmean(kappaDataMat(t,tidxdff)) nanmean(modelDataMat(t,tidxdff)) ];
        end

        % sign flip (for negatives) and threshold
        dKappaVec = dffVDKappaMat(:,1);
        negi = find(dKappaVec < 0);
        posi = find(dKappaVec >= 0);
        if(proRetOnly(e) == 'p') ; dKappaVec(posi) = nan; end
        if(proRetOnly(e) == 'r') ; dKappaVec(negi) = nan; end

        % store for plotting
        dffVDKappaMat(:,1) = dKappaVec;

        [irr sorti] = sort(abs(dffVDKappaMat(:,1)), 'ascend');
        sorti(find(isnan(dffVDKappaMat(sorti,1)))) = [];

        dffVDKappa{e} = dffVDKappaMat;
    end
    cd(oDir);
