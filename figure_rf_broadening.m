%
% Panels:
%   a) 3d map for several mice showing single/dual-whisker encoding strength
%   b,c) bar plot showing density of sw/mw vs. depth across mice
%   d, e) mean, net score by layer/depth
%
function figure_rf_broadening
    single_call('dirboth') ; % can be dirw1 dirw2 as well
    single_call('w1w2');

function single_call(plot_mode);
    dat = get_two_whisker_data(0);

    % --- prelims
    switch plot_mode
        case 'w1w2'
            w1c = dat.settings.colors.w1Color;
            w2c = dat.settings.colors.w2Color;
            mwc = dat.settings.colors.mwColor;
            frac_max = [0.125 0.025];
            xtick_pos{1} = 0:.05:.125;
            xtick_pos{2} = 0:.01:.025;

        case {'dirboth','dirw1','dirw2'}
            w1c = dat.settings.colors.proColor;
            w2c = dat.settings.colors.retColor;
            mwc = dat.settings.colors.bidiColor;
            frac_max = [0.2 0.04];
            xtick_pos{1} = 0:.05:.2;
            xtick_pos{2} = 0:.01:.04;
    end


    % --- setup figure
    fsize = 10;
    fh = figure ('Position',[0 0 1600 1200]);
    example_3d_tw_ax = axes('Position',[.1 .5 .2 .4]);
    summary_tw_mw_ax = axes('Position',[.3 .1 .2 .15/2]);
    summary_tw_sw_ax = axes('Position',[.3 .35 .2 .45/2]);

    depth_v_frac = axes('Position',[.55 .1 .1 .6/2]);
    depth_v_mean_encoding = axes('Position',[.7 .1 .1 .6/2]);
    depth_v_norm_frac = axes('Position',[.85 .1 .1 .6/2]);


    % --- 3d plots TINA CODE GOT IT
    plot_3d_two_whiskers(dat, 'an280201', example_3d_tw_ax, [], [], .7, [2 20]);

    % --- count by depth plots

    % gather data
    for ani=1:length(dat.anim_data)
        [mw_vali all_vali w1_vali w2_vali sw_vali] = get_subgroups(plot_mode, dat, ani);

        sw1_vali = intersect(sw_vali, w1_vali);
        sw2_vali = intersect(sw_vali, w2_vali);

        L2_vali = find(ismember(dat.anim_data(ani).ids, dat.anim_data(ani).types_by_id.L2));
        L3_vali = find(ismember(dat.anim_data(ani).ids, dat.anim_data(ani).types_by_id.L3));
        L4_vali = find(ismember(dat.anim_data(ani).ids, dat.anim_data(ani).types_by_id.L4));

        n_L2(ani) = length(L2_vali);
        n_L3(ani) = length(L3_vali);
        n_L4(ani) = length(L4_vali);

        n_L2_T(ani) = length(intersect(L2_vali,all_vali));
        n_L3_T(ani) = length(intersect(L3_vali,all_vali));
        n_L4_T(ani) = length(intersect(L4_vali,all_vali));

        n_L2_T_mw(ani) = length(intersect(L2_vali,mw_vali));
        n_L3_T_mw(ani) = length(intersect(L3_vali,mw_vali));
        n_L4_T_mw(ani) = length(intersect(L4_vali,mw_vali));

        n_L2_T_sw(ani) = length(intersect(L2_vali,sw_vali));
        n_L3_T_sw(ani) = length(intersect(L3_vali,sw_vali));
        n_L4_T_sw(ani) = length(intersect(L4_vali,sw_vali));

        n_L2_T_sw1(ani) = length(intersect(L2_vali,sw1_vali));
        n_L3_T_sw1(ani) = length(intersect(L3_vali,sw1_vali));
        n_L4_T_sw1(ani) = length(intersect(L4_vali,sw1_vali));

        n_L2_T_sw2(ani) = length(intersect(L2_vali,sw2_vali));
        n_L3_T_sw2(ani) = length(intersect(L3_vali,sw2_vali));
        n_L4_T_sw2(ani) = length(intersect(L4_vali,sw2_vali));

    end

    disp(['Summary, n = ' num2str(ani)]);
    disp(sprintf('L2 T: %0.3fd SW: %0.3f MW: %0.3f frac T: %0.3f frac SW: %0.3f frac MW: %0.3f n_SW/n_T: %0.3f n_MW/n_T %0.3f', ...
         nanmean(n_L2_T), nanmean(n_L2_T_sw), nanmean(n_L2_T_mw), nanmean(n_L2_T./n_L2), nanmean(n_L2_T_sw./n_L2), nanmean(n_L2_T_mw./n_L2), nanmean(n_L2_T_sw./n_L2_T), nanmean(n_L2_T_mw./n_L2_T)));
    disp(sprintf('L2 T: %0.3fd SW: %0.3f MW: %0.3f frac T: %0.3f frac SW: %0.3f frac MW: %0.3f n_SW/n_T: %0.3f n_MW/n_T %0.3f', ...
         nanstd(n_L2_T), nanstd(n_L2_T_sw), nanstd(n_L2_T_mw), nanstd(n_L2_T./n_L2), nanstd(n_L2_T_sw./n_L2), nanstd(n_L2_T_mw./n_L2), nanstd(n_L2_T_sw./n_L2_T), nanstd(n_L2_T_mw./n_L2_T)));

    disp(sprintf('L3 T: %0.3f SW: %0.3f MW: %0.3f frac T: %0.3f frac SW: %0.3f frac MW: %0.3f n_SW/n_T: %0.3f n_MW/n_T %0.3f', ...
         nanmean(n_L3_T), nanmean(n_L3_T_sw), nanmean(n_L3_T_mw), nanmean(n_L3_T./n_L3), nanmean(n_L3_T_sw./n_L3), nanmean(n_L3_T_mw./n_L3), nanmean(n_L3_T_sw./n_L3_T), nanmean(n_L3_T_mw./n_L3_T)));
    disp(sprintf('L3 T: %0.3f SW: %0.3f MW: %0.3f frac T: %0.3f frac SW: %0.3f frac MW: %0.3f n_SW/n_T: %0.3f n_MW/n_T %0.3f', ...
         nanstd(n_L3_T), nanstd(n_L3_T_sw), nanstd(n_L3_T_mw), nanstd(n_L3_T./n_L3), nanstd(n_L3_T_sw./n_L3), nanstd(n_L3_T_mw./n_L3), nanstd(n_L3_T_sw./n_L3_T), nanstd(n_L3_T_mw./n_L3_T)));

    disp(sprintf('L4 T: %0.3f SW: %0.3f MW: %0.3f frac T: %0.3f frac SW: %0.3f frac MW: %0.3f n_SW/n_T: %0.3f n_MW/n_T %0.3f', ...
         nanmean(n_L4_T), nanmean(n_L4_T_sw), nanmean(n_L4_T_mw), nanmean(n_L4_T./n_L4), nanmean(n_L4_T_sw./n_L4), nanmean(n_L4_T_mw./n_L4), nanmean(n_L4_T_sw./n_L4_T), nanmean(n_L4_T_mw./n_L4_T)));
    disp(sprintf('L4 T: %0.3f SW: %0.3f MW: %0.3f frac T: %0.3f frac SW: %0.3f frac MW: %0.3f n_SW/n_T: %0.3f n_MW/n_T %0.3f', ...
         nanstd(n_L4_T), nanstd(n_L4_T_sw), nanstd(n_L4_T_mw), nanstd(n_L4_T./n_L4), nanstd(n_L4_T_sw./n_L4), nanstd(n_L4_T_mw./n_L4), nanstd(n_L4_T_sw./n_L4_T), nanstd(n_L4_T_mw./n_L4_T)));
 
    [h p_sw] = ttest((n_L4_T_sw./n_L4), (n_L3_T_sw./n_L3));
    [h p_mw] = ttest((n_L4_T_mw./n_L4), (n_L3_T_mw./n_L3));
    disp(sprintf(' L4 vs L3 frac sw: %0.3f mw: %0.3f', p_sw, p_mw));
 
    [h p_sw] = ttest((n_L2_T_sw./n_L2), (n_L3_T_sw./n_L3));
    [h p_mw] = ttest((n_L2_T_mw./n_L2), (n_L3_T_mw./n_L3));
    disp(sprintf(' L3 vs L2 frac sw: %0.3f mw: %0.3f', p_sw, p_mw));
    
    [h p_sm] = ttest((n_L2_T_sw1./n_L2), (n_L2_T_sw2./n_L2));
    disp(sprintf(' L2 sw1 v sw2: %0.3f mu1: %0.3f SD: %0.3f mu2: %0.3f SD: %0.3f', p_sm, nanmean(n_L2_T_sw1./n_L2), nanstd(n_L2_T_sw1./n_L2), nanmean(n_L2_T_sw2./n_L2), nanstd(n_L2_T_sw2./n_L2)));;
    [h p_sm] = ttest((n_L3_T_sw1./n_L3), (n_L3_T_sw2./n_L3));
    disp(sprintf(' L3 sw1 v sw2: %0.3f mu1: %0.3f SD: %0.3f mu2: %0.3f SD: %0.3f', p_sm, nanmean(n_L3_T_sw1./n_L3), nanstd(n_L3_T_sw1./n_L3), nanmean(n_L3_T_sw2./n_L3), nanstd(n_L3_T_sw2./n_L3)));;
    [h p_sm] = ttest((n_L4_T_sw1./n_L4), (n_L4_T_sw2./n_L4));
    disp(sprintf(' L4 sw1 v sw2: %0.3f mu1: %0.3f SD: %0.3f mu2: %0.3f SD: %0.3f', p_sm, nanmean(n_L4_T_sw1./n_L4), nanstd(n_L4_T_sw1./n_L4), nanmean(n_L4_T_sw2./n_L4), nanstd(n_L4_T_sw2./n_L4)));;

    % ploty!
    h = .8;

    offs = linspace(-1*h/4, h/4, length(dat.anim_data));
    msize = 8;
    w1cc = w1c/2;
    w2cc = w2c/2;
    mwcc = mwc/2;
    swc = [0 0 0];
    symbols = {'^','o','v','d','s','p','h'};

    sqN = sqrt(length(dat.anim_data));

    axes(summary_tw_sw_ax);
    hold on;
    vvec = (n_L2_T_sw1+n_L2_T_sw2)./n_L2;
    rectangle('Position',[0 10-h/2 nanmean((n_L2_T_sw1+n_L2_T_sw2)./n_L2) h],  'EdgeColor','none', 'FaceColor', swc);
    plot(nanmean(vvec)+[1 -1]*nanstd(vvec)/sqN, [10 10], 'k-', 'LineWidth', 1);
    vvec = (n_L3_T_sw1+n_L3_T_sw2)./n_L3;
    rectangle('Position',[0 6-h/2 nanmean((n_L3_T_sw1+n_L3_T_sw2)./n_L3) h],  'EdgeColor','none', 'FaceColor', swc);
    plot(nanmean(vvec)+[1 -1]*nanstd(vvec)/sqN, [6 6], 'k-', 'LineWidth', 1);
    vvec =(n_L4_T_sw1+n_L4_T_sw2)./n_L4;
    rectangle('Position',[0 2-h/2 nanmean((n_L4_T_sw1+n_L4_T_sw2)./n_L4) h],  'EdgeColor','none', 'FaceColor', swc); 
    plot(nanmean(vvec)+[1 -1]*nanstd(vvec)/sqN, [2 2], 'k-', 'LineWidth', 1);

    vvec = n_L2_T_sw1./n_L2;
    rectangle('Position',[0 8-h/2 nanmean(n_L2_T_sw1./n_L2) h],  'EdgeColor','none', 'FaceColor', w1c);
    plot(nanmean(vvec)+[1 -1]*nanstd(vvec)/sqN, [8 8], 'k-', 'LineWidth', 1);
    vvec = n_L3_T_sw1./n_L3;
    rectangle('Position',[0 4-h/2 nanmean(n_L3_T_sw1./n_L3) h],  'EdgeColor','none', 'FaceColor', w1c);
    plot(nanmean(vvec)+[1 -1]*nanstd(vvec)/sqN, [1 1]*4, 'k-', 'LineWidth', 1);
    vvec = n_L4_T_sw1./n_L4;
    rectangle('Position',[0 0-h/2 nanmean(n_L4_T_sw1./n_L4) h],  'EdgeColor','none', 'FaceColor', w1c); 
    plot(nanmean(vvec)+[1 -1]*nanstd(vvec)/sqN, [1 1]*0, 'k-', 'LineWidth', 1);

    vvec = n_L2_T_sw2./n_L2;
    rectangle('Position',[0 9-h/2 nanmean(n_L2_T_sw2./n_L2) h],  'EdgeColor','none', 'FaceColor', w2c);
    plot(nanmean(vvec)+[1 -1]*nanstd(vvec)/sqN, [1 1]*9, 'k-', 'LineWidth', 1);
    vvec = n_L3_T_sw2./n_L3;
    rectangle('Position',[0 5-h/2 nanmean(n_L3_T_sw2./n_L3) h],  'EdgeColor','none', 'FaceColor', w2c);
    plot(nanmean(vvec)+[1 -1]*nanstd(vvec)/sqN, [1 1]*5, 'k-', 'LineWidth', 1);
    vvec = n_L4_T_sw2./n_L4;
    rectangle('Position',[0 1-h/2 nanmean(n_L4_T_sw2./n_L4) h],  'EdgeColor','none', 'FaceColor', w2c);
    plot(nanmean(vvec)+[1 -1]*nanstd(vvec)/sqN, [1 1]*1, 'k-', 'LineWidth', 1);

    for a=1:length(dat.anim_data)
        plot ((n_L2_T_sw1(a)+n_L2_T_sw2(a))/n_L2(a), 10+offs(a), symbols{a}, 'MarkerFaceColor', swc, ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);
        plot ((n_L3_T_sw1(a)+n_L3_T_sw2(a))/n_L3(a), 6+offs(a), symbols{a}, 'MarkerFaceColor', swc, ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);
        plot ((n_L4_T_sw1(a)+n_L4_T_sw2(a))/n_L4(a), 2+offs(a), symbols{a}, 'MarkerFaceColor', swc, ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);        

        plot (n_L2_T_sw1(a)/n_L2(a), 8+offs(a), symbols{a}, 'MarkerFaceColor', w1cc, ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);
        plot (n_L3_T_sw1(a)/n_L3(a), 4+offs(a), symbols{a}, 'MarkerFaceColor', w1cc, ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);
        plot (n_L4_T_sw1(a)/n_L4(a), 0+offs(a), symbols{a}, 'MarkerFaceColor', w1cc, ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);        

        plot (n_L2_T_sw2(a)/n_L2(a), 9+offs(a), symbols{a}, 'MarkerFaceColor', w2cc, ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);
        plot (n_L3_T_sw2(a)/n_L3(a), 5+offs(a), symbols{a}, 'MarkerFaceColor', w2cc, ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);
        plot (n_L4_T_sw2(a)/n_L4(a), 1+offs(a), symbols{a}, 'MarkerFaceColor',w2cc, ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);        

    end

    axis(summary_tw_sw_ax, [0 frac_max(1) -0.5 11.5]);
    set(summary_tw_sw_ax,'TickDir','out','XTick', xtick_pos{1}, 'YTick',[1 5 9]-0.5, 'YTickLabel',{'L4','L3','L2'}, 'FontSize', fsize);

    axes(summary_tw_mw_ax);
    hold on;
    rectangle('Position',[0 2-h/2 nanmean(n_L2_T_mw./n_L2) h],  'EdgeColor','none', 'FaceColor', mwc);
    rectangle('Position',[0 1-h/2 nanmean(n_L3_T_mw./n_L3) h],  'EdgeColor','none', 'FaceColor', mwc);
    rectangle('Position',[0 0-h/2 nanmean(n_L4_T_mw./n_L4) h],  'EdgeColor','none', 'FaceColor', mwc); 

    for a=1:length(dat.anim_data)
        plot (n_L2_T_mw(a)/n_L2(a), 2+offs(a), symbols{a}, 'MarkerFaceColor',mwcc , ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);
        plot (n_L3_T_mw(a)/n_L3(a), 1+offs(a), symbols{a}, 'MarkerFaceColor',mwcc , ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);
        plot (n_L4_T_mw(a)/n_L4(a), 0+offs(a), symbols{a}, 'MarkerFaceColor',mwcc , ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);        
    end

    vvec = n_L2_T_mw./n_L2;
    plot(nanmean(vvec)+[1 -1]*nanstd(vvec)/sqN, [1 1]*2, 'k-', 'LineWidth', 1);
    vvec = n_L3_T_mw./n_L3;
    plot(nanmean(vvec)+[1 -1]*nanstd(vvec)/sqN, [1 1]*1, 'k-', 'LineWidth', 1);
    vvec = n_L4_T_mw./n_L4;
    plot(nanmean(vvec)+[1 -1]*nanstd(vvec)/sqN, [1 1]*0, 'k-', 'LineWidth', 1);

    

    axis(summary_tw_mw_ax, [0 frac_max(2) -0.5 2.5]);
    set(summary_tw_mw_ax,'TickDir','out','XTick', xtick_pos{2}, 'YTick',[0 1 2], 'YTickLabel',{'L4','L3','L2'}, 'FontSize', fsize);
    xlabel('Fraction of cells');

    % --- score by depth plots
    bin_size = 0.25; % in normalized units, which go from 2 to 5 and correspond to layer.
    min_cells_in_bin = 1; % at least this many cells or stay nan
    bins = 1:bin_size:5;
    bin_centers = bins(1:end-1)+(bin_size/2);
    
    w1cc = [1 0.5 0.5];
    w2cc = [.5 .5 1];
    mwcc = [1 .5 1];

    comboc = [0 0 0];

    mean_score_by_depth_mw = nan*zeros(length(dat.anim_data), length(bin_centers));
    frac_by_depth_mw = nan*zeros(length(dat.anim_data), length(bin_centers));

    mean_score_by_depth_w1 = nan*zeros(length(dat.anim_data), length(bin_centers));
    frac_by_depth_w1 = nan*zeros(length(dat.anim_data), length(bin_centers));

    mean_score_by_depth_w2 = nan*zeros(length(dat.anim_data), length(bin_centers));
    frac_by_depth_w2 = nan*zeros(length(dat.anim_data), length(bin_centers));

    hold(depth_v_mean_encoding, 'on');
    hold(depth_v_frac, 'on');
    hold(depth_v_norm_frac, 'on');
    lw = 1;   

    [prefTouchScores prefBalancedTouchScores] = get_preferred_scores(dat); % new thing - we use the MEAN of preferred scores for each cell (ones where it had valid response)
    for ani=1:length(dat.anim_data)
        ts = struct_nan_blanks(dat.anim_data(ani).touch_scores);
        [mw_vali all_vali w1_vali w2_vali sw_vali] = get_subgroups(plot_mode, dat, ani);

        sw1_vali = intersect(sw_vali, w1_vali);
        sw2_vali = intersect(sw_vali, w2_vali);

        % for fairest comparison use same score -- OVERALL
        w1_scores = [dat.anim_data(ani).touch_scores.allTrials];
        w2_scores = [dat.anim_data(ani).touch_scores.allTrials];
        mw_scores = [dat.anim_data(ani).touch_scores.allTrials];

        mw_scores = prefTouchScores{ani};
        w1_scores = prefTouchScores{ani};
        w2_scores = prefTouchScores{ani};

        for b=1:length(bin_centers)
            rel_depth_vali = find(dat.anim_data(ani).z_normalized > bins(b) & dat.anim_data(ani).z_normalized <= bins(b+1));

            if (length(intersect(rel_depth_vali,mw_vali)) >= min_cells_in_bin) 
                mean_score_by_depth_mw(ani,b) = nanmean(mw_scores(intersect(rel_depth_vali,mw_vali)));
                frac_by_depth_mw(ani,b) = length(intersect(rel_depth_vali,mw_vali))/length(rel_depth_vali);
            end
                
            if (length(intersect(rel_depth_vali,sw1_vali)) >= min_cells_in_bin) 
                mean_score_by_depth_w1(ani,b) = nanmean(w1_scores(intersect(rel_depth_vali,sw1_vali)));
                frac_by_depth_w1(ani,b) = length(intersect(rel_depth_vali,sw1_vali))/length(rel_depth_vali);
            end

            if (length(intersect(rel_depth_vali,sw2_vali)) >= min_cells_in_bin) 
                mean_score_by_depth_w2(ani,b) = nanmean(w2_scores(intersect(rel_depth_vali,sw2_vali)));
                frac_by_depth_w2(ani,b) = length(intersect(rel_depth_vali,sw2_vali))/length(rel_depth_vali);
            end
        end     
    end

    lw = 8;

    plot(depth_v_mean_encoding, nanmean(mean_score_by_depth_mw) , bin_centers, '-', 'LineWidth', lw, 'Color', mwc);
    plot(depth_v_mean_encoding, nanmean(mean_score_by_depth_w1) , bin_centers, '-', 'LineWidth', lw, 'Color', w1c);
    plot(depth_v_mean_encoding, nanmean(mean_score_by_depth_w2) , bin_centers, '-', 'LineWidth', lw, 'Color', w2c);
    plot(depth_v_mean_encoding, nanmean(0.5*mean_score_by_depth_w1 + 0.5*mean_score_by_depth_w2) , bin_centers, '-', 'LineWidth', lw, 'Color', comboc);
    set(depth_v_mean_encoding,'TickDir','out','YTick',[2 3 4], 'YTickLabel',{'L2','L3','L4'}, 'FontSize', fsize, 'YDir','reverse');
    xlabel(depth_v_mean_encoding, 'Mean encoding score');
    ylabel(depth_v_mean_encoding,'Normalized depth');
    aa = axis(depth_v_mean_encoding);
    for y = 2:4 ; plot(depth_v_mean_encoding,[0 1], [y y], 'k:', 'LineWidth', 5); end
    axis(depth_v_mean_encoding, [0 0.5 2 5]);

    plot(depth_v_frac, nanmean(frac_by_depth_mw) , bin_centers, '-', 'LineWidth', lw, 'Color', mwc);
    plot(depth_v_frac, nanmean(frac_by_depth_w1) , bin_centers, '-', 'LineWidth', lw, 'Color', w1c);
    plot(depth_v_frac, nanmean(frac_by_depth_w2) , bin_centers, '-', 'LineWidth', lw, 'Color', w2c);
    plot(depth_v_frac, nanmean(frac_by_depth_w1 + frac_by_depth_w2) , bin_centers, '-', 'LineWidth', lw, 'Color', comboc);
    set(depth_v_frac,'TickDir','out','YTick',[2 3 4], 'YTickLabel',{'L2','L3','L4'}, 'FontSize', fsize, 'YDir','reverse');
    xlabel(depth_v_frac, 'Fraction with type');
    ylabel(depth_v_frac,'Normalized depth');
    aa = axis(depth_v_frac);
    max_frac = 0.1;
    axis(depth_v_frac, [0 max_frac 2 5]);
    for y = 2:4 ; plot(depth_v_frac,[0 max_frac], [y y], 'k:', 'LineWidth', 5); end

    norm_frac_by_depth_mw = nanmean(frac_by_depth_mw);
    norm_frac_by_depth_w1 = nanmean(frac_by_depth_w1);
    norm_frac_by_depth_w2 = nanmean(frac_by_depth_w2);
    norm_frac_by_depth_mw = norm_frac_by_depth_mw/nanmax(norm_frac_by_depth_mw(3:end-2)); % bc of degeneracy
    norm_frac_by_depth_w1 = norm_frac_by_depth_w1/nanmax(norm_frac_by_depth_w1);
    norm_frac_by_depth_w2 = norm_frac_by_depth_w2/nanmax(norm_frac_by_depth_w2);

    norm_frac_by_depth_combo = nanmean(frac_by_depth_w1 + frac_by_depth_w2);
    norm_frac_by_depth_combo = norm_frac_by_depth_combo/nanmax(norm_frac_by_depth_combo);

    plot(depth_v_norm_frac, norm_frac_by_depth_mw , bin_centers, '-', 'LineWidth', lw, 'Color', mwc);
    plot(depth_v_norm_frac, norm_frac_by_depth_w1 , bin_centers, '-', 'LineWidth', lw, 'Color', w1c);
    plot(depth_v_norm_frac, norm_frac_by_depth_w2 , bin_centers, '-', 'LineWidth', lw, 'Color', w2c);
    plot(depth_v_norm_frac, norm_frac_by_depth_combo , bin_centers, '-', 'LineWidth', lw, 'Color', comboc);
    set(depth_v_norm_frac,'TickDir','out','YTick',[2 3 4], 'YTickLabel',{'L2','L3','L4'}, 'FontSize', fsize, 'YDir','reverse');
    xlabel(depth_v_norm_frac,'Peak-normalized fraction');
    ylabel(depth_v_norm_frac,'Normalized depth');
    aa = axis(depth_v_norm_frac);
    axis(depth_v_norm_frac, [0 1 2 5]);
    for y = 2:4 ; plot(depth_v_norm_frac,[0 1], [y y], 'k:', 'LineWidth', 5); end


