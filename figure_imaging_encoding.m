function fh = figure_imaging_encoding
%  - (draw : cartoon of volumetric 2 column imaging)
%  - 3d stack -- as in roi gui 
%  - plane examples: one plane per layer, w/ segmentation
%  - cell count per layer histogram
%  - (TV: touch aligned responses for 3 force partitions for a few example cells)
%  - (schematic of encoding model)
%  - touch aligned responses for major types - w1, w2, 2 combo cells
%  - model RFs for these cells
%  - count of all w1p/w1r/w2p/w2r (16) possible combos - histo

    dat = get_two_whisker_data(0);

    %  --- setup figure
    fsize = 10;
    symbols = {'^','o','v','d','s','p','h'};
    
    fh = figure ('Position',[0 0 1600 500]);
    blank_cartoon_ax(1) = subplot('Position',[.05 .5 .15 .4]); % for volume cartoon
    for i=1:3
        example_planes_ax(i) = subplot('Position',[.45+(i-1)*.17 .5 .15 .4]);
    end
    blank_cartoon_ax(2) = subplot('Position',[.05 .05 .2 .4]); % for volume cartoon
    depth_count_ax = subplot('Position',[.275 .125 .15 .25]);
    type_frac_ax = subplot('Position',[.5 .175 .2 .3]);
    primitive_excess_ax(1) = subplot('Position',[.75 .175 .08 .3]);
    primitive_excess_ax(2) = subplot('Position',[.85 .175 .06 .3]);
    primitive_excess_ax(3) = subplot('Position',[.95 .175 .02 .3]);
 
    for i=1:length(blank_cartoon_ax)
        set(blank_cartoon_ax(i), 'XTick',[], 'YTick',[]);
    end

    % --- layer count

    % settings
    mfc = [1 1 1]*.5;
    msize = 8;
    h = .8;
    col = [0 0 0];

    % tally it up
    n_L2 = [];
    n_L3 = [];
    n_L4 = [];
    hold (depth_count_ax, 'on');
    for ani=1:length(dat.anim_data)
        n_L2(ani) = length(dat.anim_data(ani).types_by_id.L2);
        n_L3(ani) = length(dat.anim_data(ani).types_by_id.L3);
        n_L4(ani) = length(dat.anim_data(ani).types_by_id.L4);
    end

    % plotting itself
    Nsq = sqrt(sum(length(n_L2)));
    rectangle(depth_count_ax, 'Position',[0 2-h/2 nanmean(n_L2) h],  'EdgeColor','none', 'FaceColor', col);
    plot(depth_count_ax, nanmean(n_L2)+[1 -1]*nanstd(n_L2)/Nsq,[1 1]*2,  'k-','LineWidth',1); 
    rectangle(depth_count_ax, 'Position',[0 3-h/2 nanmean(n_L3) h],  'EdgeColor','none', 'FaceColor', col);
    plot(depth_count_ax, nanmean(n_L3)+[1 -1]*nanstd(n_L3)/Nsq,[1 1]*3,  'k-','LineWidth',1); 
    rectangle(depth_count_ax, 'Position',[0 4-h/2 nanmean(n_L4) h], 'EdgeColor','none', 'FaceColor', col);
    plot(depth_count_ax, nanmean(n_L4)+[1 -1]*nanstd(n_L4)/Nsq, [1 1]*4, 'k-','LineWidth',1); 
    
    % individual mice
    offs = linspace(-0.2,0.2, length(dat.anim_data));
    for ani=1:length(dat.anim_data)
        plot(depth_count_ax, n_L2(ani), 2+offs(ani), symbols{ani},'MarkerFaceColor', mfc, ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);
        plot(depth_count_ax, n_L3(ani), 3+offs(ani), symbols{ani},'MarkerFaceColor', mfc, ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);
        plot(depth_count_ax, n_L4(ani), 4+offs(ani), symbols{ani},'MarkerFaceColor', mfc, ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);
    end

    set(depth_count_ax,'TickDir','out','FontSize',fsize, 'YTick', [2 3 4], 'YTickLabel',{'L2','L3','L4'}, 'YDir','reverse');
    xlabel(depth_count_ax,'Cell count');    

    % --- type count

    % settubgs
    tstr = {'w1p','w1r','w2p','w2r'};
    score_field = { 'w1ProExclusiveTouchTrials',  'w1RetExclusiveTouchTrials', ...
                  'w2ProExclusiveTouchTrials', 'w2RetExclusiveTouchTrials'};
    base_color = {dat.settings.colors.w1Color, dat.settings.colors.w1Color, ...
                 dat.settings.colors.w2Color, dat.settings.colors.w2Color };
    types = {'1000','0100','1100','0010','0001','0011','1010','1001','0110','0101','1110','1101','1011','0111','1111'};
    primitives = [1 2 4 5]; % where is 1000 0100 0010 and 0001
    plot_ax = [0 0 1 0 0 1 1 1 1 1 2 2 2 2 3];
    type_x = [0 0 1 0 0 2 3 4 5 6 1 2 3 4 1];
    categories = [1 1 2 1 1 2 3 3 3 3 3 3 3 3 3]; % 1: uni single 2: bi single 3: multi

    % prebuild labels and colors
    tags = {};
    t_colors = {};
    for t=1:length(types)
        tags{t} = [];
        for f=1:4
            if (str2num(types{t}(f))) ; if(isempty(tags{t})) ; tags{t} = tstr{f} ; else ; tags{t} = [tags{t} '/' tstr{f}] ; end; end
        end
        if (strcmp(types{t}(3:4),'00'))
            t_colors{t} = dat.settings.colors.w1Color;
        elseif (strcmp(types{t}(1:2),'00'))
            t_colors{t} = dat.settings.colors.w2Color;
        else
            t_colors{t} = dat.settings.colors.mwColor;
        end

    end

    % main counting, and individual mouse plotting
    type_counts = zeros(length(dat.anim_data), length(types));
    hold(type_frac_ax,'on');

    for ani=1:length(dat.anim_data)
        vali_f = {};
        all = dat.anim_data(ani).ids;
        if (isfield(dat.anim_data(ani).types_by_id, 'excluded_always')) ; all = setdiff(dat.anim_data(ani).ids, dat.anim_data(ani).types_by_id.excluded_always) ; end
        cell_counts(ani) = length(all);
        all_vali = [];
        for f=1:length(score_field)
            vali_f{f} = dat.anim_data(ani).types_by_id.(score_field{f});
            all_vali = [all_vali vali_f{f}];
        end

        for t=1:length(types)
            vali = [];
            for f=1:4
                if (str2num(types{t}(f)) == 1)
                    if (isempty(vali))
                        vali = vali_f{f};
                    else
                        vali = intersect(vali, vali_f{f});
                    end
                end
            end  
            for f=1:4
                if (str2num(types{t}(f)) == 0)
                    vali = setdiff(vali, vali_f{f});
                end
            end
            type_counts(ani,t) = length(vali);
        end
    end
    frac_divider = repmat(cell_counts', 1,size(type_counts,2));
    type_frac = type_counts./frac_divider;
    mu_frac = nanmean(type_frac);

    % info for reporting
    frac_vec = nansum(type_frac(:, find(categories == 1))');
    disp(sprintf('Single-whisker unidirectional total percent +- SD: %0.2f +- %0.2f n=%d mice',  100*nanmean(frac_vec), 100*nanstd(frac_vec), length(dat.anim_data)));
    frac_vec = nansum(type_frac(:, find(categories == 2))');
    disp(sprintf('Single-whisker bidirectional total percent +- SD: %0.2f +- %0.2f n=%d mice',  100*nanmean(frac_vec), 100*nanstd(frac_vec), length(dat.anim_data)));
    frac_vec = nansum(type_frac(:, find(categories == 3))');
    disp(sprintf('Multi-whisker total percent +- SD: %0.2f +- %0.2f n=%d mice',  100*nanmean(frac_vec), 100*nanstd(frac_vec), length(dat.anim_data)));
    frac_vec = nansum(type_frac');
    disp(sprintf('Touch total percent +- SD: %0.2f +- %0.2f n=%d mice',  100*nanmean(frac_vec), 100*nanstd(frac_vec), length(dat.anim_data)));

    % finally plot
    w = .8;
    msize = 3;
    Nsq = sqrt(size(type_frac,1));
    for t=1:length(types)
        rectangle(type_frac_ax, 'Position',[t-w/2 0 w mu_frac(t)],  'EdgeColor','none', 'FaceColor', t_colors{t});
        vvec = type_frac(:,t);
        plot(type_frac_ax, [1 1]*t, nanmean(vvec)+[1 -1]*nanstd(vvec)/Nsq, 'k-','LineWidth',1); 
    end

     % single animal symbols
    offs = linspace(-0.2,0.2, length(dat.anim_data));
    for ani=1:length(dat.anim_data) % plot single mice
        plot(type_frac_ax,(1:length(types))+offs(ani), type_frac(ani,:), symbols{ani},'MarkerFaceColor', mfc, ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);
    end

    set(type_frac_ax,'TickDir','out','FontSize',fsize, 'XTick', 1:length(types), 'XTickLabel',tags, 'XTickLabelRotation', 45);
    ylabel (type_frac_ax,'Fraction of cells responsive');

    % --- primitives -- depends heavily on above stuff (Types)
    nonprimi = setdiff(1:length(types),primitives);
    actual_versus_predicted = [];
    for ti=1:length(nonprimi)
        t = nonprimi(ti);

        % figure out its primitives
        primi = find(types{t} == '1');

        frac_predicted = 1;
        for p=1:length(primi)
            frac_predicted = frac_predicted * mu_frac(primitives(primi(p)));
        end

        actual_versus_predicted(t) = mu_frac(t)/frac_predicted;

        % plot
        fcol = 'none';
        ecol = [0 0 0];
        ecol = t_colors{t};
        lw = 2;
        rectangle(primitive_excess_ax(plot_ax(t)), 'Position',[type_x(t)-w/2 0 w actual_versus_predicted(t)],  'EdgeColor', ecol, 'FaceColor', fcol, 'LineWidth', lw);
    end

    for i=1:3
        vali = find(plot_ax == i);
        set(primitive_excess_ax(i),'TickDir','out','FontSize',fsize, 'XTick', 1:length(vali), 'XTickLabel',tags(vali), 'XTickLabelRotation', 45);
        if (i == 1)
            ylabel (primitive_excess_ax(i),{'Ratio actual to predicted','assuming independence'});
        end
        aa = axis(primitive_excess_ax(i));
        axis(primitive_excess_ax(i), [0 length(vali)+1 0 aa(4)]);
    end


function plot_example_plane(filename, ax)        
    plot_rois = 0;

    load(filename);

    obj.fillToBorder;
    obj.colorByBaseScheme(4);
    im = obj.generateImage(0, 0, [], []);

    obj.workingImageSettings.chanColor{1} = [1 1 1];
    obj.workingImageSettings.pixelRange{1} = '[0 1]';
    obj.workingImage = 0*obj.workingImage + 0.85;
    obj.workingImage(1,1) = 0;
    obj.workingImage(end,end) = 1;
    im_with_rois = obj.generateImage(0, 1);

    % overlay two images
    im_combo = im;
    if (plot_rois)
        lm = repmat(obj.getLabMat, 1, 1, 3);
        roii = find(lm > 0);
        im_combo(roii) = im_with_rois(roii);
        im_combo(roii) = 0.5*im(roii)+0.5*im_with_rois(roii);
    end

    axes(ax);
    imshow(im_combo);
    set(gca,'XTick',[], 'YTick',[]);

function fh = plot_example_cell(obj, ts, dat, tstr)
    fh = figure('Position', [0 0 1200 600], 'Name', tstr);

    w1p_ax = subplot('Position',[.15 .8 .1 .15]);
    w1r_ax = subplot('Position',[.15 .55 .1 .15]);
    w2p_ax = subplot('Position',[.15 .3 .1 .15]);
    w2r_ax = subplot('Position',[.15 .05 .1 .15]);

    w1_rf_ax = subplot('Position',[.275 .55 .1 .2]);
    w2_rf_ax = subplot('Position',[.275 .05 .1 .2]);

    w2p1p_ax = subplot('Position',[.65 .8 .1 .15]);
    w2r1p_ax = subplot('Position',[.65 .55 .1 .15]);
    w1r2r_ax = subplot('Position',[.65 .3 .1 .15]);
    w1p2r_ax = subplot('Position',[.65 .05 .1 .15]);

    w2p1p_rf_real_ax = subplot('Position',[.775 .8 .1 .2]);
    w2r1p_rf_real_ax = subplot('Position',[.775 .55 .1 .2]);
    w1r2r_rf_real_ax = subplot('Position',[.775 .3 .1 .2]);
    w1p2r_rf_real_ax = subplot('Position',[.775 .05 .1 .2]);

    w2p1p_rf_model_ax = subplot('Position',[.9 .8 .1 .2]);
    w2r1p_rf_model_ax = subplot('Position',[.9 .55 .1 .2]);
    w1r2r_rf_model_ax = subplot('Position',[.9 .3 .1 .2]);
    w1p2r_rf_model_ax = subplot('Position',[.9 .05 .1 .2]);


    [trials touchES] = obj.getWhiskerTouchTrialTypes;
    scores = obj.getTwoWhiskerModelScores(ts.id);

    timeWindow = [-2 5];
    excludeTimeWindow = [-1 1];
    timeUnit = pldo.timeSeries.second;
    allowOverlap = 0;
    
    modelColor = [0 0 0];

    w1Name = obj.whiskerTag{1};
    w2Name = obj.whiskerTag{2};

    dKappaTS{1} = obj.derivedDataTSA.getTimeSeriesById(20112);
    dKappaTS{2} = obj.derivedDataTSA.getTimeSeriesById(20113); 

    ri = find(obj.caTSA.ids == ts.id);
    fitField = sprintf('dffBased_%sAbsMaxKappaZeroNotouch-%sAbsMaxKappaZeroNotouchtwoWhiskerModelAllDataFit', w1Name, w2Name);
    twmFit = obj.caTSA.cellFeatures.get(fitField);
    modelFit = twmFit(ri);

    % now let's plot dff-aligned for ALL touch types, looping thru types
    usedES = {touchES.w1ProExclusiveES, touchES.w1RetExclusiveES, touchES.w2ProExclusiveES, touchES.w2RetExclusiveES, ...
              touchES.w2pw1pES, touchES.w2rw1pES, touchES.w1pw2rES, touchES.w1rw2rES };
    ax = {w1p_ax, w1r_ax, w2p_ax, w2r_ax, w2p1p_ax, w2r1p_ax, w1r2r_ax, w1p2r_ax};
    tstr = {'w1p','w1r','w2p','w2r','w2pw1p','w2rw1p','w1pw2r','w1rw2r'};
    score_val = { scores.w1ProExclusiveTouchTrialsBalanced,  scores.w1RetExclusiveTouchTrialsBalanced, ...
                  scores.w2ProExclusiveTouchTrialsBalanced,  scores.w2RetExclusiveTouchTrialsBalanced, ...
                  scores.w2pw1pTrials,  scores.w2rw1pTrials,...
                  scores.w1pw2rTrials,  scores.w1rw2rTrials   };
    baseColor = {dat.settings.colors.w1Color, dat.settings.colors.w1Color, ...
                 dat.settings.colors.w2Color, dat.settings.colors.w2Color, ...
                 dat.settings.colors.w2w1Color, dat.settings.colors.w2w1Color, ...
                 dat.settings.colors.w1w2Color, dat.settings.colors.w1w2Color};
    scale_factor = {nan, nan, nan, nan, ...
                    modelFit.w2pw1p.kWaWb, modelFit.w2rw1p.kWaWb, modelFit.w1pw2r.kWaWb, modelFit.w1rw2r.kWaWb};

    for e=1:length(usedES)
        [dataMat timeMat idxMat plotTimeVec ieIdxVec] = ts.getValuesAroundEvents(usedES{e}, timeWindow, timeUnit, allowOverlap, [], excludeTimeWindow);

        for i=1:size(dataMat,1)
            line('Parent', ax{e}, 'XData', plotTimeVec, 'YData', dataMat(i,:), 'Color', min([1 1 1], baseColor{e}+0.5), 'LineWidth', 1);
        end
        line('Parent', ax{e}, 'XData', plotTimeVec, 'YData', nanmean(dataMat), 'Color', baseColor{e}, 'LineWidth', 3);
        axis(ax{e}, 'auto');
    end

    % cross axis harmonization for dff
    minVal = -0.5;
    axi = axis(ax{1});
    maxVal = axi(4);
    for i=1:length(ax);
        axiSize = axis(ax{i});
        maxVal = max(maxVal, axiSize(4));
    end

    for i=1:length(ax);
        axiSize = axis(ax{i});
        axis(ax{i}, [timeWindow(1) timeWindow(2) minVal maxVal]);
        line('Parent', ax{i}, 'XData', [0 0], 'YData', [minVal maxVal], 'Color', [0 0 0], 'LineWidth', 1);

        if (isnan(scale_factor{i}))
            set(ax{i}, 'TickDir','out', 'YLabel', text('String','\Delta F/F'), ...
                'Title', text('String', [tstr{i} ' exc R:' sprintf('%0.2f', score_val{i})]),'XLabel',text('String','Time to touch (s)'));
        else
            set(ax{i}, 'TickDir','out', 'YLabel', text('String','\Delta F/F'), ...
                'Title', text('String', [tstr{i} ' exc R:' sprintf('%0.2f', score_val{i}) ' sf: ' sprintf('%0.2f',scale_factor{i})]),'XLabel',text('String','Time to touch (s)'));
        end
     end

    % dff v force
    dffPredicted = ts.predictDffFromTwoWhiskerModelFit(obj.caTSA.trialIndices, trials, dKappaTS{1}, dKappaTS{2}, modelFit);

    threshold = [nan nan modelFit.w1.oPro modelFit.w1.oRet modelFit.w2.oPro modelFit.w2.oRet];
    useWhiskerForDffVDKappa = [1 2 1 1 2 2 ];
    proRetOnly = 'bbpprr'; % b = both; p = pro ; r = ret
    dffVDKappa = {};
    modelDffVDKappa = {};
    correctedThreshold = [];
 
    usedES = {touchES.w1ExclusiveTouchES, touchES.w2ExclusiveTouchES, touchES.w2pw1pES , touchES.w2rw1pES, touchES.w1pw2rES, touchES.w1rw2rES};
    baseColor = {dat.settings.colors.w1Color, dat.settings.colors.w2Color, dat.settings.colors.mwColor, dat.settings.colors.mwColor, ...
                 dat.settings.colors.mwColor, dat.settings.colors.mwColor};
    
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
        for t=1:size(dataMat,1)
            dffVDKappaMat(t,:) = [nanmean(kappaDataMat(t,tidxdff)) nanmean(dataMat(t,tidxdff)) ];
            modelDffVDKappaMat(t,:) = [nanmean(kappaDataMat(t,tidxdff)) nanmean(modelDataMat(t,tidxdff)) ];
        end

        % apply logarithmic transform, with sign flip (for negatives) and threshold
        dKappaVec = dffVDKappaMat(:,1);
        negi = find(dKappaVec < 0);
        posi = find(dKappaVec >= 0);
        if(proRetOnly(e) == 'p') ; dKappaVec(posi) = nan; end
        if(proRetOnly(e) == 'r') ; dKappaVec(negi) = nan; end

        dKappaVec = get_log_adjusted_dkappas (dKappaVec);

        if (isnan(threshold(e)) || threshold(e) == 0) 
            correctedThreshold(e) = nan;
        else
            if (proRetOnly(e) == 'p') 
                correctedThreshold(e) = get_log_adjusted_dkappas(-1*threshold(e));
            else
                correctedThreshold(e) = get_log_adjusted_dkappas(threshold(e));
            end
        end

        % store for plotting
        dffVDKappaMat(:,1) = dKappaVec;
        modelDffVDKappaMat(:,1) = dKappaVec;
        dffVDKappa{e} = dffVDKappaMat;
        modelDffVDKappa{e} = modelDffVDKappaMat;
    end

    % plot with real dFF uni whisker touch, then mw
    ax = {nan , nan, w2p1p_rf_real_ax, w2r1p_rf_real_ax, w1r2r_rf_real_ax, w1p2r_rf_real_ax};
    
    edge_color = [1 1 1]*0.75;
    msize = 40;
    lw = 1;

    maxAbsKappaValue = [0 0];
    additionalESPlot = [nan nan 1 1 2 2];
    for e=1:length(usedES)
        if (~ishandle(ax{e})); continue ; end
        hold(ax{e},'on');

        if (~isnan(additionalESPlot(e)))
            dKappaMat = dffVDKappa{additionalESPlot(e)};
            negi = find(dKappaMat(:,1) < 0);
            posi = find(dKappaMat(:,1) >= 0);
            if(proRetOnly(e) == 'r') ; dKappaMat(posi,:) = nan; end
            if(proRetOnly(e) == 'p') ; dKappaMat(negi,:) = nan; end 

            if (additionalESPlot(e) == 1) 
                circ_color = dat.settings.colors.w1Color;
            elseif (additionalESPlot(e) == 2) 
                circ_color = dat.settings.colors.w2Color;
            else
                circ_color = [.5 .5 .5];
            end

            sc = scatter(ax{e}, dKappaMat(:,1), dKappaMat(:,2), msize, 'MarkerEdgeColor',  circ_color, 'LineWidth', lw);
            alpha(sc,.5);
        end
        sc = scatter(ax{e}, dffVDKappa{e}(:,1), dffVDKappa{e}(:,2), msize, 'MarkerEdgeColor',  baseColor{e}, 'LineWidth', lw);
        axis(ax{e}, 'auto');
        maxAbsKappaValue(useWhiskerForDffVDKappa(e)) = max(maxAbsKappaValue(useWhiskerForDffVDKappa(e)),  max(abs(dffVDKappa{e}(:,1))));
    end

    % cross axis harmonization 
    for i=1:length(ax);
        if (~ishandle(ax{i})); continue ; end
        axiSize = axis(ax{i});
        axis(ax{i}, [-5 5 minVal maxVal]);
        line('Parent', ax{i}, 'XData', [0 0], 'YData', [minVal maxVal], 'Color', [0 0 0], 'LineWidth', 1);
        set(ax{i}, 'TickDir','out', 'YLabel', text('String','\Delta F/F'), 'XLabel', text('String','Adjusted log_1_0(\Delta\kappa)'));
    end
    fixAdjustedLogXLabels(ax(3:6));

    % plot with model dFF
    ax = {w1_rf_ax, w2_rf_ax, w2p1p_rf_model_ax, w2r1p_rf_model_ax, w1r2r_rf_model_ax, w1p2r_rf_model_ax};
    maxAbsKappaValue = [0 0 0 0 0 0];
    
    for e=1:length(usedES)
        hold(ax{e},'on');
        sc = scatter(ax{e}, dffVDKappa{e}(:,1), dffVDKappa{e}(:,2), msize, 'MarkerEdgeColor',  baseColor{e}, 'LineWidth', 1);
        sc = scatter(ax{e}, modelDffVDKappa{e}(:,1), modelDffVDKappa{e}(:,2), msize, 'MarkerEdgeColor', modelColor, 'LineWidth', 1);
        alpha(sc,.25);
        axis(ax{e}, 'auto');
        maxAbsKappaValue(useWhiskerForDffVDKappa(e)) = max(maxAbsKappaValue(useWhiskerForDffVDKappa(e)),  max(abs(dffVDKappa{e}(:,1))));
    end

    % cross axis harmonization 
    for i=1:length(ax);
        axiSize = axis(ax{i});
        axis(ax{i}, [-5 5 minVal maxVal]);
        line('Parent', ax{i}, 'XData', [0 0], 'YData', [minVal maxVal], 'Color', [0 0 0], 'LineWidth', 1);
        set(ax{i}, 'TickDir','out', 'YLabel', text('String','\Delta F/F'), 'XLabel', text('String','Adjusted log_1_0(\Delta\kappa)'));
    end
    fixAdjustedLogXLabels(ax);


% log-adjusted axis plotting
function fixAdjustedLogXLabels(ax);
    minLogValue = -7; % in logarithmic space this is minimal
    for i=1:length(ax)
        zti = find(ax{i}.XTick == 0);
        ax{i}.XTickLabel{zti} = num2str(minLogValue);
        nti = find(ax{i}.XTick < 0);
        for n=1:length(nti)
            ax{i}.XTickLabel{nti(n)} = num2str(-1*ax{i}.XTick(nti(n)) + minLogValue);
        end      
        pti = find(ax{i}.XTick > 0);
        for p=1:length(pti)
            ax{i}.XTickLabel{pti(p)} = num2str(ax{i}.XTick(pti(p)) + minLogValue);
        end           
    end

%% this is for supplement - plot w1 v w2 for inclusive and exclusive model - basically same
%% aggregate across all miceto make sure you aren't missing something
function junk
    animal_id = 'an280201';
    data_dir = [dat.settings.data_path_root filesep animal_id filesep 'session_neuropilone_final_merged'];
    cd (data_dir);
    fl = dir('an*sess.mat');
    for f=1:length(fl)
        load(fl(f).name);
    
        sce =  s.caTSA.cellFeatures.get('dffBased_c2AbsMaxKappaZeroNotouch-c3AbsMaxKappaZeroNotouchtwoWhiskerModelWithCrossvalScore');
        sci =  s.caTSA.cellFeatures.get('dffBased_inclusive_c2AbsMaxKappaZeroNotouch-c3AbsMaxKappaZeroNotouchtwoWhiskerModelWithCrossvalScore');
        w1i = [sci.w1ExclusiveTouchTrials];
        w2i = [sci.w2ExclusiveTouchTrials];
        w1e = [sce.w1ExclusiveTouchTrials];
        w2e = [sce.w2ExclusiveTouchTrials];

        figure;
        hold on; 
        plot(w1i, w2i, 'r.'); 
        plot(w1e, w2e, 'b.'); 
        plot ([-0.5 1],[-0.5 1],'k-');
        axis([-0.5 1 -0.5 1]); 
        title(strrep(fl(f).name, '_', '-'));

        for i=1:length(w1e) 
            plot([w1e(i) w1i(i)], [w2e(i) w2i(i)], 'k-'); 
        end
    end
