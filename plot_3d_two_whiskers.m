function plot_3d_two_whiskers(dat, animal, ax, mode, show_rejects, sat_val, size_range)
    if (nargin < 3 || isempty(ax)) ; figure; ax = axes ;end
    if (nargin < 4 || isempty(mode)) ; mode = 'two_whisker' ;end
    if (nargin < 5 || isempty(show_rejects)) ; show_rejects = 0 ; end
    if (nargin < 6 || isempty(sat_val)) ; sat_val = 0.75 ; end
    if (nargin < 7 || isempty(size_range)) ; size_range = [2 15]; end
    
    % settings
    w1_color = dat.settings.colors.w1Color;
    w2_color =  dat.settings.colors.w2Color;
    labels_shown = 1 ; 

    % clear axeds
    axes(ax);
    cla;
    hold on;

    % data
    ani = find(strcmp(animal, {dat.settings.animals.name}));
    if (strcmp(mode, 'two_whisker'))
        w1_score = [dat.anim_data(ani).touch_scores.w1ExclusiveTouchTrials];
        w2_score = [dat.anim_data(ani).touch_scores.w2ExclusiveTouchTrials];
        vali = find(ismember(dat.anim_data(ani).ids, union(dat.anim_data(ani).types_by_id.w1ExclusiveTouchTrials, dat.anim_data(ani).types_by_id.w2ExclusiveTouchTrials)));

    elseif (length(strfind(mode, 'directional')))
        w1_color = dat.settings.colors.proColor;
        w2_color = dat.settings.colors.retColor;
        sat_val = .75;    

        % trickier -- compute directionality among all cells above ; 
        %    for SW cells, directionality is just for the whisker that works ; for MW cells, use mean of both
        touch = struct_nan_blanks(dat.anim_data(ani).touch_scores);
        w1_score = [touch.w1ExclusiveTouchTrials];
        w2_score = [touch.w2ExclusiveTouchTrials];
        
        pw1_score = abs([touch.w1ProExclusiveTouchTrials]);
        pw2_score = abs([touch.w2ProExclusiveTouchTrials]);
        rw1_score = abs([touch.w1RetExclusiveTouchTrials]);
        rw2_score = abs([touch.w2RetExclusiveTouchTrials]);

        if (strcmp(mode, 'directional_w1_only'))
            pw2_score = 0*pw2_score;
            rw2_score = 0*rw2_score;
        elseif (strcmp(mode, 'directional_w2_only'))
            pw1_score = 0*pw1_score;
            rw1_score = 0*rw1_score;
        end

        vali = find(ismember(dat.anim_data(ani).ids, union(dat.anim_data(ani).types_by_id.w1ExclusiveTouchTrials, dat.anim_data(ani).types_by_id.w2ExclusiveTouchTrials)));
        w1_vali = find(ismember(dat.anim_data(ani).ids, dat.anim_data(ani).types_by_id.w1ExclusiveTouchTrials));
        w2_vali = find(ismember(dat.anim_data(ani).ids, dat.anim_data(ani).types_by_id.w2ExclusiveTouchTrials));
        mw_vali = find(ismember(dat.anim_data(ani).ids, intersect(dat.anim_data(ani).types_by_id.w1ExclusiveTouchTrials, dat.anim_data(ani).types_by_id.w2ExclusiveTouchTrials)));
        sw_w1_vali = setdiff(w1_vali, mw_vali);
        sw_w2_vali = setdiff(w2_vali, mw_vali);

        % DI for SW cells - note that w1 now is actually pro ; w2 is ret
        w1_score = nan * w1_score;
        w2_score = nan * w2_score;

        w1_score(mw_vali) = 0.5*pw1_score(mw_vali) + 0.5*pw2_score(mw_vali);
        w2_score(mw_vali) = 0.5*rw1_score(mw_vali) + 0.5*rw2_score(mw_vali);

        w1_score(sw_w1_vali) = pw1_score(sw_w1_vali);
        w2_score(sw_w1_vali) = rw1_score(sw_w1_vali);

        w1_score(sw_w2_vali) = pw2_score(sw_w2_vali);
        w2_score(sw_w2_vali) = rw2_score(sw_w2_vali);

        w1_score(find(w1_score < 0)) = nan; % set retraction scores to 0
        w2_score(find(w2_score < 0)) = nan; % set protraction scores to 0

        % who is plotted
        w1_score(find(isnan(w1_score))) = 0;
        w2_score(find(isnan(w2_score))) = 0;
        vali = union(find(w1_score > 0), find(w2_score > 0));
    end

    w1_score = w1_score/sat_val;
    w1_score(find(w1_score > 1)) = 1;
    w1_score(find(w1_score < 0)) = 0;
    w2_score = w2_score/sat_val;
    w2_score(find(w2_score > 1)) = 1;
    w2_score(find(w2_score < 0)) = 0;
    
    xyz = nan*zeros(length(w1_score),3);
    xyz(:,1:2) = dat.anim_data(ani).xy_micron'; 
    xyz(:,3) = dat.anim_data(ani).z_micron_corrected'; 

    % non-categorical
    if (show_rejects)
        invali = setdiff(1:length(xyz), vali);
        plot3(xyz(invali,1),xyz(invali,2),xyz(invali,3),'.','Color',  [1 1 1]*0.5, 'MarkerSize', 0.5);
    end

    disp(['Max depth: ' num2str(max(xyz(:,3)))]);

    cmap_w1 = [linspace(0,w1_color(1),100); linspace(0,w1_color(2),100); linspace(0,w1_color(3),100)]';
    cmap_w2 = [linspace(0,w2_color(1),100); linspace(0,w2_color(2),100); linspace(0,w2_color(3),100)]';

    % categorical
    for v=1:length(vali)
        i = vali(v);

        this_w1_color = cmap_w1(max(1,round(w1_score(i)*100)),:);
        this_w2_color = cmap_w2(max(1,round(w2_score(i)*100)),:);
          
        col = this_w1_color + this_w2_color;
        col = min(col,[1 1 1]);
        max_score = max(w1_score(i),w2_score(i));
        msize = round(max(size_range(1), (max_score^2)*size_range(2)));

        plot3(xyz(i,1),xyz(i,2),xyz(i,3),'o','Color',col,'MarkerFaceColor',col,'MarkerSize',msize);
    end
    
    ax_w = 750;
    ax_h = 500;

    % pedestal
    pedestal_color = [1 1 1]*0.85;
    fill3([0 ax_w ax_w 0 ], [0 0 ax_w ax_w], [ 1 1 1 1 ]*500, pedestal_color, 'EdgeColor', pedestal_color);

    % laminar bounds
    l12_border = min(dat.anim_data(ani).z_micron_corrected);
    plot3([0 ax_w ax_w], [0 0 ax_w ], [1 1 1 ]*l12_border, ':', 'Color', [0 0 0], 'LineWidth',10);
    text(0, 0, l12_border-30, 'L1','FontSize', 15);
    l23_border = dat.settings.animals(ani).l2l3_border;
    plot3([0 ax_w ax_w], [0 0 ax_w ], [1 1 1 ]*l23_border, ':', 'Color', [0 0 0], 'LineWidth',10)
    text(0, 0, l23_border-30, 'L2','FontSize', 15);
    l34_border = dat.settings.animals(ani).l3l4_border;
    plot3([0 ax_w ax_w], [0 0 ax_w ], [1 1 1 ]*l34_border, ':', 'Color', [0 0 0], 'LineWidth',10)
    text(0, 0, l34_border-30, 'L3','FontSize', 15);
    text(0, 0, l34_border+30, 'L4','FontSize', 15);

    % touchups
    set(gca,'XTick',[],'YTick',[],'YDir','reverse', 'ZDir','reverse','FontSize',15);
    %set(gca,'CameraPosition', [-3000 4000 -750]); % good for 1/2 as tall
    set(gca,'CameraPosition', [-3000 5000 0]);
    ax = gca;

    axis([0 ax_w 0 ax_w 0 ax_h ]);
    set(ax.YAxis,'visible','off');
    set(ax.XAxis,'visible','off');

    if (labels_shown)
        ax.ZTick= [100 200 300 400];
        ax.ZLabel.String = 'Depth \mum';
        ax.Title.String = sprintf('Two whisker 3d plot, n=%d cells, %s', length(w1_score), animal);
    else
        set(ax.ZAxis,'visible','off');
        ax.ZTick= [];
    end

