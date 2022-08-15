%% Figure for task stuff and kinematics; panel list:
%  - (draw : cartoon of task )
%  - trial-by-trial for one mouse
%  - average day performance trends, all mice
%  - whisking correlation to mean across days, all mice, non-touch trials - NOT REALLY
%  - (example : frame and kinematic parameters drawn over it)
%  - (example : trace of parameters for said frame, with touches colored) [280201 9/15 tr114 R P 9/6 tr290 ]
%  - frequency for all touch types - all 1 anim, by subvol by anim
function figure_task_kinematics
    dat = get_two_whisker_data(0);

    %  --- setup figure
    fsize = 10;
    fh = figure ('Position',[0 0 800 400]);
    frac_anim_ax = subplot('Position',[.1 .1 .4 .8]);
    frac_subvol_ax = subplot('Position',[.55 .1 .4 .8]);

    % --- frequency for all touch types
    touch_dat = get_touch_data(dat);
    touch_types = {'noTouchTrials', 'w1ExclusiveTouchTrials','w1ProExclusiveTouchTrials','w1RetExclusiveTouchTrials', ...
                   'w2ExclusiveTouchTrials','w2ProExclusiveTouchTrials','w2RetExclusiveTouchTrials', ...
                   'multiWhiskerTrials', 'w1pw2rTrials', 'w1rw2rTrials', 'w2pw1pTrials', 'w2rw1pTrials'};
    touch_color = {[0 0 0], [1 0.2 0.2], [1 .5 .5], [1 .5 .5], [0 0.5 1] , [0.5 .5 1], [.5 .5 1], ...
                   [1 0.2 1], [1 0.2 1], [1 0.2 1], [1 0.2 1], [1 0.2 1]};
    touch_x =[-1 1 2 3 5 6 7 9 10 11 12 13];                   
    symbols = {'^','o','v','d','s','p','h'};
    symbols = {'o','o','o','o','o','o','o'};

    animali = [touch_dat.anim_idx];
    subvoli = [touch_dat.plane_idx];
    touch_counts_animal = zeros(length(dat.settings.animals), length(touch_types));
    touch_counts_subvol = zeros(length(animali), length(touch_types));
    for a=1:length(dat.settings.animals)
        for t=1:length(touch_types)
            vali = find(animali == a);
            usv = unique(subvoli(vali));
            for s=1:length(usv)
                svi = find(subvoli == usv(s) & animali == a);
                if (strcmp(touch_types{t}, 'multiWhiskerTrials'))
                    for tt=t+1:length(touch_types)
                      touch_counts_animal(a,t) = touch_counts_animal(a,t) +  length(touch_dat(svi).trials.(touch_types{tt}));
                      touch_counts_subvol(svi, t) = touch_counts_subvol(svi,t) +  length(touch_dat(svi).trials.(touch_types{tt}));
                    end
                else
                    touch_counts_animal(a,t) = touch_counts_animal(a,t) +  length(touch_dat(svi).trials.(touch_types{t}));
                    touch_counts_subvol(svi, t) = touch_counts_subvol(svi,t) +  length(touch_dat(svi).trials.(touch_types{t}));
                end
            end
        end
    end

    % for each animal, aggregate bar plot and per-subvol plot
    Nsq = sqrt(length(dat.settings.animals));
    for a=1:length(dat.settings.animals)
        ai = find(animali == a);
        usv = unique(subvoli(ai));

        mean_counts = nanmean(touch_counts_subvol(ai,:));
        if (a == 4)
            axes(frac_subvol_ax);
        else
            figure('Position',[0 0 600 600]);
            subplot('Position', [.15 .15 .4 .8]);
        end
        hold on;
        half_width = 0.4;
        lab = {};
        msize = 8;
        for t=1:length(touch_types)
            lab{t} = touch_types{t}(1:6);
            rectangle('Position',[touch_x(t)-half_width 0 2*half_width mean_counts(t)],  'EdgeColor','none', 'FaceColor', touch_color{t});

            vvec = touch_counts_subvol(ai,t);
            plot(touch_x(t)*[1 1], nanmean(vvec)+[1 -1]*nanstd(vvec)/Nsq, 'k-','LineWidth', 1);
        end
        sv_offset = linspace(-0.3, 0.3, length(usv));
        for s=1:length(usv)
            svi = find(subvoli == usv(s) & animali == a);
            plot(touch_x+sv_offset(s), touch_counts_subvol(svi,:), symbols{s} , 'MarkerFaceColor', [0 0 0], ...
                    'Color',[1 1 1]*0.75, 'MarkerSize', msize);
        end    
        set(gca,'TickDir','out','FontSize',fsize, 'XTick', touch_x, 'XTickLabel',lab, 'XTickLabelRotation', 45);
        title([dat.settings.animals(a).name ' subvolumes']);
        ylabel('trials');
        
        if (a ~= 4)
            subplot('Position', [.6 .15 .4 .8]);
            hold on;
            msize = 8;
            for s=1:length(usv)
                svi = find(subvoli == usv(s) & animali == a);
                plot(touch_x, touch_counts_subvol(svi,:), [symbols{s} ':'] , 'MarkerFaceColor', [0 0 0], ...
                        'Color',[.5 .5 .5], 'MarkerSize', msize, 'LineWidth',1);
            end
            set(gca,'TickDir','out','FontSize',fsize, 'XTick', touch_x, 'XTickLabel',lab, 'XTickLabelRotation', 45);    
        end
    end   

    % overall agrregate plot
    mean_counts = nanmean(touch_counts_animal);
    axes(frac_anim_ax);
    hold on;
    half_width = 0.4;
    lab = {};
    msize = 8;
    Nsq = sqrt(length(dat.settings.animals));
    for t=1:length(touch_types)
        lab{t} = touch_types{t}(1:6);
        rectangle('Position',[touch_x(t)-half_width 0 2*half_width mean_counts(t)],  'EdgeColor','none', 'FaceColor', touch_color{t});

        vvec = touch_counts_animal(:,t);
        plot(touch_x(t)*[1 1], nanmean(vvec)+[1 -1]*nanstd(vvec)/Nsq, 'k-','LineWidth', 1);
        
    end
    ani_offset = linspace(-0.3, 0.3, length(dat.settings.animals));
    for a=1:length(dat.settings.animals)
        plot(touch_x+ani_offset(a), touch_counts_animal(a,:), symbols{a} , 'MarkerFaceColor', [0 0 0], ...
                'Color',[1 1 1]*0.75, 'MarkerSize', msize);
    end    
    set(gca,'TickDir','out','FontSize',fsize, 'XTick', touch_x, 'XTickLabel',lab, 'XTickLabelRotation', 45);
    title('All animals');
    ylabel('trials');

    if (0) % an extra one that is not so fun
        figure('Position',[0 0 600 600]);
        ax = subplot('Position', [.6 .15 .4 .8]);
        hold on;
        msize = 8;
        for a=1:length(dat.settings.animals)
            plot(touch_x, touch_counts_animal(a,:), [symbols{a} ':'], 'MarkerFaceColor', [0 0 0], ...
                'Color',[.5 .5 .5], 'MarkerSize', msize, 'LineWidth',1);
        end
        set(gca,'TickDir','out','FontSize',fsize, 'XTick', touch_x, 'XTickLabel',lab, 'XTickLabelRotation', 45);
    end

function touch_dat = get_touch_data(dat)
    dat_file = [dat.settings.summary_data_root filesep 'touch_dat.mat'];
    redo = 0; % every once in a while . . . 

    if (exist(dat_file,'file') & ~redo)
        load(dat_file);
    else % generate - THIS REQUIRES ACCESS TO RAW DATA; CONTACT PERON LAB
        touch_dat = [];
        t = 1;
        for a=1:length(dat.settings.animals)
            data_dir = [dat.settings.data_path_root filesep dat.settings.animals(a).name filesep 'session_neuropilone_final_merged'];
            cd(data_dir);
            fl = dir('an*sess.mat');
            for f=1:length(fl)
                load(fl(f).name);
                thisData = s.getWhiskerTouchTrialTypes;

                touch_dat(t).source_file = fl(f).name;
                touch_dat(t).trials = thisData;
                touch_dat(t).anim_idx = a;
                touch_dat(t).plane_idx = f; %oops -- should be like subvol_idx or sumtin oh well

                t = t+1;
            end
        end
        save(dat_file, 'touch_dat');
    end


