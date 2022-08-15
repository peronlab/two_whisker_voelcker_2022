%% this will return structures of heavily processed data ; normally load from file, but also contains create-de-novo stuff
function dat = get_two_whisker_data(params)
    if (nargin < 1) 
        force_redo = 0;
        params = [];
    elseif (~isstruct(params)) ; 
        force_redo = params ; 
    else 
        if (isfield(params, 'force_redo'))
            force_redo = params.force_redo;
        end
    end
    if (force_redo == 0) ; disp('get_two_whisker_data :: force_redo off.  If you have new data, be sure to call get_two_whisker_data(1)'); end

    exclude_matrix_name = 'excluded_always'; 
    if (isstruct(params) & isfield(params, 'exclude_matrix_name'))
        exclude_matrix_name = params.exclude_matrix_name;
    end

    settings = get_two_whisker_settings;
    summary_data_file = [settings.summary_data_root filesep 'two_whisker_main_summary_data.mat'];
    first_time = 0; % some things only should happen once

    % the core loop that regenrates summary file if requested/missing
    % --- THIS WILL NOT WORK WITHOUT RAW DATA ; RAW DATA IS MASSIVE (~500GB), CONTACT PERONLAB IF YOU WOULD LIKE IT
    if (force_redo | ~exist(summary_data_file, 'file'))
        if (force_redo) ; disp('force redo enabled ; doing it all de novo.  go make a cup of coffee.'); end

        disp([summary_data_file ' DOES NOT EXIST']);
        disp(' ');
        disp('WARNING: Could not find summary data. Your path is almost certainly not set properly.');
        disp('This code will only work if you have raw data; this is massive (~500GB) and will be');
        disp('made available upon request if you truly want it.  Contact PeronLab!');
        pause;
        
        % dir?
        if (~exist(settings.summary_data_root)) 
            disp([settings.summary_data_root ' does not exist; creating']);
            mkdir(settings.summary_data_root);
        end

        % animal loop for loading then assigning thresholds
        for a=1:length(settings.animals)
            data_dir = [settings.data_path_root filesep settings.animals(a).name filesep 'session_neuropilone_final_merged'];
            roi_dir = [settings.data_path_root filesep settings.animals(a).name filesep 'rois'];

            anim_data = get_animal_data(data_dir, roi_dir, settings.animals(a));
            anim_data = assign_types(anim_data);
            dat.anim_data(a) = anim_data;
        end

        % write to file
        dat.settings = settings;
        dat.generate_time_stamp = datestr(now);

        % break up
        inc_dat = [];
        inc_dat.settings = dat.settings;
        inc_dat.generate_time_stamp = dat.generate_time_stamp;

        for a=1:length(dat.settings.animals)
            inc_dat.anim_data(a).inc_touch_scores = dat.anim_data(a).inc_touch_scores;
            inc_dat.anim_data(a).inc_touch_scores_sh = dat.anim_data(a).inc_touch_scores_sh;
        end

        % saveZ
        for a=1:length(dat.anim_data)
            anim_data = dat.anim_data(a);
            save(strrep(summary_data_file, '.mat',sprintf ('_%02d.mat', a)), 'anim_data');
        end
        tmp_dat = dat;
        dat = rmfield(dat,'anim_data');
        save(summary_data_file, 'dat');
        dat = tmp_dat;
        clear tmp_dat;

        first_time = 1;
    else
        disp([summary_data_file ' exists; force redo if you want - call with parameter 1']);
        global dat;
        if (isempty (dat))
            load(summary_data_file);
            for a=1:length(dat.settings.animals)
                fname = strrep(summary_data_file, '.mat',sprintf ('_%02d.mat', a));
                disp(['Loading ' fname]);
                load(fname);
                dat.anim_data(a) = anim_data;
            end
            first_time = 1;
        end
    end

    % exclusion data - loop is here bc this was added later
    if (first_time)
        for a=1:length(dat.anim_data)
            data_dir = [settings.data_path_root filesep settings.animals(a).name filesep 'analysis'];

            fl = dir([data_dir filesep 'flaggedRoisList*mat']);
            [settings.data_path_root filesep settings.animals(a).name  filesep 'session_neuropilone_final_merged' filesep 'find_duplicate_rois_exclude_list.mat'];
            fl2 = dir([settings.data_path_root filesep settings.animals(a).name  filesep 'session_neuropilone_final_merged' filesep 'find_duplicate_rois_exclude_list.mat']);
            if (length(fl) == 1)
                disp(['Processing ' data_dir ' with exclusion file flaggedRoisList*mat format']);
                load([data_dir filesep fl(1).name]);
                excluded_always_ids = [];
                excluded_for_correlations_only_ids = [];

                for i=1:length(flaggedRoisList.keySVs)
                    excluded_always_ids = [excluded_always_ids flaggedRoisList.falsePositiveList_byRoiId{i}];
                    excluded_for_correlations_only_ids = [excluded_for_correlations_only_ids flaggedRoisList.contaminatedList_byRoiId{i}];
                end

                dat.anim_data(a).types_by_id.excluded_always = excluded_always_ids;
                dat.anim_data(a).types_by_id.excluded_for_correlations = union(excluded_for_correlations_only_ids, excluded_always_ids);
            elseif (length(fl2) == 1)
                disp(['Processing ' data_dir ' with exclusion file find_duplicate_rois_exclude_list.mat format']);
                load([settings.data_path_root filesep settings.animals(a).name  filesep 'session_neuropilone_final_merged' filesep fl2(1).name]);

                dat.anim_data(a).types_by_id.excluded_always = exclude_ids;
                dat.anim_data(a).types_by_id.excluded_for_correlations = exclude_ids;
            else
                disp(['No exclusion roi data found for ' settings.animals(a).name]);
            end
        end
    end

    %%% --- updates post-loading
    dat.settings = settings; % update in case chg

    % rerun type assigns in case code changes
    for a=1:length(settings.animals)
        anim_data = assign_types(dat.anim_data(a) );
        dat.anim_data(a) = anim_data;
    end

    % lamniar stuff
    layer_padding = 25;
    disp(['Layer padding is ' num2str(layer_padding) ' um, which is added and subtracted from lalb borders to assign laminar type - so multiply by 2']);
    dat.settings.layer_padding_microns_plus_minus = layer_padding;
    for a=1:length(dat.settings.animals)


        % assign layers
        dat.anim_data(a).types_by_id.L2 = dat.anim_data(a).ids(find(dat.anim_data(a).z_micron_corrected < (dat.settings.animals(a).l2l3_border-layer_padding)));
        dat.anim_data(a).types_by_id.L3 = dat.anim_data(a).ids(find(dat.anim_data(a).z_micron_corrected > (dat.settings.animals(a).l2l3_border+layer_padding) & ...
                                                                    dat.anim_data(a).z_micron_corrected < (dat.settings.animals(a).l3l4_border-layer_padding) ));
        dat.anim_data(a).types_by_id.L4 = dat.anim_data(a).ids(find(dat.anim_data(a).z_micron_corrected > (dat.settings.animals(a).l3l4_border+layer_padding) & ...
                                                                    dat.anim_data(a).z_micron_corrected < (dat.settings.animals(a).l4l5_border-layer_padding) ));
    end

    % normalized depth : assign L2/L3 border a value of 3, L3/L4 4, and then use the micron-to-normalized conversion and apply it to everything outside L3
    for a=1:length(dat.settings.animals)
        linear_scaling = 1/(dat.settings.animals(a).l3l4_border - dat.settings.animals(a).l2l3_border);

        dat.anim_data(a).z_normalized = (dat.anim_data(a).z_micron_corrected - dat.settings.animals(a).l2l3_border)*linear_scaling+3;
    end

    % apply exclusion (cells we omit b/c copycat/overlap or other problems) to all types_by_id
    fn = fieldnames(dat.anim_data(a).types_by_id);
    for a=1:length(dat.settings.animals)
        if (isfield(dat.anim_data(a).types_by_id, exclude_matrix_name))
            excluded_ids = dat.anim_data(a).types_by_id.(exclude_matrix_name);
            disp(['For ' dat.settings.animals(a).name ' removed: ' num2str(length(excluded_ids)) ' out of ' num2str(length(dat.anim_data(a).ids)) ' ; frac: ' num2str(length(excluded_ids)/length(dat.anim_data(a).ids))]);
            for f=1:length(fn)
                if (isempty(strfind(fn{f},'excluded_')))
                    n_pre = length(dat.anim_data(a).types_by_id.(fn{f}));
                    dat.anim_data(a).types_by_id.(fn{f}) = setdiff(dat.anim_data(a).types_by_id.(fn{f}), excluded_ids);
                    n_post = length(dat.anim_data(a).types_by_id.(fn{f}));
                end
            end
        else
            disp(['No exclusion data for ' settings.animals(a).name]);  
        end
    end

% session loop for individual animal
function dat = get_animal_data(data_dir, roi_dir, settings, wcard)

    if (~exist(data_dir,'dir'))
        dat = [];
        disp([data_dir ' does not exist; not loading']);
        return;
    end

    if (nargin < 4) ; wcard = 'an*sess.mat' ; end
    
    % blank
    dat.ids = [];
    dat.subvol_idx = [];

    dat.xy_micron = [];
    dat.z_micron_corrected = [];

    dat.touch_scores = [];
    dat.touch_scores_sh = [];
    dat.inc_touch_scores = [];
    dat.inc_touch_scores_sh = [];
    dat.event_rate = [];
    dat.burst_rate = [];

    % session object stuff
    cd(data_dir);
    fl = dir(wcard);
    for f=1:length(fl)
        load(fl(f).name);
        disp(['Processing ' fl(f).name]);

        % basics
        dat.ids = [dat.ids s.caTSA.ids];
        dat.subvol_idx = [dat.subvol_idx f+0*s.caTSA.ids];

        % encoding
        score_feat = sprintf('dffBased_c%01dAbsMaxKappaZeroNotouch-c%01dAbsMaxKappaZeroNotouchtwoWhiskerModelWithCrossvalScore', settings.whiskers(1), settings.whiskers(2));
        score_feat_sh = sprintf('dffBased_c%01dAbsMaxKappaZeroNotouch-c%01dAbsMaxKappaZeroNotouchtwoWhiskerModelShuffleWithCrossvalScore', settings.whiskers(1), settings.whiskers(2));
        score_structs = s.caTSA.cellFeatures.get(score_feat);
        score_structs_sh = s.caTSA.cellFeatures.get(score_feat_sh);
        if (isempty(score_structs))
            disp('EMPTY score_feat ; aborting');
            return;
        end
        
        dat.touch_scores = [dat.touch_scores score_structs];
        dat.touch_scores_sh = [dat.touch_scores_sh score_structs_sh];

        score_inc_feat = sprintf('dffBased_inclusive_c%01dAbsMaxKappaZeroNotouch-c%01dAbsMaxKappaZeroNotouchtwoWhiskerModelWithCrossvalScore', settings.whiskers(1), settings.whiskers(2));
        score_inc_feat_sh = sprintf('dffBased_inclusive_c%01dAbsMaxKappaZeroNotouch-c%01dAbsMaxKappaZeroNotouchtwoWhiskerModelShuffleWithCrossvalScore', settings.whiskers(1), settings.whiskers(2));
        score_inc_structs = s.caTSA.cellFeatures.get(score_inc_feat);

        score_inc_structs_sh = s.caTSA.cellFeatures.get(score_inc_feat_sh);

        dat.inc_touch_scores = [dat.inc_touch_scores score_inc_structs];
        dat.inc_touch_scores_sh = [dat.inc_touch_scores_sh score_inc_structs_sh];

        event_rate = s.caTSA.cellFeatures.get('eventRateHz');
        burst_rate = s.caTSA.cellFeatures.get('burstRateHz');
        dat.event_rate = [dat.event_rate event_rate];
        dat.burst_rate = [dat.burst_rate burst_rate];

        % xyz - subvolumes
        if (~isempty(roi_dir))
            cd (roi_dir);
            bl = load([roi_dir filesep 'basisList.mat']);
            coms = nan*zeros(2, length(s.caTSA.ids));
            depths = nan*zeros(1, length(s.caTSA.ids));
            for p=1:length(s.caTSA.genParams.roiArray)
                fov_str = s.caTSA.genParams.roiArray{p}(end-9:end-4);
                i = find(strcmp(fov_str, bl.fovIds));
                if (length(i) == 1)
                    o = load(bl.basisFiles{i});
                    rA = o.obj;
                    [x y d] = rA.getPositionForRoisRelativeRefStack(rA.roiIds);

                    vali = find(ismember(rA.roiIds, s.caTSA.ids));
                    valii = find(ismember(s.caTSA.ids, rA.roiIds));

                    depths(valii) = d(vali);
                    coms(:,valii) = [x(vali) ; y(vali)];
                end
            end

            dat.xy_micron = [dat.xy_micron coms];
            dat.z_micron_corrected = [dat.z_micron_corrected depths];

            % back
            cd(data_dir);
        else
            disp(['No ROI dir specified for ' data_dir ' ; data structure will lack position information']);
        end
    end

    
    if (~isempty(roi_dir))
        dir_dat = compute_directionality_index;
        if (sum(dir_dat.ids ~= dat.ids) > 0) 
            disp('PROBLEM in get_two_whikser data -- directionality id mismatch');
            pause;
        end
        dat.w1DI = dir_dat.w1DI;
        dat.w2DI = dir_dat.w2DI;
    end
    
    % clean up -- replace empties with nans
    if (~isempty(dat.touch_scores))
        fnames = fieldnames(dat.touch_scores(1));
        for f=1:length(fnames)
            for i=1:length(dat.touch_scores)
                if (isempty(dat.touch_scores(i).(fnames{f})))
                    dat.touch_scores(i).(fnames{f}) = nan;
                end
                if (isempty(dat.touch_scores_sh(i).(fnames{f})))
                    dat.touch_scores_sh(i).(fnames{f}) = nan;
                end
            end
        end    
    end
    
    if (length(dat.inc_touch_scores) > 0)
        fnames = fieldnames(dat.inc_touch_scores(1));
        for f=1:length(fnames)
            for i=1:length(dat.inc_touch_scores)
                if (isempty(dat.inc_touch_scores(i).(fnames{f})))
                    dat.inc_touch_scores(i).(fnames{f}) = nan;
                end
                if (isempty(dat.inc_touch_scores_sh(i).(fnames{f})))
                    dat.inc_touch_scores_sh(i).(fnames{f}) = nan;
                end
            end
        end    
    end

