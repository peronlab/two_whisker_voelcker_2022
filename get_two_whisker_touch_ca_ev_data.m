% This will gather # ca events per cell per touch
%   - mean max(Abs(dKappa)) for both whiskers for every trial's pole epoch
%   - dff for each cell for each trial
% This should be used for dFF v. dKappa plots and noise correlation analysis
function te_dat = get_two_whisker_touch_ca_ev_data (force_redo)
    if (nargin < 1) ; force_redo = 0 ;end
    
    settings = get_two_whisker_settings;
    fname = [settings.summary_data_root filesep 'two_whisker_ca_ev_data_summary.mat'];

    if (~exist(fname, 'file') | force_redo)
        for ai=1:length(settings.animals)
            te_dat(ai).sv_dat = get_single_animal([settings.data_path_root filesep settings.animals(ai).name filesep 'session_neuropilone_final_merged']);
        end
        save(fname, 'te_dat')
    else
        global te_dat;
        if (isempty(te_dat))
            load(fname);
        end
    end
    
% name says it all
function sv_dat = get_single_animal(sessdir)
    odir = pwd;
    cd (sessdir);
    fl = dir ('an*sess.mat');

    sv_dat = {};

    time_window_msec = [0 2000]; % in seconds, rel touch

    for f=1:length(fl)
        disp(['Processing ' fl(f).name]);
        load(fl(f).name);
        obj = s;
            
        % 3 types consider - w1, w2, both, none
        w1_t_trials = obj.whiskerBarContactESA.esa{1}.eventTrials;
        w2_t_trials = obj.whiskerBarContactESA.esa{2}.eventTrials;
        mw_t_trials = intersect(w1_t_trials, w2_t_trials);
        w1_t_trials = setdiff(w1_t_trials, mw_t_trials);
        w2_t_trials = setdiff(w2_t_trials, mw_t_trials);
        all_touch_trials = union(w1_t_trials, union(w2_t_trials, mw_t_trials));
        non_touch_trials = setdiff(unique(obj.caTSA.trialIndices), all_touch_trials);

        all_trials = union(all_touch_trials, non_touch_trials);
        trial_type_vec = 0*all_trials;
        for t=1:length(trial_type_vec)
            if (ismember(all_trials(t), non_touch_trials))
                trial_type_vec(t) = 0;
            elseif (ismember(all_trials(t), w1_t_trials))
                trial_type_vec(t) = 1;
            elseif (ismember(all_trials(t), w2_t_trials))
                trial_type_vec(t) = 2;
            else
                trial_type_vec(t) = 3;
            end
        end

        % for every trial, time of first touch OR time of pole appearance
        [irr ia ib] = union( obj.whiskerBarContactESA.esa{1}.eventTimes, obj.whiskerBarContactESA.esa{2}.eventTimes); 
        t_start_trials = [obj.whiskerBarContactESA.esa{1}.eventTrials(ia) obj.whiskerBarContactESA.esa{2}.eventTrials(ib)];
        t_start_times= [obj.whiskerBarContactESA.esa{1}.eventTimes(ia) obj.whiskerBarContactESA.esa{2}.eventTimes(ib)];

        p_start_times = obj.behavESA.getEventSeriesByIdStr('Pole in').eventTimes(1:2:end);
        p_start_trials = obj.behavESA.getEventSeriesByIdStr('Pole in').eventTrials(1:2:end);
        
        eval_start_times = nan*trial_type_vec;
        for t=1:length(all_trials)
            if (trial_type_vec(t) == 0)
                candi_ti = find(p_start_trials == all_trials(t));
                eval_start_times(t) = min(p_start_times(candi_ti));
            else
                candi_ti = find(t_start_trials == all_trials(t));
                eval_start_times(t) = min(t_start_times(candi_ti));
            end
        end

        % cell by cell
        cell_ids = obj.caTSA.caPeakEventSeriesArray.ids;
        ev_per_trial = zeros(length(cell_ids), length(all_trials));
        ev_per_trial_sh = zeros(length(cell_ids), length(all_trials));
        for c=1:length(cell_ids)
            e_times = obj.caTSA.caPeakEventSeriesArray.esa{c}.eventTimes(1:2:end);
            e_times_sh = sort(obj.caTSA.time(randperm(length(obj.caTSA.time), length(e_times))), 'ascend'); % pick same number of times, randomly
            for t=1:length(all_trials)
                ev_per_trial(c,t) = length(find(e_times >= eval_start_times(t)+time_window_msec(1) & e_times < eval_start_times(t)+time_window_msec(2)));
                ev_per_trial_sh(c,t) = length(find(e_times_sh >= eval_start_times(t)+time_window_msec(1) & e_times_sh < eval_start_times(t)+time_window_msec(2)));
            end
        end

        % output
        sv_dat{f}.source_file = fl(f).name;
        sv_dat{f}.cell_ids = cell_ids;
        sv_dat{f}.trial_type_vec = trial_type_vec;
        sv_dat{f}.readme = 'trial type 0 = no touch ; 1 = w1 ; 2 = w2 ; 3 = mw';
        sv_dat{f}.ev_per_trial = ev_per_trial;
        sv_dat{f}.ev_per_trial_sh = ev_per_trial_sh;
        sv_dat{f}.overall_burst_rate_hz = obj.caTSA.cellFeatures.get('burstRateHz');
        sv_dat{f}.overall_event_rate_hz = obj.caTSA.cellFeatures.get('eventRateHz');
    end

    cd (odir);

