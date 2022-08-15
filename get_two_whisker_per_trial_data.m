% This will gather per trial data, meaning
%   - mean max(Abs(dKappa)) for both whiskers for every trial's pole epoch
%   - dff for each cell for each trial
% This should be used for dFF v. dKappa plots and noise correlation analysis
function pt_dat = get_per_trial_two_whisker_data (force_redo)
    if (nargin < 1) ; force_redo = 0 ;end

    settings = get_two_whisker_settings;
    fname = [settings.summary_data_root filesep 'two_whisker_per_trial_summary.mat'];

    if (~exist(fname, 'file') | force_redo)
        for ai=1:length(settings.animals)
            pt_dat(ai) = get_single_animal([settings.data_path_root filesep settings.animals(ai).name filesep 'session_neuropilone_final_merged']);
        end
        save(fname, 'pt_dat')
    else
        global pt_dat;
        if (isempty(pt_dat))
            load(fname);
        end
    end
    

function dat = get_single_animal(sessdir)
    odir = pwd;
    cd (sessdir);
    fl = dir ('an*sess.mat');

    dat = [];
    dat.ids = [];
    dat.subvol = [];

    timeWindow = [-3 5];
    timeWindowForTouchAndKappaData = [0 2];
    timeWindowForPreTouchData = [-3 0];
    excludeTimeWindow = [-1 1];
    timeUnit = pldo.timeSeries.second;
    allowOverlap = 0;

    for f=1:length(fl)
        disp(['Processing ' fl(f).name]);

        % the meat
        load(fl(f).name);
        obj = s;
        [trials touchES] = obj.getWhiskerTouchTrialTypes;
        es = obj.behavESA.getEventSeriesByIdStr('Pole movement');
        nes = es.copy();
        nes.type = 1;
        L = length(es.eventTrials);
        nes.eventTrials = es.eventTrials(1:2:L);
        nes.eventTimes = es.eventTimes(1:2:L);
        nes.eventTimesRelTrialStart = es.eventTimesRelTrialStart(1:2:L);

        w1Name = obj.whiskerTag{1};
        w2Name = obj.whiskerTag{2};

        dKappaTS{1} = obj.derivedDataTSA.getTimeSeriesById(20112);
        dKappaTS{2} = obj.derivedDataTSA.getTimeSeriesById(20113); 

        % administrative stuff -- cell ids, subvolume id
        dat.ids = [dat.ids obj.caTSA.ids];
        dat.subvol = [dat.subvol 0*obj.caTSA.ids+f];

        dat.per_subvol(f).trials = trials;
        dat.per_subvol(f).kappa_w1 = nan*zeros(1, nes.length);
        dat.per_subvol(f).kappa_w2 = nan*zeros(1, nes.length);
        dat.per_subvol(f).dff_pole_epoch = nan*zeros(obj.caTSA.length, nes.length);
        dat.per_subvol(f).dff_pre_pole_epoch = nan*zeros(obj.caTSA.length, nes.length);

        % dkappa for w1, w2
        [dataMat timeMat idxMat plotTimeVec ieIdxVec] = dKappaTS{1}.getValuesAroundEvents(nes, timeWindow, timeUnit, allowOverlap, [], excludeTimeWindow);
        tidxdff = find(plotTimeVec >= timeWindowForTouchAndKappaData(1) & plotTimeVec < timeWindowForTouchAndKappaData(2));
        tidxdffNotouch = find(plotTimeVec >= timeWindowForPreTouchData(1) & plotTimeVec < timeWindowForPreTouchData(2));
        vali = find(idxMat > 0);
        for t=1:size(dataMat,1)
            dat.per_subvol(f).kappa_w1(t) = nanmean(dataMat(t, tidxdff));
        end

        kappaDataMat = nan*dataMat;
        kappaDataMat(vali) = dKappaTS{2}.value(idxMat(vali));
        for t=1:size(kappaDataMat,1)
            dat.per_subvol(f).kappa_w2(t) = nanmean(kappaDataMat(t, tidxdff));
        end

        for c=1:obj.caTSA.length
            dffDataMat = nan*dataMat;
            dffDataMat(vali) = obj.caTSA.dffTimeSeriesArray.valueMatrix(c,idxMat(vali));
            for t=1:size(dffDataMat,1)
                % stim period mean dFF for each trial/cell
                dat.per_subvol(f).dff_pole_epoch(c,t) = nanmean(dffDataMat(t, tidxdff));

                % prestim period mean dFF for each trial/cell
                dat.per_subvol(f).dff_pre_pole_epoch(c,t) = nanmean(dffDataMat(t, tidxdffNotouch));
             end
         end

         if (length(dat.per_subvol(f).trials.allTrials) ~= size(dat.per_subvol(f).kappa_w1,2))
            disp('  PROBLEM : trial length not same as allTrials . . . look into this');
            pause;
         end
    end

    cd (odir);

