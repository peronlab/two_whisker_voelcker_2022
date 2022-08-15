
function dec_dat = get_two_whisker_decode_data (force_redo)
    if (nargin < 1) ; force_redo = 0 ; end
    settings = get_two_whisker_settings;
    fname = [settings.summary_data_root filesep 'two_whisker_per_decode_summary.mat'];
       
    if (~exist(fname, 'file') | force_redo)
        for ai=1:length(settings.animals)
            dec_dat(ai) = get_decode_data_single_animal([settings.data_path_root filesep settings.animals(ai).name filesep 'session_neuropilone_final_merged']);
        end
        save(fname, 'dec_dat')
    else
        global dec_dat;
        if (isempty(dec_dat))
            load(fname);
        end
    end

    % some touchups
    for a=1:length(dec_dat)
        dec_dat(a).scores_fields = {'W1P','W1R','W2P','W2R'};
    end

function dec_dat = get_decode_data_single_animal(datapath)
    odir = pwd;
    cd (datapath);

    cf_fields = {'W1vW2','W1ProW1Ret','W2ProW2Ret','W1vMW','W2vMW'};

    % loop over session objects
    fl = dir('an*sess.mat');
    for f=1:length(fl)
        load(fl(f).name);

        % d-prime between top 1/3 and lowest 1/3 of forces
%disp('GATHER_DPRIME_FORCE_HILO TEMPORARILY DISABLED');        
%ids = s.caTSA.ids ; dffVDKappa = [];  scores = [];
        [ids scores dffVDKappa] = gather_dprime_force_hilo(s);

        dec_dat.subvol(f).ids = ids;
        dec_dat.subvol(f).scores = scores;
        dec_dat.subvol(f).dffVDKappa = dffVDKappa;

        % discrim fields
        for c=1:length(cf_fields)
            cf_name = sprintf('dffBased__twoWhisker%sDiscrimScores', cf_fields{c});
            cf_name_sh = sprintf('dffBased__twoWhisker%sDiscrimScoresShuffled', cf_fields{c});

            new_name = sprintf('%sDiscrimAUC', cf_fields{c});
            new_name_sh = sprintf('%sDiscrimAUCShuffled', cf_fields{c});

            dec_dat.subvol(f).(new_name) = s.caTSA.cellFeatures.get(cf_name);
            dec_dat.subvol(f).(new_name_sh) = s.caTSA.cellFeatures.get(cf_name_sh);
        end
    end

    
function [ids scores dffVDKappa] = gather_dprime_force_hilo(obj)
    dat = get_two_whisker_data(0);
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

    % collect data
    disp(['Processing dataset: ' obj.baseFileName]);
    for i=1:length(obj.caTSA)
        ts = obj.caTSA.dffTimeSeriesArray.getTimeSeriesByIdx(i);
        fitField = sprintf('dffBased_%sAbsMaxKappaZeroNotouch-%sAbsMaxKappaZeroNotouchtwoWhiskerModelAllDataFit', w1Name, w2Name);
        twmFit = obj.caTSA.cellFeatures.get(fitField);
        modelFit = twmFit(i);
        if (i == 1) ; disp('      '); end
        disp([8 8 8 8 8 8 8 sprintf('%06d', i)])

        % dff v force
        try
            dffPredicted = ts.predictDffFromTwoWhiskerModelFit(obj.caTSA.trialIndices, trials, dKappaTS{1}, dKappaTS{2}, modelFit);
        catch me
            disp('FAILED ; skip');
            continue;
        end

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

            if (length(sorti) > 10) % drop degenerates
                g1 = sorti(1:round(length(sorti)/3));
                g2 = sorti(round(2*length(sorti)/3):end);

                d1 = dffVDKappaMat(g1,2);
                d2 = dffVDKappaMat(g2,2);

                dprime_force_hilo = abs(nanmean(d1)-nanmean(d2))/sqrt(0.5*(nanvar(d1) + nanvar(d2)));

                scores(i,e) = dprime_force_hilo;
            else
                disp(['TOO FEW TRIALS: ' usedES{e}.idStr]);
            end

            dffVDKappa{i,e} = dffVDKappaMat;
        end
     end
