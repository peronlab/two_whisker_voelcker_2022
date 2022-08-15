% HERE IS HOLY SCRIPTURE THAT ASSIGNS TYPE - BEHOLD 
function anim_data = assign_types(anim_data)

    % SETTINGS
    n_event_rate_bins = 10; % how many bins to use fo event rate matching
    shuffle_bin_quantile = .99; % this is threshold within-bin
    encoding_thresh_overall = 0.1; % must meet this criteria *and* shuffle criteria
    encoding_thresh_hi_overall = 0.25; % must meet this criteria *and* shuffle criteria

    % prepare touch encoding scores
    enc_fields = {'allTrials', 'w1ExclusiveTouchTrials', 'w2ExclusiveTouchTrials', ...
                      'w1ProExclusiveTouchTrials', 'w1RetExclusiveTouchTrials', 'w2ProExclusiveTouchTrials', 'w2RetExclusiveTouchTrials'};
    for f=1:length(enc_fields)
        anim_data.types_by_id.(enc_fields{f}) = [];
        anim_data.types_by_id.([enc_fields{f} '_hi_thresh']) = [];
    end
    
    % thins requiring partitioning on event rates & subvolume
    e_rates = anim_data.event_rate;
    subvoli = anim_data.subvol_idx;
    usv = unique(subvoli);
    for s=1:length(usv)
        % --- erate partitioning
        sv_vali = find(subvoli == usv(s));
        bini = nan*sv_vali;
        [irr sorti] = sort(e_rates(sv_vali), 'descend');
        d_bin = floor(length(sorti)/n_event_rate_bins);
        bin_idx = 1;
        for bi=1:d_bin:(length(sorti)-d_bin)
            bini(sorti(bi:(bi+d_bin-1))) = bin_idx;
            bin_idx = bin_idx+1;
        end

        % --- significant encoding scores assign
        for f=1:length(enc_fields)
            scores = [anim_data.touch_scores.(enc_fields{f})];
            sh_scores = [anim_data.touch_scores_sh.(enc_fields{f})];
            scores = scores(sv_vali);
            sh_scores = sh_scores(sv_vali);
            sv_ids = anim_data.ids(sv_vali);

            for bi = 1:n_event_rate_bins
                binii = find(bini == bi);
                bin_thresh = quantile(sh_scores(binii),shuffle_bin_quantile);

                overall_thresh = max(bin_thresh, encoding_thresh_overall);
                vali = binii(find (scores(binii) > overall_thresh));
                anim_data.types_by_id.(enc_fields{f}) = [anim_data.types_by_id.(enc_fields{f}) sv_ids(vali)];

                overall_thresh_hi = max(bin_thresh, encoding_thresh_hi_overall);
                vali = binii(find (scores(binii) > overall_thresh_hi));
                field_name = [enc_fields{f} '_hi_thresh'];
                anim_data.types_by_id.(field_name) = [anim_data.types_by_id.(field_name) sv_ids(vali)];
            end
        end
    end

    % "hardcoded" types 
% old - very similar to what we have but did not apply exclusion; keeping for posterity    
%    w1CellIds = anim_data.types_by_id.w1ExclusiveTouchTrials;
%    w2CellIds = anim_data.types_by_id.w2ExclusiveTouchTrials;
%    anim_data.types_by_id.mwCellIds = intersect(w1CellIds,w2CellIds);
%    anim_data.types_by_id.w1CellIds = setdiff(w1CellIds, anim_data.types_by_id.mwCellIds);
%    anim_data.types_by_id.w2CellIds = setdiff(w2CellIds, anim_data.types_by_id.mwCellIds);
%    anim_data.types_by_id.swCellIds = union(anim_data.types_by_id.w1CellIds,  anim_data.types_by_id.w2CellIds);

    w1pCellIds = anim_data.types_by_id.w1ProExclusiveTouchTrials;
    w1rCellIds = anim_data.types_by_id.w1RetExclusiveTouchTrials;
    w2pCellIds = anim_data.types_by_id.w2ProExclusiveTouchTrials;
    w2rCellIds = anim_data.types_by_id.w2RetExclusiveTouchTrials;

    w1OnlyCellIds = setdiff(union(w1pCellIds, w1rCellIds), union(w2pCellIds, w2rCellIds));
    w2OnlyCellIds = setdiff(union(w2pCellIds, w2rCellIds), union(w1pCellIds, w1rCellIds));
    anim_data.types_by_id.w1OnlyCellIds = w1OnlyCellIds;
    anim_data.types_by_id.w2OnlyCellIds = w2OnlyCellIds;
    anim_data.types_by_id.swCellIds = union(w1OnlyCellIds, w2OnlyCellIds);
    anim_data.types_by_id.mwCellIds = intersect(union(w1pCellIds, w1rCellIds), union(w2pCellIds, w2rCellIds));

    % isolate single-whisker cells
    w1pCellIds = intersect(w1pCellIds,  anim_data.types_by_id.swCellIds);
    w1rCellIds = intersect(w1rCellIds,  anim_data.types_by_id.swCellIds);
    w2pCellIds = intersect(w2pCellIds,  anim_data.types_by_id.swCellIds);
    w2rCellIds = intersect(w2rCellIds,  anim_data.types_by_id.swCellIds);

    anim_data.types_by_id.w1pCellIds = w1pCellIds;
    anim_data.types_by_id.w1rCellIds = w1rCellIds;
    anim_data.types_by_id.w2pCellIds = w2pCellIds;
    anim_data.types_by_id.w2rCellIds = w2rCellIds;

    anim_data.types_by_id.w1BidirCellIds = intersect(w1OnlyCellIds, intersect(w1pCellIds, w1rCellIds));
    anim_data.types_by_id.w1UnidirCellIds = intersect(w1OnlyCellIds, setdiff(union(w1pCellIds, w1rCellIds),intersect(w1pCellIds, w1rCellIds)));
    anim_data.types_by_id.w2BidirCellIds = intersect(w2OnlyCellIds, intersect(w2pCellIds, w2rCellIds));
    anim_data.types_by_id.w2UnidirCellIds = intersect(w2OnlyCellIds, setdiff(union(w2pCellIds, w2rCellIds),intersect(w2pCellIds, w2rCellIds)));

    anim_data.types_by_id.unidirCellIds = union(anim_data.types_by_id.w1UnidirCellIds, anim_data.types_by_id.w2UnidirCellIds);
    anim_data.types_by_id.bidirCellIds = union(anim_data.types_by_id.w1BidirCellIds, anim_data.types_by_id.w2BidirCellIds);


