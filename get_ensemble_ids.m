function ens_ids = get_ensemble_ids
    min_touch_cell_frac_touch_ensemble = 0.5; % @ least 50%
    params.exclude_matrix_name = 'excluded_for_correlations';
    params.force_redo = 0;
    global dat;
    global ens_dat;
    if (isempty(dat))
        dat = get_two_whisker_data(params);
    end    
    if (isempty(ens_dat))
        ens_dat = get_two_whisker_ensemble_data;
    end
    [subvol_idx l_tag] = get_best_subvol_list (dat);

    for ani=1:length(ens_dat.anim)
        l_ids = {dat.anim_data(ani).types_by_id.L2, dat.anim_data(ani).types_by_id.L3, dat.anim_data(ani).types_by_id.L4};
        t_ids = union(dat.anim_data(ani).types_by_id.mwCellIds, union(dat.anim_data(ani).types_by_id.unidirCellIds,dat.anim_data(ani).types_by_id.bidirCellIds));
        
        ens_ids(ani).touch_ensemble_ids = [];
        ens_ids(ani).non_touch_ensemble_ids = [];
        
        for li=1:length(l_tag)
            svi = subvol_idx(ani, li);
            
            ens_member_mat = ens_dat.anim(ani).subvol(svi).all_corr.ens_matrix;
            idx_touch = find(ismember(ens_dat.anim(ani).subvol(svi).all_corr.ids, t_ids));

            touch_ensi = [];
            nontouch_ensi = [];
            for e=1:size(ens_member_mat,1)
                idx_ens = find(ens_member_mat(e,:));
                frac_t_cells_per_ens= length(intersect(idx_touch,idx_ens))/length(idx_ens);

                if (frac_t_cells_per_ens >= min_touch_cell_frac_touch_ensemble)
                    touch_ensi = [touch_ensi idx_ens];
                else
                    nontouch_ensi = [nontouch_ensi idx_ens];
                end
            end

            % store for now n later
            ens_ids(ani).touch_ensemble_ids = [ens_ids(ani).touch_ensemble_ids ens_dat.anim(ani).subvol(svi).all_corr.ids(touch_ensi)];
            ens_ids(ani).non_touch_ensemble_ids = [ens_ids(ani).non_touch_ensemble_ids ens_dat.anim(ani).subvol(svi).all_corr.ids(nontouch_ensi)];
        end
    end

