function [mw_vali all_vali g1_vali g2_vali sw_vali] = get_subgroups(mode, dat, ani)
    if (strcmp(mode, 'w1w2'))
        all_vali = find(ismember(dat.anim_data(ani).ids, union(dat.anim_data(ani).types_by_id.mwCellIds,dat.anim_data(ani).types_by_id.swCellIds)));

        mw_vali = find(ismember(dat.anim_data(ani).ids, dat.anim_data(ani).types_by_id.mwCellIds));
        g1_vali = find(ismember(dat.anim_data(ani).ids, dat.anim_data(ani).types_by_id.swCellIds));
        g2_vali = find(ismember(dat.anim_data(ani).ids, union( dat.anim_data(ani).types_by_id.w2pCellIds,  dat.anim_data(ani).types_by_id.w2rCellIds))) ;
        sw_vali = union(g1_vali, g2_vali);

    elseif (strcmp(mode, 'dirboth'))

        all_vali = find(ismember(dat.anim_data(ani).ids, union(dat.anim_data(ani).types_by_id.mwCellIds,dat.anim_data(ani).types_by_id.swCellIds)));
        mw_vali = find(ismember(dat.anim_data(ani).ids, dat.anim_data(ani).types_by_id.bidirCellIds));
        g1_vali = find(ismember(dat.anim_data(ani).ids, union( dat.anim_data(ani).types_by_id.w1pCellIds,  dat.anim_data(ani).types_by_id.w2pCellIds))) ;
        g2_vali = find(ismember(dat.anim_data(ani).ids, union( dat.anim_data(ani).types_by_id.w1rCellIds,  dat.anim_data(ani).types_by_id.w2rCellIds))) ;
        sw_vali = find(ismember(dat.anim_data(ani).ids, dat.anim_data(ani).types_by_id.unidirCellIds));
    end
