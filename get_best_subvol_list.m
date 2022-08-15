% returns subvolume that has most cells for selected layer; returns animal X layer matrix
function [subvol_idx l_tag] = get_best_subvol_list (dat)
    l_tag = {'L2','L3','L4'};
    subvol_idx = nan*zeros(length(dat.anim_data), length(l_tag));

    for ani=1:length(dat.anim_data)
        l_ids = {dat.anim_data(ani).types_by_id.L2, dat.anim_data(ani).types_by_id.L3, dat.anim_data(ani).types_by_id.L4};
        for l=1:length(l_tag)
            % figure out subvol
            usv = unique(dat.anim_data(ani).subvol_idx);
            n_cells_per_subvol = zeros(1,length(usv));
            for s=1:length(usv)
                svix = find(dat.anim_data(ani).subvol_idx == usv(s));
                subvol_ids = dat.anim_data(ani).ids(svix);

                n_cells_per_subvol(s) = length(intersect(subvol_ids, l_ids{l}));
            end
%disp([dat.settings.animals(ani).name ' ' l_tag{l} ' ' num2str(n_cells_per_subvol)]);            
            [irr svi] = max (n_cells_per_subvol);

            % overrides
            if (strcmp(dat.settings.animals(ani).name, 'an283544') & strcmp(l_tag{l},'L2')) ; svi = 2; disp('OVERRIDE an283544 L2 subvolume - use 2, not 1'); end
            if (strcmp(dat.settings.animals(ani).name, 'an274688') & strcmp(l_tag{l},'L3')) ; svi = 4; disp('OVERRIDE an274688 L3 subvolume - use 4, not 5'); end

            subvol_idx(ani,l) = svi;
        end
    end

