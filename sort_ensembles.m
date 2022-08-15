
% sorts ensembles given binary matrix eT where each dim 1 is an ensemble, dim 2 neuron
function [eT sort_eT] = sort_ensembles(eT, M, sort_mode)
    sort_eT = ones(1,size(eT,2));
   
    M = get_full_corrmat(M);

    % pre-arrange by size
    if (strcmp(sort_mode, 'biggest-first'))
        num_entries = sum(eT');
        [irr sidx] = sort(num_entries, 'descend');
        eT = eT(sidx,:);
    elseif (strcmp(sort_mode, 'smallest-first'))
        num_entries = sum(eT');
        num_entries(find(num_entries < 1)) = max(num_entries)+1;
        [irr sidx] = sort(num_entries, 'ascend');
        eT = eT(sidx,:); 
    elseif (strcmp(sort_mode, 'maxcorr-first'))
        for e=1:size(eT,1)
            ensi = find(eT(e,:));
            if (~isempty(ensi))
                mu_corr(e) = nanmean(reshape(M(ensi,ensi),[],1));
            else
                mu_corr(e) = nan;
            end
        end
        [irr sidx] = sort(mu_corr, 'descend');
        eT = eT(sidx,:);
    elseif (strcmp(sort_mode, 'sumcorr-first'))
        for e=1:size(eT,1)
            ensi = find(eT(e,:));
            if (~isempty(ensi))
                sum_corr(e) = nansum(reshape(M(ensi,ensi),[],1));
            else
                sum_corr(e) = nan;
            end
        end
        [irr sidx] = sort(sum_corr, 'descend');
        eT = eT(sidx,:);
    end

    % calculate best ensemble for each neuron based on its correlation to each ensemble's members
    corrs = nan*ones(1,size(eT,1));
    best_ens_per_nrn = nan*ones(1,size(eT,2));
    best_corr = nan*ones(1,size(eT,2));
    for n=1:size(eT,2)
        for e=1:size(eT,1)
            ensi = find(eT(e,:)); 
            ensi = setdiff(ensi,n); % self-exclude
            corrs(e) = nanmean(M(n,ensi));
        
            % if not part of ensemble, cannot be your best
            if (~eT(e,n)) ; corrs(e) = 0 ; end
        end
        [mcorr idx] = max(corrs); 
        if (mcorr > 0)
            best_corr(n) = mcorr;
            best_ens_per_nrn(n) = idx;
        end
    end

    % arrange by ensemble
    si_offs = 0;
    for e=1:size(eT,1)
        ii = find(best_ens_per_nrn == e);

        % within ensemble, sort by best_corr
        [irr sidx] = sort(best_corr(ii),'descend');
        ii = ii(sidx);

        % now put the vector in final output, sort_eT
        sort_eT(si_offs + (1:length(ii))) = ii;
        si_offs = si_offs + length(ii);
    end

    % rejects
    rejecti = find(isnan(best_ens_per_nrn));
    if (~isempty(rejecti))
        sort_eT(si_offs + (1:length(rejecti))) = rejecti;
    end

