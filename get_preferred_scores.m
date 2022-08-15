% 
% returns a touch score where you take the mean for each touch cell across all trial types on which it was responsive significantly (CALL assign_types!)
%
function [prefTouchScores prefBalancedTouchScores] = get_preferred_scores(dat)

    prefTouchScores = {}; % 1/anim
    prefBalanceTouchScores = {};
    for ai=1:length(dat.anim_data)
        w1pi = find(ismember(dat.anim_data(ai).ids, dat.anim_data(ai).types_by_id.w1ProExclusiveTouchTrials));
        w1ri = find(ismember(dat.anim_data(ai).ids, dat.anim_data(ai).types_by_id.w1RetExclusiveTouchTrials));
        w2pi = find(ismember(dat.anim_data(ai).ids, dat.anim_data(ai).types_by_id.w2ProExclusiveTouchTrials));
        w2ri = find(ismember(dat.anim_data(ai).ids, dat.anim_data(ai).types_by_id.w2RetExclusiveTouchTrials));

        cell_anim_scores = 0*dat.anim_data(ai).ids;
        cell_anim_balanced_scores = 0*dat.anim_data(ai).ids;

        for c=1:length(cell_anim_scores)
            n = 0;
            if (ismember(c, w1pi))
                cell_anim_scores(c) = cell_anim_scores(c) + dat.anim_data(ai).touch_scores(c).w1ProExclusiveTouchTrials;
                cell_anim_balanced_scores(c) = cell_anim_balanced_scores(c) + dat.anim_data(ai).touch_scores(c).w1ProExclusiveTouchTrialsBalanced;
                n=n+1;
            end
            if (ismember(c, w1ri))
                cell_anim_scores(c) = cell_anim_scores(c) + dat.anim_data(ai).touch_scores(c).w1RetExclusiveTouchTrials;
                cell_anim_balanced_scores(c) = cell_anim_balanced_scores(c) + dat.anim_data(ai).touch_scores(c).w1RetExclusiveTouchTrialsBalanced;
                n=n+1;
            end
            if (ismember(c, w2pi))
                cell_anim_scores(c) = cell_anim_scores(c) + dat.anim_data(ai).touch_scores(c).w2ProExclusiveTouchTrials;
                cell_anim_balanced_scores(c) = cell_anim_balanced_scores(c) + dat.anim_data(ai).touch_scores(c).w2ProExclusiveTouchTrialsBalanced;
                n=n+1;
            end
            if (ismember(c, w2ri))
                cell_anim_scores(c) = cell_anim_scores(c) + dat.anim_data(ai).touch_scores(c).w2ProExclusiveTouchTrials;
                cell_anim_balanced_scores(c) = cell_anim_balanced_scores(c) + dat.anim_data(ai).touch_scores(c).w2RetExclusiveTouchTrialsBalanced;
                n=n+1;
            end

            if (n > 1)
                cell_anim_scores(c) = cell_anim_scores(c)/n;
                cell_anim_balanced_scores(c) = cell_anim_balanced_scores(c)/n;
            end
        end

        prefTouchScores{ai} = cell_anim_scores;
        prefBalancedTouchScores{ai} = cell_anim_balanced_scores;
    end



