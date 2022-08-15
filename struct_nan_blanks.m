function s = struct_nan_blanks (s)
% 
% This will remove [] values in array of structures, replacing with nan
%
% USAGE:
%
%   s = struct_nan_blanks(s);
%
    fn = fieldnames(s);
    N = length(s);

    for n=1:N
        for f=1:length(fn)
            if (isempty(s(n).(fn{f}))) ; s(n).(fn{f}) = nan ; end
            
        end
    end


