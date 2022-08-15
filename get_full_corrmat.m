function [M nM] = get_full_corrmat(M)
    for i1=1:size(M,1)    
        M(i1,i1) = 1;
        for i2=i1+1:size(M,2)    
            M(i2,i1) = M(i1,i2);
        end
    end
    nM = M; 
    for i1=1:size(M,1) % NaN the diagonal
        nM(i1,i1) = nan;
    end

