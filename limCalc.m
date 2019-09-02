function [lim] = limCalc(vals,coeff,N)
%limCalc Calculates the limit by choice of the user
%   vals is a vector or matrix with graphs points
%   coeff is a 2 index vector with the limit coeff. For example [0.95 1.1]
%   N is the number of digits after the decimal point
%   if there's an error in the limits it returns the minium value in lim(1)
%   and lim(2) equals lim(1)+1
    
    coeff = coeff*10^N;
    lim = [round(coeff(1)*min(vals,[],'all'))/10^N round(coeff(2)*max(vals,[],'all'))/10^N];
    if length(lim)>1
        if lim(1) >= lim(2)
            lim(1) = min(vals,[],'all');
            lim(2) = lim(1) + 1;
        end
        if lim(1)>=0 && lim(1) <= 0.2
            lim(1) = 0;
        end
    else
        lim = [0 1];
    end
end

