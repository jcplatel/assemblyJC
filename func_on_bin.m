function [out,index_to_avg] = func_on_bin(var, var_to_bin, bin, func)

index_to_avg=zeros(1,length(bin));
out=zeros(1,length(bin));
% here we first need to convert to position to bin data. Instead of the
% vlaue of the position you will get the value of the bin for a given data
% point. Here 1cm -> 1bins 100 cm -> 50bin (because we have 2cm bins)
% this is similar than taking floor(position/binsize) but this is nicer
index_to_avg = discretize(var_to_bin, bin)+1;
index_to_avg(isnan(index_to_avg))=1;
% remove out of bound value (that are below your min or above your max bin)
% that are affected to nan.
% m_nan = ~isnan(index_to_avg);

% accumarray is cool, absically you tell you it to apply a function on var
% for each index_to_avg value. What is nice is that as you can see you can
% very easily change the function you apply here we do the same function
% and we just tell it to calculate a mean or std for each bins

% out = accumarray(index_to_avg(m_nan)', var(m_nan)',[], func);
out = accumarray(index_to_avg', var',[], func);
end