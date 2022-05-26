 function [mu, error, bins2] = KS_bindata_mean_20140916(x, y, numbins);
 %function [mu, bins] = bindata(x, y, numbins);
 % bins the data y according to x and returns the bins and the average
 % value of y for that bin

maxval=max(x);
 bins = linspace(0, maxval, numbins);
  bins2=bins+0.5*(bins(2)-bins(1));
 %bins = linspace(3000, 20000, numbins);
 
 [n,bin] = histc(x, bins);
 mu = NaN*zeros(size(bins));
 error = NaN*zeros(size(bins));
 for k = [1:numbins],
 ind = find(bin==k);
 %if (~isempty(ind))
 if (sum(bin==k)>=2)
 mu(k) = nanmean(y(ind));
 error(k)=nanstd(y(ind))/sqrt(sum(~isnan(y(ind))));
%error(k)=nanstd(y(ind));
 end
 end