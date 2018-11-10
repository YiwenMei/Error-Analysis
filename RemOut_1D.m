% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 02/02/2018

%% Functionality
% Remove outliers (extremely large/small values) from time series and perform
% linear interpolation to fill the gaps.
% The removal of outliers is based on two evaluation processes:
%  1) Global outlier evaluation: find outliers using the entire time series;
%  2) Local outlier evaluation: find outliers based on a local time window.

%% Input:
%  dn: date number of the time series;
% TSo: time series inputted for outlier removal (T-by-N, where T is the length
%      of the time series and N is the number of time series);
%  Lw: lenght of time window to perform the local outlier evaluation (n time step);
%  a : a ratio used to define the upper bound (mean+a*SD) of TS above which values
%      are regarded as outlier.

%% Output:
% TSi: time series with the removal of outlier and linear interpolation to fill the gap;
% TSn: time series with the removal of outlier only.

function [TSi,TSn]=RemOut_1D(dn,TSo,Lw,a)
%% Global outlier evaluation
MTh=quantile(TSo,.5,1); % median
TSn=abs(TSo-repmat(MTh,length(dn),1)); % detrended TS (flip all values to + side)

xTh=quantile(TSn,.9,1); % maximum without outliers
TSn1=TSn;
TSn1(TSn1>repmat(xTh,length(dn),1))=NaN; % Remove maximum

mTS=nanmean(TSn1,1); % Mean of the stable time series
vTS=nanstd(TSn1,1); % SD of the stable time series
xTS=mTS+a*vTS; % Bound of TS

TSo(TSn-repmat(xTS,length(dn),1)>0)=NaN;

%% Local outlier evaluation
mTS=movmean(TSo,Lw,'omitnan');
xTS=movmax(TSo,Lw,'omitnan');
iTS=movmin(TSo,Lw,'omitnan');
mTS=(mTS*Lw-xTS-iTS)/(Lw-2);
TSn=abs(TSo-mTS); % detrended TS

vTS=movstd(TSo,Lw,'omitnan');
xTS=mTS+a*vTS; % Bound of TS

TSo(TSn-xTS>0)=NaN; % New time series with the removal of extremes
TSn=TSo;

%% Interpolation
TSi=nan(size(TSo));
for vi=1:size(TSo,2)
  v=TSo(:,vi);
  if length(find(isnan(v)))<=size(TSo,1)*.4;
    v(isnan(v))=interp1(dn(~isnan(v)),v(~isnan(v)),dn(isnan(v)));
    TSi(:,vi)=v;
  end
end
end
