% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/13/2019

%% Functionality
% This code calculate statistics of matched target and reference value pairs and
%  their error metrics.

%% Input
% TS : matched data pairs (N-by-2 array where N is the number of matched pair
%      and the first/second column is for the target/reference data);
% Thr: thresholds used to calculate the contigency statistics (set it to [] if
%      no need to calculate).

%% Output
% sts: statistics of target and reference time series;
% ems: error metrics for the target time series;
% css: contigency statistics for the target time series (if Thr is [], ccs is
%      []).

%% Additional note
% sts includes
%  1)sample size of matched target and reference time series pair,
%  2)mean of target time series,        3)mean of reference time series,
%  4)variance of target time series and 5)variance of reference time series;
% ems includes
%  1)root mean square error,  2)centered root mean square error,
%  3)correlation coefficient, 4)Nash-Sutcliffe efficiency and
%  5)Kling-Gupta efficiency;
% css includes rate of 1)hit, 2)missing, 3)false alarm, and 4)correct negative.

function [sts,ems,css]=errM(TS,Thr)
%% Statistics
N=size(TS,1); % Sample size
mn=mean(TS,1); % Mean of target and reference
vr=var(TS,1,1); % Variance of target and reference
sts=[N mn vr];

%% Error metrics
RMS=sqrt(mean(diff(TS,[],2).^2)); % Root mean square error
CRMS=std(diff(TS,[],2),1); % Centered root mean square error
CC=corr(TS); % Correlation coefficient
NSE=1-mean(diff(TS,[],2).^2)/vr(2); % Nash Sutcliff efficiency
KGE=1-sqrt((CC-1).^2+(mn(1)/mn(2)-1)^2+(sqrt(vr(1)/vr(2))-1)^2); % Kling-Gupta efficiency
ems=[RMS CRMS CC(2,1) NSE KGE(2,1)];

%% Contigency statistics
if ~isempty(Thr)
  [th,i]=min(TS,[],2);
  [th1,~]=max(TS,[],2);

  r_h=length(find(th>Thr))/N; % Hit rate
  r_n=length(find(th1<=Thr))/N; % Correct negative rate
  r_m=length(find(th<=Thr & th1>Thr & i==1))/N; % Missing rate
  r_f=length(find(th<=Thr & th1>Thr & i==2))/N; % False alarm rate
  css=[r_h r_m r_f r_n];
else
  css=[];
end
end
