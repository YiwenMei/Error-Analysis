% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 11/20/2019

%% Functionality
% This code calculate statistics and error metrics for matched time series pairs
%  and conduct statistical significant tests for some of the metrics.

%% Input
% TS : matched data pairs (N-by-2 array where N is the number of matched pair
%       and the first/second column is for the target/reference data);
%  a : significant level for the statistical significant tests;
% Thr: thresholds used to calculate the contigency statistics.

%% Output
% sts: statistics of target(s) and reference time series including common sample
%       size, mean of target(s) and reference, variance of target(s) and reference;
% ems: error metrics for the target(s) time series including RMS(s), CRMS(s),
%       CC(s), NSE(s), and KGE(s);
% sgs: significant tests results for test of ME/MRE, CRMS/NCRMS, and CC with
%       1 means statistical significant and 0 otherwise (for single target, the
%       result means ME/MRE is significantly different than 0/0, CRMS/NCRMS is
%       significantly different than the standard deviation of reference/1, and
%       CC is significantly different than 0; for multiple targets, result shows
%       the pair-wise comparison, meaning that error metrics of target i is significantly
%       different than that of target j);

% css: contigency statistics for the target time series including percentage
%       of H(s), M(s), F(s), and N(s).

%% Additional note
% sts includes
%  1)size of common samples,    2)mean of target(s), 3)mean of reference,
%  4)variance of target(s), and 5)variance of reference;
% ems includes
%  1)root mean square error,  2)centered root mean square error,
%  3)correlation coefficient, 4)Nash-Sutcliffe efficiency, and
%  5)Kling-Gupta efficiency;
% sgs includes significant tests results of
%  1)mean (relative) error; 2)(normalized) centered root mean square error, and
%  3)correlation coefficient;
% css includes percentage of
%  1)hit, 2)missing, 3)false alarm, and 4)correct negative.

function [sts,ems,sgs,css]=errM(TS,a,Thr)
%% Statistics
N=size(TS,1); % Sample size
mn=mean(TS,1); % Mean of target and reference
vr=var(TS,1,1); % Variance of target and reference
sts=[N mn vr];

%% Error metrics
ems=nan(1,5*(size(TS,2)-1));
css=nan(1,4*(size(TS,2)-1));
sgs=nan(1,3*(size(TS,2)-1));
for i=1:size(TS,2)-1
  Ts=TS(:,[end i]);
  RMS=sqrt(mean(diff(Ts,[],2).^2)); % Root mean square error
  CRMS=std(diff(Ts,[],2),1); % Centered root mean square error
  CC=corr(Ts); % Correlation coefficient
  NSE=1-mean(diff(Ts,[],2).^2)/vr(end); % Nash Sutcliff efficiency
 % Kling-Gupta efficiency
  KGE=1-sqrt((CC-1).^2+(mn(i)/mn(end)-1)^2+(sqrt(vr(i)/vr(end))-1)^2);
  ems(i:(size(TS,2)-1):end)=[RMS CRMS CC(2,1) NSE KGE(2,1)];

% Contigency statistics
  if ~isempty(Thr)
    [th,j]=min(Ts,[],2);
    [th1,~]=max(Ts,[],2);

    r_h=length(find(th>Thr))/N; % Hit rate
    r_n=length(find(th1<=Thr))/N; % Correct negative rate
    r_f=length(find(th<=Thr & th1>Thr & j==1))/N; % False alarm rate
    r_m=length(find(th<=Thr & th1>Thr & j==2))/N; % Missing rate
    css(i:(size(TS,2)-1):end)=[r_h r_m r_f r_n];
  end

%% Significant tests
% Single significant tests
  H1=ttest(Ts(:,2),Ts(:,1),'Alpha',a); % M(R)E sig<> 0 (1) or not (0);
  H2=vartest(diff(Ts,1,2),vr(end),'Alpha',a,'Tail','left'); % (N)CRMS sig< SD_ref/1 (1) or not (0);
  [~,H3]=corr(Ts,'Tail','right'); % CC sig> 0 (1) or not (0)
  H3=H3(2,1)<a;
  sgs(i:(size(TS,2)-1):end)=[H1 H2 H3];
end

% Paired significant tests
if size(TS,2)>2
  K=nchoosek(1:size(TS,2)-1,2);
  H1=nan(1,size(K,1));
  H2=nan(1,size(K,1));
  H3=nan(1,size(K,1));
  for i=1:size(K,1)
    Ts=TS(:,K(i,:))-TS(:,end); % Tag - ref
    H1(i)=ttest2(Ts(:,1),Ts(:,2),'Alpha',a); % M(R)E_a sig<> M(R)E_b (1) or not (0);
    H2(i)=vartest2(Ts(:,1),Ts(:,2),'Alpha',a); % (N)CRMS_a sig<> (N)CRMS_b (1) or not (0);
% CC_a sig<> CC_b (1) or not (0)
    H3(i)=Cor2OverlapTTest(TS(:,end),TS(:,K(i,1)),TS(:,K(i,2)),a);
  end
  sgs=[sgs H1 H2 H3];
end
end

function [H,zs,pv,msg]=Cor2OverlapTTest(Sref,S1,S2,alp)
%% Compute the correlation
cc=corr([Sref S1 S2],'rows','complete');
cr1=cc(1,2); % sample 1 to ref
cr2=cc(1,3); % sample 2 to ref
c12=cc(2,3); % sample 1 to sample 2
N=length(find(~isnan(Sref) & ~isnan(S1) & ~isnan(S2))); % Sample size

%% Hypothesis test
% Compute the test statistic (Meng et al. 1992)
zr1=log((1+cr1)/(1-cr1))/2; % Fisher transform
zr2=log((1+cr2)/(1-cr2))/2;
c_2m=(cr1^2+cr2^2)/2;
f=min([(1-c12)/(2*(1-c_2m)) 1]);
h=1+c_2m*(1-f)/(1-c_2m);
zs=(zr1-zr2)*sqrt((N-3)/(2*(1-c12)*h));

% Compare to the normal distribution
pv=1-normcdf(abs(zs),0,1)+normcdf(-abs(zs),0,1);
H=pv<alp;
if H
  msg='Reject H_0: the two /rho are equaled';
else
  msg='Fail to reject H_0: the two /rho are equaled';
end
end
