% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 02/12/2019

%% Functionality
% This code performs error analysis between target and reference data. Target
% and reference data are inputted as 2-D image. 

%% Input
% Ftg/Frf: list of full name of target/reference images;
% Ntg/Nrf: field name of target/reference image for .nc, .nc4, .hdf, .hdf5 and
%          .mat format files (for .mat file, f_tg is the Matlab variable name);
% Vtg/Vrf: no data value of target/reference image;
% Gtg/Grf: x and y boundary and resolution of the target/referece image (2-by-3
%          array where row 1 stores the left, right boundary and resolution of
%          x and row 2 stores the top, bottom boundary and resolution of y);
% Ttg/Trf: date number of the target/reference images;
% Rtg/Rrf: time resolution of target/reference image;
% Ctg/Crf: time convention of target/reference image;
%   cf   : unit conversion factor for target data;
%  Z_thr : a threshold value used to calculate the contigency statistics (set
%          it to "[]" if no need to output the contigency statistics);
%   thr  : percentage of nodata value for a block where block with percent of
%          nodata value greater than this will be assigned the nodata value.

%% Output
% TS.STs: statistics of target and reference map of every time step;
% RS.STs: statistics of target and reference time series of every pixel;
% TS.EMs: error metrics derived between target and reference map of every time
%         step;
% RS.EMs: error metrics derived between target and reference time series of every
%         pixel;
% TS.CSs: contigency statistics given threshold value "thr" of target and reference
%         map of every time step;
% RS.CSs: contigency statistics given threshold value "thr" of target and reference
%         time series of every pixel;

%% Additional note
% STs includes time series and maps of
%  1)sample size of matched target and reference time series pair,
%  2)mean of target time series,        3)mean of reference time series,
%  4)variance of target time series and 5)variance of reference time series;
% EMs includes time series and maps of
%  1)root mean square error,  2)centered root mean square error,
%  3)correlation coefficient, 4)Nash-Sutcliffe efficiency and
%  5)Kling-Gupta efficiency;
% CSs includes time series and maps of rate of
%  1)hit, 2)missing, 3)false alarm, and 4)correct negative.

function [TS,RS]=comp_RS(Ftg,Ntg,Vtg,Gtg,Ttg,Rtg,Ctg,Frf,Nrf,Vrf,Grf,Trf,Rrf,Crf,fRtg,cf,Z_thr,thr,wkpth)
fprintf('\nStart');

Rrf=Rrf/24; % Convert from hour to day
Rtg=Rtg/24;
if strcmp(Ctg,'f') % Make the target record to the centered time convention
  Ttg=Ttg-Rtg/2;
elseif strcmp(Ctg,'b')
  Ttg=Ttg+Rtg/2;
end

rs=Grf(1,3)/Gtg(1,3);

Na=[];
Za_tg=[];
Za_rf=[];
Za2_tg=[];
Za2_rf=[];
Na_h=[];
Na_m=[];
Na_f=[];
Na_n=[];
pZa=[];
dZa=[];
dZa2=[];

STs=nan(size(Frf,1),5);
EMs=nan(size(Frf,1),5);
CSs=nan(size(Frf,1),4);
for t=1:size(Frf,1)
%% Reference image
% Read the image
  Z_rf=read2Dvar({Frf(t,:),Vrf,Inf,-Inf,Nrf});
% Match reference resolution to target
  if rs>1 % If target is finer fRtg is needed
    if ~isempty(fRtg) % Match the ref resolution to target
      Z_rf=imresize(Z_rf,rs,fRtg);
      Grf(:,3)=Gtg(:,3);
    end
  elseif rs<1 % If reference is finer
    Z_rf(isnan(Z_rf))=Vrf;
    idn=fullfile(wkpth,sprintf('id_rs%d_%d.mat',Grf(1,3),Gtg(1,3)));
    Z_rf=resizeimg(Z_rf,Vrf,Gtg(1,:),Gtg(2,:),idn,thr,Grf(1,:),Grf(2,:));
    Z_rf(Z_rf==Vrf)=NaN;
    Grf=Gtg;
  end

%% Target image
  if strcmp(Crf,'f') % Forward
    ftg=Ftg(Ttg>=Trf(t)-Rrf & Ttg<Trf(t),:);
  elseif strcmp(Crf,'c') % Centered
    ftg=Ftg(Ttg>=Trf(t)-Rrf/2 & Ttg<Trf(t)+Rrf/2,:);
  elseif strcmp(Crf,'b') % Backward
    ftg=Ftg(Ttg>=Trf(t) & Ttg<Trf(t)+Rrf,:);
  end

  if ~isempty(ftg)
% Read the image
    Z_tg=[];
    N=[];
    for ti=1:size(ftg,1)
      z_tg=read2Dvar({ftg(ti,:),Vtg,Inf,-Inf,Ntg});
      Z_tg=nansum(cat(3,Z_tg,z_tg),3);
      k=double(~isnan(z_tg));
      N=nansum(cat(3,N,k),3); % All NaN means 0
    end

% Aggregate target image
    Z_tg=cf*Z_tg./N; % Convert to unit of reference
    clear z_tg N
    Z_tg(isnan(Z_tg))=Vtg;
    idn=fullfile(wkpth,sprintf('id_rs%d_%d.mat',Gtg(1,3),Grf(1,3)));
    Z_tg=resizeimg(Z_tg,Vtg,Grf(1,:),Grf(2,:),idn,thr,Gtg(1,:),Gtg(2,:));
    Z_tg(Z_tg==Vtg)=NaN;

  else
    Z_tg=nan(size(Z_rf));
  end

%% Error analysis
% Time series
  k=~isnan(Z_rf) & ~isnan(Z_tg);
  N=length(find(k)); % Sample size

  if N>size(Z_rf,1)*size(Z_rf,2)*(1-thr)
    [sts,ems,css]=errM_TS(Z_rf,Z_tg,Z_thr);

    STs(t,:)=sts; % Statistics
    EMs(t,:)=ems; % Error metrics
    if ~isempty(Z_thr)
      CSs(t,:)=css; % Contigency statistics
    end
  end

% 2D map
  Z_rf(~k)=NaN;
  Z_tg(~k)=NaN;

  Na=sum(cat(3,Na,double(k)),3);
  Za_rf=nansum(cat(3,Za_rf,Z_rf),3);
  Za_tg=nansum(cat(3,Za_tg,Z_tg),3);
  Za2_rf=nansum(cat(3,Za2_rf,Z_rf.^2),3);
  Za2_tg=nansum(cat(3,Za2_tg,Z_tg.^2),3);
  pZa=nansum(cat(3,pZa,Z_rf.*Z_tg),3);
  dZa=nansum(cat(3,dZa,Z_tg-Z_rf),3);
  dZa2=nansum(cat(3,dZa2,(Z_tg-Z_rf).^2),3);

  if ~isempty(Z_thr)
    Na_h=sum(cat(3,Na_h,double(Z_tg>Z_thr & Z_rf>Z_thr)),3);
    Na_m=sum(cat(3,Na_m,double(Z_tg<=Z_thr & Z_rf>Z_thr)),3);
    Na_f=sum(cat(3,Na_f,double(Z_tg>Z_thr & Z_rf<=Z_thr)),3);
    Na_n=sum(cat(3,Na_n,double(Z_tg<=Z_thr & Z_rf<=Z_thr)),3);
  end

% Progress indicator
  if rem(t,floor(size(Frf,1)/100))==0
    fprintf('--');
    if rem(t,25*floor(size(Frf,1)/100))==0
      fprintf('%i%%\n',round(t/size(Frf,1)*100));
    end
  end
end
k=Na>=20/Rrf;

% Statistics
Na(~k)=NaN; % Sample size
Za_rf(~k)=NaN;
Za_tg(~k)=NaN;
Za2_rf(~k)=NaN;
Za2_tg(~k)=NaN;
m_rf=Za_rf./Na; % Mean of reference time series
m_tg=Za_tg./Na; % Mean of target time series
v_rf=Za2_rf./Na-m_rf.^2; % Variance of reference time series
v_tg=Za2_tg./Na-m_tg.^2; % Variance of target time series

% Error metrics
dZa2(~k)=NaN;
dZa(~k)=NaN;
pZa(~k)=NaN;
RMS=sqrt(dZa2./Na); % Root mean square error
CRMS=sqrt(dZa2./Na-(dZa./Na).^2); % Centered root mean square error
CC=(pZa./Na-m_tg.*m_rf)./sqrt(v_rf.*v_tg); % Correlation coefficient
NSE=1-dZa2./(Na.*v_rf); % Nash Sutcliff efficiency
KGE=1-sqrt((CC-1).^2+(m_tg./m_rf-1).^2+(sqrt(v_tg)./sqrt(v_rf).*m_rf./m_tg-1).^2); % Kling-Gupta efficiency

% Contigency statistics
if ~isempty(Z_thr)
  Na_h(~k)=NaN;
  Na_m(~k)=NaN;
  Na_f(~k)=NaN;
  Na_n(~k)=NaN;
  Na_h=Na_h./Na; % Hit rate
  Na_m=Na_m./Na; % Missing rate
  Na_f=Na_f./Na; % False alarm rate
  Na_n=Na_n./Na; % Correct negative rate

  RS=struct('STs',cat(3,Na,m_tg,m_rf,v_tg,v_rf),'EMs',cat(3,RMS,CRMS,CC,NSE,KGE),'CSs',cat(3,Na_h,Na_m,Na_f,Na_n));
  TS=struct('STs',STs,'EMs',EMs,'CSs',CSs);
else
  RS=struct('STs',cat(3,Na,m_tg,m_rf,v_tg,v_rf),'EMs',cat(3,RMS,CRMS,CC,NSE,KGE));
  TS=struct('STs',STs,'EMs',EMs);
end

delete([wkpth 'id*.mat']);
fprintf('Done!\n');
end
