% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 08/20/2018

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
%          nodata value greater than this will be assigned the nodata value;

function [STs,EMs,CSs]=comp_RS(Ftg,Ntg,Vtg,Gtg,Ttg,Rtg,Ctg,Frf,Nrf,Vrf,Grf,Trf,Rrf,Crf,cf,Z_thr,thr)
Rrf=Rrf/24; % Convert from hour to day
Rtg=Rtg/24;

if strcmp(Ctg,'f') % Make the target record to the centered time convention
  Ttg=Ttg-Rtg/2;
elseif strcmp(Ctg,'b')
  Ttg=Ttg+Rtg/2;
end

% Target/reference coordinate
xtg=Gtg(1,1)+Gtg(1,3)/2:Gtg(1,3):Gtg(1,2)-Gtg(1,3)/2;
ytg=Gtg(2,1)-Gtg(2,3)/2:-Gtg(2,3):Gtg(2,2)+Gtg(2,3)/2;

xrf=Grf(1,1)+Grf(1,3)/2:Grf(1,3):Grf(1,2)-Grf(1,3)/2;
yrf=Grf(2,1)-Grf(2,3)/2:-Grf(2,3):Grf(2,2)+Grf(2,3)/2;

Na=zeros(length(yrf),length(xrf));
Za_tg=zeros(length(yrf),length(xrf));
Za_rf=zeros(length(yrf),length(xrf));
Za2_tg=zeros(length(yrf),length(xrf));
Za2_rf=zeros(length(yrf),length(xrf));

Na_h=zeros(length(yrf),length(xrf));
Na_m=zeros(length(yrf),length(xrf));
Na_f=zeros(length(yrf),length(xrf));
Na_n=zeros(length(yrf),length(xrf));

pZa=zeros(length(yrf),length(xrf));
dZa=zeros(length(yrf),length(xrf));
dZa2=zeros(length(yrf),length(xrf));

for t=1:size(Frf,1)
%% Reference image
% Read the image
  [wpth,~,fex]=fileparts(Frf(t,:));
  if strncmp(fex,'.tif',4) % compatable for .tiff
    Z_rf=double(imread(Frf(t,:)));
  elseif strncmp(fex,'.nc4',3) % compatable for .nc
    Z_rf=double(ncread(Frf(t,:),Nrf));
  elseif strncmp(fex,'.hdf',4) % compatable for .hdf5
    Z_rf=double(hdfread(Frf(t,:),Nrf));
  elseif strcmp(fex,'.asc') || strcmp(fex,'.txt')
    Z_rf=double(dlmread(Frf(t,:),'',6,0));
  else
    load(Frf(t,:),Nrf);
    Z_rf=Nrf;
  end
  Z_rf(Z_rf==Vrf)=NaN;

%% Target image
  if strcmp(Crf,'f') % Forward
    ftg=Ftg(Ttg>=Trf(t)-Rrf & Ttg<Trf(t),:);
  elseif strcmp(Crf,'c') % Centered
    ftg=Ftg(Ttg>=Trf(t)-Rrf/2 & Ttg<Trf(t)+Rrf/2,:);
  elseif strcmp(Crf,'b') % Backward
    ftg=Ftg(Ttg>=Trf(t) & Ttg<Trf(t)+Rrf,:);
  end

% Read the image
  Z_tg=nan(length(ytg),length(xtg),1);
  N=nan(length(ytg),length(xtg),1);
  for ti=1:size(ftg,1)
    [~,~,fex]=fileparts(ftg(ti,:));
    if strncmp(fex,'.tif',4) % compatable for .tif & .tiff
      z_tg=double(imread(ftg(ti,:)));
    elseif strncmp(fex,'.nc4',3) % compatable for .nc & .nc4
      z_tg=double(ncread(ftg(ti,:),Ntg));
    elseif strncmp(fex,'.hdf',4) % compatable for .hdf & .hdf5
      z_tg=double(hdfread(ftg(ti,:),Ntg));
    elseif strcmp(fex,'.asc') || strcmp(fex,'.txt')
      z_tg=double(dlmread(ftg(ti,:),'',6,0));
    else
      load(Ftg(t,:),Ntg);
      z_tg=Ntg;
    end
    z_tg(z_tg==Vtg)=NaN;
    Z_tg=nansum(cat(3,Z_tg,z_tg),3);
    k=double(~isnan(z_tg));
    N=nansum(cat(3,N,k),3); % All NaN means 0
  end
% Aggregate target image
  Z_tg=cf*Z_tg./N; % Convert to unit of reference
  clear N k z_tg
  Z_tg(isnan(Z_tg))=Vtg;
  Z_tg=resizeimg(Z_tg,Gtg(1,:),Grf(1,:),Gtg(2,:),Grf(2,:),fullfile(wpth,'id.mat'),thr,Vtg);
  Z_tg(Z_tg==Vtg)=NaN;

%% Error analysis
% Time series
  k=~isnan(Z_rf) & ~isnan(Z_tg);

% Statistics
  N=length(find(k)); % Sample size
  if N>=size(Z_rf,1)*size(Z_rf,2)*(1-thr)
    m_tg=mean(Z_tg(k)); % Mean of target map
    m_rf=mean(Z_rf(k)); % Mean of reference map
    v_tg=var(Z_tg(k),1); % Variance of target map
    v_rf=var(Z_rf(k),1); % Variance of reference map
    STs.TS(t,:)=[N m_tg m_rf v_tg v_rf];

% Error metrics
    RMS=sqrt(mean((Z_tg(k)-Z_rf(k)).^2)); % Root mean square error
    CRMS=std(Z_tg(k)-Z_rf(k),1); % Centered root mean square error
    CC=corr(Z_tg(k),Z_rf(k)); % Correlation coefficient
    NSE=1-sum((Z_tg(k)-Z_rf(k)).^2)/(N*v_rf); % Nash Sutcliff efficiency
    KGE=1-sqrt((CC-1)^2+(m_tg/m_rf-1)^2+(sqrt(v_tg)/sqrt(v_rf)*m_rf/m_tg-1)^2); % Kling-Gupta efficiency
    EMs.TS(t,:)=[RMS CRMS CC NSE KGE];

% Contigency statistics
    if ~isempty(Z_thr)
      r_h=length(find(Z_tg(k)>Z_thr & Z_rf(k)>Z_thr))/N; % Hit rate
      r_m=length(find(Z_tg(k)<=Z_thr & Z_rf(k)>Z_thr))/N; % Missing rate
      r_f=length(find(Z_tg(k)>Z_thr & Z_rf(k)<=Z_thr))/N; % False alarm rate
      r_n=length(find(Z_tg(k)<=Z_thr & Z_rf(k)<=Z_thr))/N; % Correct negative rate
      CSs.TS(t,:)=[r_h r_m r_f r_n];
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
STs.RS=cat(3,Na,m_tg,m_rf,v_rf,v_tg);

% Error metrics
dZa2(~k)=NaN;
dZa(~k)=NaN;
pZa(~k)=NaN;
RMS=sqrt(dZa2./Na); % Root mean square error
CRMS=sqrt(dZa2./Na-(dZa./Na).^2); % Centered root mean square error
CC=(pZa./Na-m_tg.*m_rf)./sqrt(v_rf.*v_tg); % Correlation coefficient
NSE=1-dZa2./(Na.*v_rf); % Nash Sutcliff efficiency
KGE=1-sqrt((CC-1).^2+(m_tg./m_rf-1).^2+(sqrt(v_tg)./sqrt(v_rf).*m_rf./m_tg-1).^2); % Kling-Gupta efficiency
EMs.RS=cat(3,RMS,CRMS,CC,NSE,KGE);

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
  CSs.RS=cat(3,Na_h,Na_m,Na_f,Na_n); 
end
end
