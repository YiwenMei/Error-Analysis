% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 12/11/2019

%% Functionality
% This code performs error analysis between target and reference image stacks.
%  The image stacks are supplied as space-time class (V2DTCls.m) object.

% The code calculates statistics (sample size, mean, and variance), error metrics
%  (RMS, CRMS, CC, NSE, and KGE), results of three statistical significant tests
%  (for M(R)E, (N)CRMS, and CC), and contigency statistics (optional, percentage
%  of H, M, F, and N).

%% Input
% OStg : V2DTCls.m object for the target image stack;
% OSrf : V2DTCls.m object for the reference image stack;
% wkpth: working directory of the code;

% pflg: parallel flag (false/true - squential/parallel, default is false);
% Thr : thresholds used to calculate the contigency statistics for the time series
%        of different locations (default is []).
%  cf : multiplicative conversion factor to reference unit (default is 1);
%  a  : significant level for the statistical significant test (defaul is 0.05);

%% Output
% TS/RS: structure array stores the time series/spatial image of statistics,
%         error metrics, significant tests, and contigency statistics (optional);
%         the strucure array has four fields, STs, EMs, SGs, and CSs, that corresponds
%         to the four results.

%% Additional note
% STs includes time series/spatial images of
%  1)common sample size,     2)mean of target,        3)mean of reference,
%  4)variance of target, and 5)variance of reference;
% EMs includes time series/spatial images of
%  1)root mean square error,  2)centered root mean square error,
%  3)correlation coefficient, 4)Nash-Sutcliffe efficiency, and
%  5)Kling-Gupta efficiency;
% CSs includes time series/spatial images of
%  rate of 1)hit, 2)missing, 3)false alarm, and 4)correct negative.

function [TS,RS]=comp_RS(OStg,OSrf,wkpth,varargin)
%% Check the inputs
narginchk(3,7);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'OStg',@(x) validateattributes(x,{'V2DTCls','Wind2DTCls'},{'nonempty'},...
    mfilename,'OStg'));
addRequired(ips,'OSrf',@(x) validateattributes(x,{'V2DTCls','Wind2DTCls'},{'nonempty'},...
    mfilename,'OSrf'));
addRequired(ips,'wkpth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wkpth'));

addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},...
    mfilename,'pflg'));
addOptional(ips,'Thr',[],@(x) validateattributes(x,{'double'},{'nonempty'},...
    mfilename,'Thr'));
addOptional(ips,'cf',1,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'cf'));
addOptional(ips,'a',.05,@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'a'));

parse(ips,OStg,OSrf,wkpth,varargin{:});
pflg=ips.Results.pflg;
Thr=ips.Results.Thr;
cf=ips.Results.cf;
a=ips.Results.a;
clear ips varargin

%% Align the images
% Canonical time line with begin time convention
rt=OStg.TmR/OSrf.TmR;
Tcn=OSrf.TimeCls('begin');
TRcn=OSrf.TmR;

if rt>1 % Reference is finer
  TRcn=TRcn*rt;
  Tcn=min(Tcn):TRcn/24:max(Tcn);
end

% Canonical grids
rs=OStg.GIf(3,:)/OSrf.GIf(3,:);
[xcn,ycn]=OSrf.GridCls;
SRcn=OSrf.GIf(3,:);

if all(rs>=1)
  SRcn=SRcn*rs;
  xcn=imresize(xcn,1/rs,'bilinear');
  ycn=imresize(ycn,1/rs,'bilinear');

elseif any(rs<1) && any(rs>=1)
  error('Target and reference x-y resolution must be finer/coarser than each other at the same time');
end

%% Statistics and error metrics
% Initialization
Srf=Taggr(OSrf,Tcn,TRcn,1);
Stg=Taggr(OStg,Tcn,TRcn,1);

Srf=Saggr(OSrf,Srf,xcn,ycn,SRcn,fullfile(wkpth,'id_rf.mat'));
Stg=cf*Saggr(OStg,Stg,xcn,ycn,SRcn,fullfile(wkpth,'id_tg.mat'));

[N,Z_tg,Z_rf,Z2_tg,Z2_rf,pZ,dZ,dZ2,Nh,Nm,Nf,Nn,STs,EMs,SGs,CSs]=...
    comp_RS_sub(Srf,Stg,Thr,a,repmat(zeros(size(xcn)),1,1,12),[],[],[],[]);

% Main loop
switch pflg
  case true
    parfor t=2:length(Tcn)
% Map to the canonical time line and grids
      Srf=Taggr(OSrf,Tcn,TRcn,t);
      Stg=Taggr(OStg,Tcn,TRcn,t);

      Srf=Saggr(OSrf,Srf,xcn,ycn,SRcn,fullfile(wkpth,'id_rf.mat'));
      Stg=cf*Saggr(OStg,Stg,xcn,ycn,SRcn,fullfile(wkpth,'id_tg.mat'));

% Temporal statistics and error metrics - preparation
      k=~isnan(Srf) & ~isnan(Stg); % Common grids
      Srf(~k)=0;
      Stg(~k)=0;

      N=N+k;
      Z_rf=Z_rf+Srf;
      Z_tg=Z_tg+Stg;
      Z2_rf=Z2_rf+Srf.^2;
      Z2_tg=Z2_tg+Stg.^2;
      pZ=pZ+Srf.*Stg;
      dZ=dZ+(Stg-Srf);
      dZ2=dZ2+(Stg-Srf).^2;

      if ~isempty(Thr)
        Nh=Nh+(Stg>Thr & Srf>Thr);
        Nm=Nm+(Stg<=Thr & Srf>Thr);
        Nf=Nf+(Stg>Thr & Srf<=Thr);
        Nn=Nn+(Stg<=Thr & Srf<=Thr);
      end

% Spatial statistics and error metrics
      M=length(find(k)); % Sample size
      if M>4
        [sts,ems,sgs,css]=errM([Stg(k) Srf(k)],a,Thr);
        STs=[STs;sts]; % Statistics
        EMs=[EMs;ems]; % Error metrics
        SGs=[SGs;sgs]; % Error metrics
        if ~isempty(Thr)
          CSs=[CSs;css]; % Contigency statistics
        end
      end
    end

  case false
    for t=2:length(Tcn)
% Map to the canonical time line and grids
      Srf=Taggr(OSrf,Tcn,TRcn,t);
      Stg=Taggr(OStg,Tcn,TRcn,t);

      Srf=Saggr(OSrf,Srf,xcn,ycn,SRcn,fullfile(wkpth,'id_rf.mat'));
      Stg=cf*Saggr(OStg,Stg,xcn,ycn,SRcn,fullfile(wkpth,'id_tg.mat'));

      Cin=cat(3,N,Z_tg,Z_rf,Z2_tg,Z2_rf,pZ,dZ,dZ2,Nh,Nm,Nf,Nn);
      [N,Z_tg,Z_rf,Z2_tg,Z2_rf,pZ,dZ,dZ2,Nh,Nm,Nf,Nn,STs,EMs,SGs,CSs]=...
          comp_RS_sub(Srf,Stg,Thr,a,Cin,STs,EMs,SGs,CSs);
    end
    clear Cin
end
delete(fullfile(wkpth,'id_*.mat'));

%% Temporal statistics and error metrics - calculation
k=N>4;
N(~k)=NaN;

% Statistics
m_tg=Z_tg./N; % Mean of target time series
m_rf=Z_rf./N; % Mean of reference time series
v_tg=Z2_tg./N-m_tg.^2; % Variance of target time series
v_rf=Z2_rf./N-m_rf.^2; % Variance of reference time series
clear Z_tg Z_rf Z2_tg Z2_rf

% Error metrics
RMS=sqrt(dZ2./N); % Root mean square error
CRMS=sqrt(dZ2./N-(dZ./N).^2); % Centered root mean square error
CC=(pZ./N-m_tg.*m_rf)./sqrt(v_rf.*v_tg); % Correlation coefficient
NSE=1-dZ2./(N.*v_rf); % Nash Sutcliff efficiency
% Kling-Gupta efficiency
KGE=1-sqrt((CC-1).^2+(m_tg./m_rf-1).^2+(sqrt(v_tg)./sqrt(v_rf).*m_rf./m_tg-1).^2);
clear dZ2 dZ pZ

% Significant tests
ts=(m_tg-m_rf)./(CRMS./sqrt(N)); % m_tg-m_rf-uthr
p=1-tcdf(abs(ts),N-1)+tcdf(-abs(ts),N-1);
H1=p<a; % M(R)E sig<> 0 (1) or not (0)
ts=(N-1).*CRMS.^2./v_rf; % CRMS.^2./v_rf/uthr
p=chi2cdf(ts,N-1); % left-tailed
H2=p<a; % (N)CRMS sig< SD_rf (1) or not (0);
ts=.5*log((1+CC)./(1-CC)).*sqrt(N-3); % .5*log((1+CC)./(1-CC))-uthr
p=1-normcdf(ts,0,1); % right-tailed
H3=p<a; % CC sig> 0 (1) or not (0)

% Contigency statistics
if ~isempty(Thr)
  Nh=Nh./N; % Hit rate
  Nm=Nm./N; % Missing rate
  Nf=Nf./N; % False alarm rate
  Nn=Nn./N; % Correct negative rate

% Resutls
  RS=struct('STs',cat(3,N,m_tg,m_rf,v_tg,v_rf),'EMs',cat(3,RMS,CRMS,CC,NSE,KGE),...
      'SGs',cat(3,H1,H2,H3),'CSs',cat(3,Nh,Nm,Nf,Nn));
  TS=struct('STs',STs,'EMs',EMs,'SGs',SGs,'CSs',CSs);
else
  RS=struct('STs',cat(3,N,m_tg,m_rf,v_tg,v_rf),'EMs',cat(3,RMS,CRMS,CC,NSE,KGE),...
      'SGs',cat(3,H1,H2,H3));
  TS=struct('STs',STs,'EMs',EMs,'SGs',SGs);
end
end

function S=Taggr(obj,Tcn,TRcn,t)
% Find overlapped time steps
TL=obj.TimeCls('center'); % Centered time line
ID=find(TL>=Tcn(t) & TL<Tcn(t)+TRcn/24);

% Aggregate the time steps
S=[];
if length(ID)>=ceil(.5*TRcn/obj.TmR)
  N=[];
  for i=1:length(ID)
    Z=obj.readCls(ID(i));
    S=nansum(cat(3,S,Z),3);
    N=nansum(cat(3,N,~isnan(Z)),3);
  end
  S=S./N;

else
  S=nan(size(obj.readCls(t)));
end
end

function S=Saggr(obj,Z,xcn,ycn,SRcn,idn)
s=size(xcn);
thr=hypot(SRcn(1)/2,SRcn(2)/2);

if exist(idn,'file')~=2
% Original grids
  [x,y]=obj.GridCls;
  x=reshape(x,numel(x),1);
  y=reshape(y,numel(y),1);

% Canonical grids
  xcn=reshape(xcn,numel(xcn),1);
  ycn=reshape(ycn,numel(ycn),1);

% Map to canoical grids
  [id,d]=knnsearch([x y],[xcn ycn],'K',ceil(prod(SRcn./obj.GIf(3,:))));
  clear x y xcn ycn
  save(idn,'d','id');
else
  load(idn,'d','id');
end

% Aggregation
Z=reshape(Z,numel(Z),1);
S=Z(id);
clear id
S(d>thr)=NaN;
d=sum(~isnan(S),2)/size(S,2); % recycle d
S=nanmean(S,2);
S(d<.5)=NaN;
S=reshape(S,s);
end

function [N,Z_tg,Z_rf,Z2_tg,Z2_rf,pZ,dZ,dZ2,Nh,Nm,Nf,Nn,STs,EMs,SGs,CSs]=...
    comp_RS_sub(Srf,Stg,Thr,a,Cin,STs,EMs,SGs,CSs)
% Temporal statistics and error metrics - preparation
k=~isnan(Srf) & ~isnan(Stg); % Common grids
Srf(~k)=0;
Stg(~k)=0;

N=Cin(:,:,1)+k;
Z_tg=Cin(:,:,2)+Stg;
Z_rf=Cin(:,:,3)+Srf;
Z2_tg=Cin(:,:,4)+Stg.^2;
Z2_rf=Cin(:,:,5)+Srf.^2;
pZ=Cin(:,:,6)+Srf.*Stg;
dZ=Cin(:,:,7)+Stg-Srf;
dZ2=Cin(:,:,8)+(Stg-Srf).^2;

Nh=Cin(:,:,9);
Nm=Cin(:,:,10);
Nf=Cin(:,:,11);
Nn=Cin(:,:,12);
if ~isempty(Thr)
  Nh=Nh+(Stg>Thr & Srf>Thr);
  Nm=Nm+(Stg<=Thr & Srf>Thr);
  Nf=Nf+(Stg>Thr & Srf<=Thr);
  Nn=Nn+(Stg<=Thr & Srf<=Thr);
end

% Spatial statistics and error metrics
M=length(find(k)); % Sample size
if M>4
  [sts,ems,sgs,css]=errM([Stg(k) Srf(k)],a,Thr);
  STs=[STs;sts]; % Statistics
  EMs=[EMs;ems]; % Error metrics
  SGs=[SGs;sgs]; % Significant tests
  if ~isempty(Thr)
    CSs=[CSs;css]; % Contigency statistics
  end
end
end
