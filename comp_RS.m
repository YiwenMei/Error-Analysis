% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 7/10/2020

%% Functionality
% This code performs error analysis between target and reference image stacks.
%  The image stacks are supplied as space-time class (V2DTCls.m) or space-time-profile
%  class (V3DTCls.m) object.

% The code calculates statistics (sample size, mean, and variance), error metrics
%  (RMS, CRMS, CC, NSE, and KGE), results of three statistical significant tests
%  (for M(R)E, (N)CRMS, and CC), and contigency statistics (optional, percentage
%  of H, M, F, and N).

%% Input
% OStg : V2DTCls.m object for the target image stack;
% OSrf : V2DTCls.m object for the reference image stack;
% wkpth: working directory of the code;

% OSmk: V2DCls.m object for the binary (1 - basin, 0 - not basin) basin mask
%        (default is NaN);
% Tmk : a 2-element cell to specify the time range of interest (the first element
%        can be 'A', 'Y', or 'M' stands for no mask, mask based on year, or mask
%        based on month; the second element can be NaN for 'A' and vector stores
%        the year(s) or month(s) for 'Y' and 'M'; eg., {'M',[1:4 11 12]});
%  L  : a 1-by-2 vector of measurement depth/height for OStg and OSrf if they
%        are V3DTCls (default is [NaN NaN] for V2DCls);
% pflg: parallel flag (false/true - squential/parallel, default is false);
% P_N : minimum percentage of sample size for the statistics and error metrics
%        calculations;
% Sval: a singular value if a time step that both the target and reference values
%        equal to it is excluded from the analysis (default is NaN);
% Thr : threshold used to calculate the contigency statistics for the time series
%        of different locations (default is NaN);
%  cf : a 2-by-2 matrix to specify the multiplicative (column 1) and additive
%        (column 2) factors for target (row 1) and reference (row 2) data to
%        a user-specified unit (default is [1 0;1 0]);
%  a  : significant level for the statistical significant test (defaul is 0.05);

%% Output
% TS/ structure array stores the time series/spatial maps of statistics, error
% RS:  metrics, significant tests results, time/spatial mask, and contigency
%      statistics (optional); the strucure array has four fields, STs, EMs, SGs,
%      and CSs, that corresponds to the four results; it also contains the time/spatial
%      mask, Tln/Smk, for TS/RS.

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

% Tln includes two columns, column 1 for time line as date number and column
%  2 for time mask with possible values of 0, 1, and 2. 0 represents time steps
%  that calculations executed; 1 or 2 represents disgarded time step (1 - outside
%  of time range of interest and 2 - less than 50% of data points for calculations);
% Smk includes three maps for X coordinates, Y coordinates, and spatial mask
%  with possible values of 0, 1, and 2. 0 represents grids that calculations
%  executed; 1 or 2 represents disgarded grids (1 - outside of basin and land
%  mask and 2 - less than 50% of data points for calculations).

% Require V2DTCls.m, V3DTCls.m, errM.m, and V2DCls (required if OSmk is supplied).

function [TS,RS]=comp_RS(OStg,OSrf,wkpth,varargin)
%% Check the inputs
narginchk(3,11);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'OStg',@(x) validateattributes(x,{'V2DTCls','V3DTCls','Wind2DTCls'},{'nonempty'},...
    mfilename,'OStg'));
addRequired(ips,'OSrf',@(x) validateattributes(x,{'V2DTCls','V3DTCls','Wind2DTCls'},{'nonempty'},...
    mfilename,'OSrf'));
addRequired(ips,'wkpth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wkpth'));

addOptional(ips,'OSmk',NaN,@(x) validateattributes(x,{'V2DCls'},{'nonempty'},mfilename,'OSmk'));
addOptional(ips,'Tmk',{'A',NaN},@(x) validateattributes(x,{'cell'},{'numel',2},mfilename,'Tmk'));
addOptional(ips,'L',[NaN NaN],@(x) validateattributes(x,{'double'},{'vector'},mfilename,'L'));
addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));
addOptional(ips,'P_N',0,@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'P_N'));
addOptional(ips,'Sval',NaN,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Sval'));
addOptional(ips,'Thr',NaN,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Thr'));
addOptional(ips,'cf',[1 0;1 0],@(x) validateattributes(x,{'double'},{'size',[2 2]},mfilename,'cf'));
addOptional(ips,'a',.05,@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'a'));

parse(ips,OStg,OSrf,wkpth,varargin{:});
OSmk=ips.Results.OSmk;
Tmk=ips.Results.Tmk;
L=ips.Results.L;
pflg=ips.Results.pflg;
P_N=ips.Results.P_N;
Sval=ips.Results.Sval;
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
[~,~,~,rsn]=OStg.GridCls;
[xcn,ycn,~,SRcn]=OSrf.GridCls; % Reference grid as canonical grid if all(rs)<1
rs=rsn./SRcn;

if all(rs>=1) % Canonical grid is new (corners and resolution of xcn/ycn changed)
  SRcn=SRcn.*rs;
  xcn=imresize(xcn,1/mean(rs),'bilinear');
  ycn=imresize(ycn,1/mean(rs),'bilinear');
  
  sz=abs(diff(xcn(1,:),[],2));
  sz=find(abs(sz([1 end])-SRcn(1))>.5*SRcn(1));
  if ~isempty(sz)
    sz(sz==2)=size(xcn,2);
    xcn(:,sz)=[];
    ycn(:,sz)=[];
  end

  sz=abs(diff(ycn(:,1),[],1));
  sz=find(abs(sz([1 end])-SRcn(2))>.5*SRcn(2));
  if ~isempty(sz)
    sz(sz==2)=size(ycn,1);
    ycn(sz,:)=[];
    xcn(sz,:)=[];
  end
  fprintf('New reference grid created\n');

elseif any(rs<1) && any(rs>=1)
  error('Target and reference x-y resolution must be finer/coarser than each other at the same time');
end

%% Statistics and error metrics
% Initialization
fprintf('Generate the ID files\n');
Srf=Taggr(OSrf,Tcn,TRcn,1,L(2));
Stg=Taggr(OStg,Tcn,TRcn,1,L(1));
Srf=cf(2,1)*Saggr(OSrf,Srf,xcn,ycn,SRcn,fullfile(wkpth,'id_rf.mat'))+cf(2,2);
Stg=cf(1,1)*Saggr(OStg,Stg,xcn,ycn,SRcn,fullfile(wkpth,'id_tg.mat'))+cf(1,2);

if isa(OSmk,'V2DCls')
  [~,~,sz]=OSmk.GridCls;
  if all(size(Srf)==sz)
    mk=OSmk.readCls;
  else
    error('Dimensions of basin mask must be the same as reference grid');
  end
else
  mk=true(size(Srf));
end
N=zeros(size(Srf));
Z_tg=zeros(size(Srf));
Z_rf=zeros(size(Srf));
Z2_tg=zeros(size(Srf));
Z2_rf=zeros(size(Srf));
pZ=zeros(size(Srf));
dZ=zeros(size(Srf));
dZ2=zeros(size(Srf));
Nh=zeros(size(Srf));
Nm=zeros(size(Srf));
Nf=zeros(size(Srf));
Nn=zeros(size(Srf));
STs=[];
EMs=[];
SGs=[];
CSs=[];
Tln=zeros(size(Tcn));

% Main loop
fprintf('Execute the error metrics calculation\n');
switch pflg
  case true
    parfor t=1:length(Tcn)
      if ~TimeSkip(Tcn(t),Tmk)
% Map to the canonical time line and grids
        Srf=Taggr(OSrf,Tcn,TRcn,t,L(2));
        Stg=Taggr(OStg,Tcn,TRcn,t,L(1));
        Srf=cf(2,1)*Saggr(OSrf,Srf,xcn,ycn,SRcn,fullfile(wkpth,'id_rf.mat'))+cf(2,2);
        Stg=cf(1,1)*Saggr(OStg,Stg,xcn,ycn,SRcn,fullfile(wkpth,'id_tg.mat'))+cf(1,2);

% Temporal statistics and error metrics - preparation
        k=~isnan(Srf) & ~isnan(Stg) & (Srf~=Sval | Stg~=Sval) & mk; % Common grids & ~singular value
        Srf(~k)=0;                                                  %  & in basin
        Stg(~k)=0;

        N=N+k;
        Z_rf=Z_rf+Srf;
        Z_tg=Z_tg+Stg;
        Z2_rf=Z2_rf+Srf.^2;
        Z2_tg=Z2_tg+Stg.^2;
        pZ=pZ+Srf.*Stg;
        dZ=dZ+(Stg-Srf);
        dZ2=dZ2+(Stg-Srf).^2;
        if ~isnan(Thr)
          Nh=Nh+(Stg>Thr & Srf>Thr);
          Nm=Nm+(Stg<=Thr & Srf>Thr);
          Nf=Nf+(Stg>Thr & Srf<=Thr);
          Nn=Nn+(Stg<=Thr & Srf<=Thr);
        end

% Spatial statistics and error metrics
        if length(find(k))/length(find(mk))>P_N % Remaining grids over grids in basin
          [sts,ems,sgs,css]=errM([Stg(k) Srf(k)],a,Thr);
          STs=[STs;sts]; % Statistics
          EMs=[EMs;ems]; % Error metrics
          SGs=[SGs;sgs]; % Error metrics
          if ~isnan(Thr)
            CSs=[CSs;css]; % Contigency statistics
          end

        else
          Tln(t)=2; % Label for not enough data point after applying the masks
        end
      else
        Tln(t)=1; % Label for not within the time mask
      end
    end

  case false
    for t=1:length(Tcn)
      if ~TimeSkip(Tcn(t),Tmk)
% Map to the canonical time line and grids
        Srf=Taggr(OSrf,Tcn,TRcn,t,L(2));
        Stg=Taggr(OStg,Tcn,TRcn,t,L(1));

        Srf=cf(2,1)*Saggr(OSrf,Srf,xcn,ycn,SRcn,fullfile(wkpth,'id_rf.mat'))+cf(2,2);
        Stg=cf(1,1)*Saggr(OStg,Stg,xcn,ycn,SRcn,fullfile(wkpth,'id_tg.mat'))+cf(1,2);

        Cin=cat(3,N,Z_tg,Z_rf,Z2_tg,Z2_rf,pZ,dZ,dZ2,Nh,Nm,Nf,Nn);
        [N,Z_tg,Z_rf,Z2_tg,Z2_rf,pZ,dZ,dZ2,Nh,Nm,Nf,Nn,STs,EMs,SGs,CSs,tln]=...
            comp_RS_sub(Srf,Stg,mk,P_N,Sval,Thr,a,Cin,STs,EMs,SGs,CSs);

        Tln(t)=tln;
      else
        Tln(t)=1;
      end
    end
    clear Cin
end
delete(fullfile(wkpth,'id_*.mat'));
fprintf('Finishing the loop\n');

%% Temporal statistics and error metrics - calculation
k=N./sum(Tln~=1)>P_N;
mk=double(~mk); % Label for not within the time mask
mk(~k & mk==0)=2; % Label for not enough data point after applying the masks
mk=cat(3,xcn,ycn,mk);

N(~k)=NaN;

% Statistics
m_tg=Z_tg./N; % Mean of target time series
m_rf=Z_rf./N; % Mean of reference time series
v_tg=Z2_tg./(N-1)-m_tg.^2; % Variance of target time series
v_rf=Z2_rf./(N-1)-m_rf.^2; % Variance of reference time series
k=v_tg<0;
if any(k)
  fprintf('%d out of %d values are negative in the target variance\n',length(find(k)),...
      length(find(~isnan(v_tg))));
end
k=v_rf<0;
if any(k)
  fprintf('%d out of %d values are negative in the reference variance\n',length(find(k)),...
      length(find(~isnan(v_rf))));
end

k=v_tg<0 | v_rf<0;
m_tg(k)=NaN;
m_rf(k)=NaN;
v_tg(k)=NaN;
v_rf(k)=NaN;
clear Z_tg Z_rf Z2_tg Z2_rf k

% Error metrics
RMS=sqrt(dZ2./N); % Root mean square error
CRMS=sqrt(dZ2./(N-1)-(dZ./N).^2); % Centered root mean square error
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
if ~isnan(Thr)
  Nh=Nh./N; % Hit rate
  Nm=Nm./N; % Missing rate
  Nf=Nf./N; % False alarm rate
  Nn=Nn./N; % Correct negative rate

% Resutls
  RS=struct('STs',cat(3,N,m_tg,m_rf,v_tg,v_rf),'EMs',cat(3,RMS,CRMS,CC,NSE,KGE),...
      'SGs',cat(3,H1,H2,H3),'CSs',cat(3,Nh,Nm,Nf,Nn),'Smk',mk);
  TS=struct('STs',STs,'EMs',EMs,'SGs',SGs,'CSs',CSs,'Tln',[Tcn Tln]);
else
  RS=struct('STs',cat(3,N,m_tg,m_rf,v_tg,v_rf),'EMs',cat(3,RMS,CRMS,CC,NSE,KGE),...
      'SGs',cat(3,H1,H2,H3),'Smk',mk);
  TS=struct('STs',STs,'EMs',EMs,'SGs',SGs,'Tln',[Tcn Tln]);
end
fprintf('Finishing the error metrics calculation\n');
end

function S=Taggr(obj,Tcn,TRcn,t,l)
% Find overlapped time steps
TL=obj.TimeCls('center'); % Centered time line
ID=find(TL>=Tcn(t) & TL<Tcn(t)+TRcn/24);

% Aggregate the time steps
S=[];
if length(ID)>=ceil(.5*TRcn/obj.TmR)
  N=[];
  for i=1:length(ID)
    if isa(obj,'V3DTCls')
      Z=obj.readCls(ID(i),l);
    else
      Z=obj.readCls(ID(i));
    end
    S=nansum(cat(3,S,Z),3);
    N=nansum(cat(3,N,~isnan(Z)),3);
  end
  S=S./N;

else
  S=nan(size(obj.GridCls));
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
  [~,~,~,rsn]=obj.GridCls;
  [id,d]=knnsearch([x y],[xcn ycn],'K',ceil(prod(SRcn./rsn)));
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

function [N,Z_tg,Z_rf,Z2_tg,Z2_rf,pZ,dZ,dZ2,Nh,Nm,Nf,Nn,STs,EMs,SGs,CSs,tln]=...
    comp_RS_sub(Srf,Stg,mk,P_N,Sval,Thr,a,Cin,STs,EMs,SGs,CSs)
% Temporal statistics and error metrics - preparation
k=~isnan(Srf) & ~isnan(Stg) & (Srf~=Sval | Stg~=Sval) & mk; % Common grids & ~singular value & in basin
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
if ~isnan(Thr)
  Nh=Nh+(Stg>Thr & Srf>Thr);
  Nm=Nm+(Stg<=Thr & Srf>Thr);
  Nf=Nf+(Stg>Thr & Srf<=Thr);
  Nn=Nn+(Stg<=Thr & Srf<=Thr);
end

% Spatial statistics and error metrics
if length(find(k))/length(find(mk))>P_N % Remaining grids over grids in basin
  [sts,ems,sgs,css]=errM([Stg(k) Srf(k)],a,Thr);
  STs=[STs;sts]; % Statistics
  EMs=[EMs;ems]; % Error metrics
  SGs=[SGs;sgs]; % Significant tests
  if ~isnan(Thr)
    CSs=[CSs;css]; % Contigency statistics
  end

  tln=0;
else
  tln=2;
end
end

function sk=TimeSkip(tcn,Tmk)
[Y,M,~]=datevec(tcn);
switch Tmk{1}
  case 'Y'
    sk=~any(Y==Tmk{2});
  case 'M'
    sk=~any(M==Tmk{2});
  case 'A'
    sk=false;
end
end
