% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 11/20/2019

%% Functionality
% This code performs error analysis for one or two target time series with respect
%  to a reference. Specifically, it calculates statistics (sample size, mean,
%  and variance), error metrics (RMS, CRMS, CC, NSE, and KGE), results of three
%  statistical significant tests (for ME/MRE, CRMS/NCRMS, and CC), and contigency
%  statistics (optional, percentage of H, M, F, and N).

%% Input
% ofl: full name list of .mat files store the time series class (TSCls.m) for
%       different locations (ofl can be N-by-1 or N-by-2 cell array for 1 or
%       2 target and N stations);

%  a  : significant level for the statistical significant test (defaul is 0.05);
% Thr : thresholds used to calculate the contigency statistics for the time series
%        of different locations (default is []);
% pflg: parallel flag (false/true - squential/parallel, default is false).

%% Output
% STs: table stores the statistics of target(s) and reference (sample size, mean,
%       and variance);
% EMs: table stores the error metrics (RMS, CRMS, CC, NSE, and KGE);
% SGs: table stores the significant tests results (ME/MRE, CRMS/NCRMS, and CC);

% CSs: table stores the contigency statistics (percentage of H, M, F, and N).

%% Additional note
% (N)(C)RMS - (normalized) (centered) root mean squared error;
% CC - correlation coefficient; NSE - Nash-Sutcliff efficiency;
% KGE - Kling-Gupta efficiency; M(R)E - mean (relative) error;
% H - hit; M - missing; F - false alarm; N - correct negative.

% ME is calculated by subtracting mean of target from mean of reference;
% MRE is calculated by dividing ME to mean of reference.

% Require TSCLs.m and errM.m.

function [STs,EMs,SGs,CSs]=comp_TS(ofl,varargin)
%% Check the inputs
narginchk(1,4);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'ofl',@(x) validateattributes(x,{'cell'},{'nonempty'},mfilename,'ofl'));

addOptional(ips,'a',.05,@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,...
    'a'));
addOptional(ips,'Thr',[],@(x) validateattributes(x,{'double'},{'nonempty'},...
    mfilename,'Thr'));
addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},...
    mfilename,'pflg'));

parse(ips,ofl,varargin{:});
a=ips.Results.a;
Thr=ips.Results.Thr;
pflg=ips.Results.pflg;
clear ips varargin

%% Statistics and Error metrics
Stg=[];
Gid=[];
RN={};
STs=[];
EMs=[];
SGs=[];
CSs=[];

% Each station
switch pflg
  case true
    parfor n=1:length(ofl)
      [sts,ems,css,sgs,stg,rn,gid]=comp_TS_sub(ofl,n,a,Thr); % Each station

      STs=[STs;sts];
      EMs=[EMs;ems];
      SGs=[SGs;sgs];
      CSs=[CSs;css];

      Stg=[Stg;stg];
      RN=[RN;{rn}];
      if n==length(ofl)
        Gid=[Gid;gid];
      end
    end

  case false
    for n=1:length(ofl)
      [sts,ems,css,sgs,stg,rn,gid]=comp_TS_sub(ofl,n,a,Thr); % Each station

      STs=[STs;sts];
      EMs=[EMs;ems];
      SGs=[SGs;sgs];
      CSs=[CSs;css];

      Stg=[Stg;stg];
      RN=[RN;rn];
      if n==length(ofl)
        Gid=[Gid;gid];
      end
    end
end

% All stations
if size(Stg,1)>1
  [sts,ems,sgs,css]=errM(Stg,a,Thr);
  RN=[RN;sprintf('Group_%i',Gid)];
  STs=[STs;sts];
  EMs=[EMs;ems];
  SGs=[SGs;sgs];
  CSs=[CSs;css];
end

%% Label the results
if ~isempty(RN)
  if size(ofl,2)==1
    STs=array2table(STs,'VariableNames',{'N','m_tg','m_rf','v_tg','v_rf'},'RowNames',RN);
    EMs=array2table(EMs,'VariableNames',{'RMS','CRMS','CC','NSE','KGE'},'RowNames',RN);
    if ~isempty(CSs)
      CSs=array2table(CSs,'VariableNames',{'r_h','r_m','r_f','r_n'},'RowNames',RN);
    end
    
  elseif size(ofl,2)==2
    STs=array2table(STs,'VariableNames',{'N','m_tg1','m_tg2','m_rf','v_tg1',...
        'v_tg2','v_rf'},'RowNames',RN);
    EMs=array2table(EMs,'VariableNames',{'RMS1','RMS2','CRMS1','CRMS2','CC1',...
        'CC2','NSE1','NSE2','KGE1','KGE2'},'RowNames',RN);
    if ~isempty(CSs)
      CSs=array2table(CSs,'VariableNames',{'r_h1','r_h2','r_m1','r_m2','r_f1',...
        'r_f2','r_n1','r_n2'},'RowNames',RN);
    end
  end

  SGs=array2table(SGs,'VariableNames',{'ME','CRMS','CC'},'RowNames',RN);
end
end

function [sts,ems,css,sgs,stg,rn,gid]=comp_TS_sub(ofl,n,a,Thr)
%% One target
if size(ofl,2)==1
% Load the TSCls.m object
  OTS=matfile(ofl{n});
  OTS=OTS.OTS;

% Match target time resolution to reference
  [Ttg,Trf]=OTS.UniTL('rf'); % Unify the time zones

% Aggregate the target time series
  stg=[];
  for t=1:length(Trf)
    k=Ttg>=Trf(t)-OTS.TR2/24/2 & Ttg<Trf(t)+OTS.TR2/24/2;
    if any(k)
      stg=[stg;[mean(OTS.TS1(k)) OTS.TS2(t)]];
    end
  end

elseif size(ofl,2)==2
%% Two target
% Load the TSCls.m object
  OTS=matfile(ofl{n,1});
  OTS=OTS.OTS;
  OTS2=matfile(ofl{n,2});
  OTS2=OTS2.OTS;

% Match target time resolution to reference
  [Ttg1,Trf]=OTS.UniTL('rf'); % Unify the time zones
  [Ttg2,trf]=OTS2.UniTL('rf'); % Unify the time zones
  if any(trf~=Trf)
    error('Reference series time line mismatched');
  end
  clear trf

% Aggregate the target time series
  stg=[];
  for t=1:length(Trf)
    k1=Ttg1>=Trf(t)-OTS.TR2/24/2 & Ttg1<Trf(t)+OTS.TR2/24/2;
    k2=Ttg2>=Trf(t)-OTS2.TR2/24/2 & Ttg2<Trf(t)+OTS2.TR2/24/2;
    if any(k1) && any(k2)
      stg=[stg;[mean(OTS.TS1(k1)) mean(OTS2.TS1(k2)) OTS.TS2(t)]];
    end
  end
end

%% Perform the stats/error calculation
gid=OTS.gid;
if size(stg,1)>1
  [sts,ems,sgs,css]=errM(stg,a,Thr);
  rn=OTS.Gtg.ID;

else
  sts=[];
  ems=[];
  css=[];
  sgs=[];
  rn='';
end
end
