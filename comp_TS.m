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
    parfor n=1:size(ofl,1)
      [sts,ems,css,sgs,stg,rn,gid]=comp_TS_sub(ofl(n,:),a,Thr); % Each station

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
    for n=1:size(ofl,1)
      [sts,ems,css,sgs,stg,rn,gid]=comp_TS_sub(ofl(n,:),a,Thr); % Each station

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
  nl1=cellfun(@(X) sprintf('m_tg%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
  nl2=cellfun(@(X) sprintf('v_tg%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
  STs=array2table(STs,'VariableNames',['N' nl1 'm_rf' nl2 'v_rf'],'RowNames',RN);

  nl1=cellfun(@(X) sprintf('RMS%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
  nl2=cellfun(@(X) sprintf('CRMS%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
  nl3=cellfun(@(X) sprintf('CC%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
  nl4=cellfun(@(X) sprintf('NSE%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
  nl5=cellfun(@(X) sprintf('KGE%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
  EMs=array2table(EMs,'VariableNames',[nl1 nl2 nl3 nl4 nl5],'RowNames',RN);

  nl1=cellfun(@(X) sprintf('ME%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
  nl2=cellfun(@(X) sprintf('CRMS%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
  nl3=cellfun(@(X) sprintf('CC%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
  K=num2cell(nchoosek(1:size(ofl,2),2))';
  nl4=cellfun(@(X,Y) sprintf('ME%i_%i',X,Y),K(1,:),K(2,:),'UniformOutput',false);
  nl5=cellfun(@(X,Y) sprintf('CRMS%i_%i',X,Y),K(1,:),K(2,:),'UniformOutput',false);
  nl6=cellfun(@(X,Y) sprintf('CC%i_%i',X,Y),K(1,:),K(2,:),'UniformOutput',false);
  SGs=array2table(SGs,'VariableNames',[nl1 nl2 nl3 nl4 nl5 nl6],'RowNames',RN);

  if ~isempty(CSs)
    nl1=cellfun(@(X) sprintf('r_h%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
    nl2=cellfun(@(X) sprintf('r_m%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
    nl3=cellfun(@(X) sprintf('r_f%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
    nl4=cellfun(@(X) sprintf('r_n%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
    CSs=array2table(CSs,'VariableNames',[nl1 nl2 nl3 nl4],'RowNames',RN);
  end
end
end

function [sts,ems,css,sgs,TS,rn,gid]=comp_TS_sub(ofl,a,Thr)
%% Align the time series
for i=1:length(ofl)
% Load the TSCls.m object
  OTS=matfile(ofl{i});
  OTS=OTS.OTS;

% Match target time resolution to reference
  [Ttg,Trf]=OTS.UniTL('rf'); % Unify the time zones
  if i==1 % Reference time series
    TS=array2table([Trf OTS.TS2],'VariableNames',{'Dnum','Srf'});
  end

% Aggregate the target time series
  stg=[];
  for t=1:length(Trf)
    k=Ttg>=Trf(t)-OTS.TR2/24/2 & Ttg<Trf(t)+OTS.TR2/24/2;
    if any(k)
      stg=[stg;[Trf(t) mean(OTS.TS1(k))]];
    end
  end

% Join the target with reference
  if ~isempty(stg)
    stg=array2table(stg,'VariableNames',[{'Dnum'} sprintf('Stg%i',i)]);
    TS=innerjoin(stg,TS,'Keys','Dnum');
  else
    TS=[];
    break;
  end
end
clear stg k Trf Ttg

%% Perform the stats/error calculation
gid=OTS.gid;
if size(TS,1)>4
  TS=table2array(TS(:,2:end));
  [sts,ems,sgs,css]=errM(TS,a,Thr);
  rn=OTS.Gtg.ID;

else
  TS=[];
  sts=[];
  ems=[];
  css=[];
  sgs=[];
  rn='';
end
end
