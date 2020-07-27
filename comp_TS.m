% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 7/18/2020

%% Functionality
% This code performs error analysis for one or multiple target time series with
%  respect to a common reference. Specifically, it calculates statistics (sample
%  size, mean, and variance), error metrics (RMS, CRMS, CC, NSE, and KGE), results
%  of three statistical significant tests (for M(R)E, (N)CRMS, and CC), and contigency
%  statistics (optional, percentage of H, M, F, and N).

%% Input
% ofl : full name list of .mat files store the time series class (TSCls.m) for
%        different locations (ofl can be N-by-1 or N-by-M cell array for 1 or
%        M target products and N stations);
% nflg: a flag iindicating whether to calculate network-based statistics and
%        error metrics or not;

% P_N : minimum percentage of sample size for the statistics and error metrics
%        calculations;
% Tmk : a 2-element cell to specify the time range of interest (the first element
%        can be 'A', 'Y', or 'M' stands for no mask, mask based on year, or mask
%        based on month; the second element can be NaN for 'A' and vector stores
%        the year(s) or month(s) for 'Y' and 'M'; eg., {'M',[1:4 11 12]});
% Sval: a singular value if a time step that both the target and reference values
%        equal to it is excluded from the analysis (default is NaN);
% Thr : thresholds used to calculate the contigency statistics for the time series
%        of different locations (default is NaN);
%  a  : significant level for the statistical significant test (defaul is 0.05);
% pflg: parallel flag (false/true - squential/parallel, default is false).

%% Output
% STs: table stores the statistics of target(s) and reference (sample size, mean,
%       and variance);
% EMs: table stores the error metrics (RMS, CRMS, CC, NSE, and KGE);
% SGs: table stores the significant tests results (ME/MRE, CRMS/NCRMS, and CC);
% Stg: matched time series;

% CSs: table stores the contigency statistics (percentage of H, M, F, and N).

%% Additional note
% If M products are inputted by ofl following the order of product P1, P2, ...,
%  PM, then the outputted results (STs, EMs, SGs, and CSs) follows the reversed
%  order as PM, PM-1, ..., P1.

% (N)(C)RMS - (normalized) (centered) root mean squared error;
% CC - correlation coefficient; NSE - Nash-Sutcliff efficiency;
% KGE - Kling-Gupta efficiency; M(R)E - mean (relative) error;
% H - hit; M - missing; F - false alarm; N - correct negative.

% ME is calculated by subtracting mean of target from mean of reference;
% MRE is calculated by dividing ME to mean of reference.

% Require TSCls.m and errM.m.

function [STs,EMs,SGs,Stg,CSs]=comp_TS(ofl,nflg,varargin)
%% Check the inputs
narginchk(2,7);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'ofl',@(x) validateattributes(x,{'cell'},{'nonempty'},mfilename,'ofl'));
addRequired(ips,'nflg',@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'nflg'));

addOptional(ips,'P_N',0,@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'P_N'));
addOptional(ips,'Tmk',{'A',NaN},@(x) validateattributes(x,{'cell'},{'numel',2},mfilename,'Tmk'));
addOptional(ips,'Sval',NaN,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Sval'));
addOptional(ips,'Thr',NaN,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Thr'));
addOptional(ips,'a',.05,@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'a'));
addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));

parse(ips,ofl,nflg,varargin{:});
P_N=ips.Results.P_N;
Tmk=ips.Results.Tmk;
Sval=ips.Results.Sval;
Thr=ips.Results.Thr;
a=ips.Results.a;
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
      [sts,ems,css,sgs,stg,rn,gid]=comp_TS_sub(ofl(n,:),a,Thr,Tmk,P_N,Sval); % Each station

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
      [sts,ems,css,sgs,stg,rn,gid]=comp_TS_sub(ofl(n,:),a,Thr,Tmk,P_N,Sval); % Each station

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
switch nflg
  case true
    if size(Stg,1)>1 && size(ofl,1)>1
      [sts,ems,sgs,css]=errM(Stg,a,Thr);
      if isa(Gid,'double')
        RN=[RN;sprintf('Group_%i',Gid)];
      elseif iscell(Gid) && ischar(Gid{1})
        RN=[RN;sprintf('G-%s',Gid{1})];
      else
        error('Group ID must be supplied as a scalar for integer or as a cell for character');
      end
      STs=[STs;sts];
      EMs=[EMs;ems];
      SGs=[SGs;sgs];
      CSs=[CSs;css];
    end

  case false
    if isa(Gid,'double')
      fprintf('Skipping statistics and error metrics for Group_%i network\n',Gid);
    elseif iscell(Gid) && ischar(Gid{1})
      fprintf('Skipping statistics and error metrics for G-%s network\n',Gid{1});
    end
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
  if size(ofl,2)>1
    K=num2cell(nchoosek(1:size(ofl,2),2))';
    nl4=cellfun(@(X,Y) sprintf('ME%i_%i',X,Y),K(1,:),K(2,:),'UniformOutput',false);
    nl5=cellfun(@(X,Y) sprintf('CRMS%i_%i',X,Y),K(1,:),K(2,:),'UniformOutput',false);
    nl6=cellfun(@(X,Y) sprintf('CC%i_%i',X,Y),K(1,:),K(2,:),'UniformOutput',false);
    SGs=array2table(SGs,'VariableNames',[nl1 nl2 nl3 nl4 nl5 nl6],'RowNames',RN);
  else
    SGs=array2table(SGs,'VariableNames',[nl1 nl2 nl3],'RowNames',RN);
  end

  if ~isempty(CSs)
    nl1=cellfun(@(X) sprintf('r_h%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
    nl2=cellfun(@(X) sprintf('r_m%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
    nl3=cellfun(@(X) sprintf('r_f%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
    nl4=cellfun(@(X) sprintf('r_n%i',X),num2cell(1:size(ofl,2)),'UniformOutput',false);
    CSs=array2table(CSs,'VariableNames',[nl1 nl2 nl3 nl4],'RowNames',RN);
  end
end
end

function [sts,ems,css,sgs,TS,rn,gid]=comp_TS_sub(ofl,a,Thr,Tmk,P_N,Sval)
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
  [Y,M,~]=datevec(table2array(TS(:,1)));
  switch Tmk{1}
    case 'Y'
      k0=any(Y==Tmk{2},2);
    case 'M'
      k0=any(M==Tmk{2},2);
    case 'A'
      k0=true(size(Y));
  end
  TS=table2array(TS(:,2:end));
  k=all(~isnan(TS),2) & any(TS~=Sval,2) & k0;
  if sum(k)/sum(k0)>P_N
    TS=TS(k,:);
    [sts,ems,sgs,css]=errM(TS,a,Thr);
    rn=OTS.Gtg.Sid;

  else
    TS=[];
    sts=[];
    ems=[];
    css=[];
    sgs=[];
    rn='';
  end

else
  TS=[];
  sts=[];
  ems=[];
  css=[];
  sgs=[];
  rn='';
end
end
