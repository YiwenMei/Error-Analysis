% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/13/2019

%% Functionality
% This code performs error analysis between target and reference time series.

%% Input
% ofl : full name list of .mat files store the time series class (TSCls.m) for
%       different locations;

% Thr : thresholds used to calculate the contigency statistics for the time series
%       of different locations (optional);
% pflg: parallel flag (optional; false, default - squential; true - parallel).

%% Output
% STs: table stores the statistics of target and reference time series;
% EMs: table stores the error metrics for the target time series;

% CSs: table stores the contigency statistics for the target time series (optional).

%% Additional note
% Values within STs, EMs, and CSs share the same meaning with those of errM.m.
% Require TSCLs.m and errM.m.

function [STs,EMs,CSs]=comp_TS(ofl,varargin)
%% Check the inputs
narginchk(1,3);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'ofl',@(x) validateattributes(x,{'cell'},{'nonempty'},mfilename,'ofl'));

addOptional(ips,'Thr',[],@(x) validateattributes(x,{'double'},{},mfilename,'Thr'));
addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));

parse(ips,ofl,varargin{:});
Thr=ips.Results.Thr;
pflg=ips.Results.pflg;
clear ips varargin

%% Statistics and Error metrics
Stg=[];
RN={};
STs=[];
CSs=[];
EMs=[];
Gid=[];

% Each station
switch pflg
  case true
    parfor n=1:length(ofl)
      [sts,ems,css,stg,rn,gid]=comp_TS_sub(ofl,n,Thr); % Each station

      CSs=[CSs;css];
      STs=[STs;sts];
      EMs=[EMs;ems];

      Stg=[Stg;stg];
      RN=[RN;{rn}];
      if n==length(ofl)
        Gid=[Gid;gid];
      end
    end

  case false
    for n=1:length(ofl)
      [sts,ems,css,stg,rn,gid]=comp_TS_sub(ofl,n,Thr); % Each station

      CSs=[CSs;css];
      STs=[STs;sts];
      EMs=[EMs;ems];

      Stg=[Stg;stg];
      RN=[RN;rn];
      if n==length(ofl)
        Gid=[Gid;gid];
      end
    end
end

% All stations
if size(Stg,1)>1
  [sts,ems,css]=errM(Stg,Thr);
  RN=[RN;sprintf('Group_%i',Gid)];
  CSs=[CSs;css];
  STs=[STs;sts];
  EMs=[EMs;ems];
end

%% Label the results
if ~isempty(RN)
  STs=array2table(STs,'VariableNames',{'N','m_tg','m_rf','v_tg','v_rf'},'RowNames',RN);
  EMs=array2table(EMs,'VariableNames',{'RMS','CRMS','CC','NSE','KGE'},'RowNames',RN);
  if ~isempty(CSs)
    CSs=array2table(CSs,'VariableNames',{'r_h','r_m','r_f','r_n'},'RowNames',RN);
  end
end
end

function [sts,ems,css,stg,rn,gid]=comp_TS_sub(ofl,n,Thr)
%% Load the TSCls.m object
OTS=matfile(ofl{n});
OTS=OTS.OTS;
gid=OTS.gid;

%% Match target time resolution to reference
[Ttg,Trf]=OTS.UniTL('rf'); % Unify the time zones

if strcmp(OTS.TC1,'end') % Make the time lines to the centered convention
  Ttg=Ttg-OTS.TR1/2;
  Trf=Trf-OTS.TR2/2;
elseif strcmp(OTS.TC1,'begin')
  Ttg=Ttg+OTS.TR1/2;
  Trf=Trf+OTS.TR2/2;
end

% Aggregate the target time series
stg=[];
for t=1:length(Trf)
  k=Ttg>=Trf(t)-OTS.TR2/2 & Ttg<Trf(t)+OTS.TR2/2;
  if any(k)
    stg=[stg;[mean(OTS.TS1(k)) OTS.TS2(t)]];
  end
end

%% Perform the stats/error calculation
if size(stg,1)>1
  [sts,ems,css]=errM(stg,Thr);
  rn=OTS.Gtg.ID;
else
  sts=[];
  ems=[];
  css=[];
  rn='';
end
end
