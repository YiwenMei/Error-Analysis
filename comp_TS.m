% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 08/04/2018

%% Functionality
% This code performs error analysis between target data and reference record.
% Target data and reference record are inputted as time series. 

%% Input
% TS_tg : target time series for all location as 1-by-N cell array where each
%         of the N cell stores the time series of a location
%   cf  : unit conversion factor for target data;
% GI_tg : geographic information of the target data as an N-by-3 table with the
%         column name in order as "ID", "X" and "Y" where "ID" is the name of
%         the N target locations, "X" & "Y" are the x and y coordinates of the
%         locations;
% TI_tg : time information of the target data as 1-by-N cell array where each
%         of the N cells stores the date number of the time series;
% Trs_tg: time resolution of target data in number of hour;
% Tcv_tg: time convention of target data (can be "f", "b" or "c" representing
%         a forward, backward and centered time convention;
% TS_rf : reference time series as 1-by-N cell array where each of the N cell
%         stores the time series of a location of interest;
% GI_rf : geographic information of the reference data as an N-by-4 table with
%         the column name in order as "ID", "X", "Y" and "Group" where "ID" is
%         the name of the N target locations, "X" & "Y" are the x and y coordinates
%         of the locations and "Group" is the grouping of the locations;
% TI_rf : time information of the reference data as 1-by-N cell array where each
%         of the N cells stores the date number of the time series
% Trs_rf: time resolution of reference data in number of hour;
% Tcv_rf: time convention of reference data (see Tcv_tg);
%   md  : max distance between the matched target and reference location pair;
%  outn : variable name for the target time series outputted (set it to "[]"
%         if no need to output the time series;
%  opth : output path to store the matched time series with time information
%         (set it to "[]" if no need to output the time series).

%% Output
% STs: statistics of target and reference time series including 1)sample size
%      of matched target and reference time series pair, 2)mean of target time
%      series, 3)mean of reference time series, 4)variance of target time series
%      and 5)variance of reference time series;
% EMs: error metrics derived between target and reference time series including
%      1)root mean square error, 2)centered root mean square error, 3)correlation
%      coefficient, 4)Nash-Sutcliffe efficiency and 5)Kling-Gupta efficiency.

function [STs,EMs]=comp_TS(TS_tg,cf,GI_tg,TI_tg,Trs_tg,Tcv_tg,TS_rf,GI_rf,TI_rf,Trs_rf,Tcv_rf,md,outn,opth)
Trs_rf=Trs_rf/24; % Convert from hour to day
Trs_tg=Trs_tg/24; % Convert from hour to day

%% Target inputted as time series
if strcmp(Tcv_tg,'f') % Make the target record to the centered time convention
  TI_tg=cellfun(@(x) x-Trs_tg/2,TI_tg,'UniformOutput',false);
elseif strcmp(Tcv_tg,'b')
  TI_tg=cellfun(@(x) x+Trs_tg/2,TI_tg,'UniformOutput',false);
end

% Match target time series to the closest stations
[id,d]=knnsearch([GI_tg.X GI_tg.Y],[GI_rf.X GI_rf.Y],'K',1);
id(d>md)=[];
TS_tg=TS_tg(id);
TI_tg=TI_tg(id);

clear id d md GI_tg Trs_tg Tcv_tg t

%% Error metrics and statistics
Gid=unique(GI_rf.Group); % Groupping of records
STs=cell(length(Gid),1);
EMs=cell(length(Gid),1);

for g=1:length(Gid)
  GIg_rf=GI_rf(GI_rf.Group==g-1,:);
  TIg_rf=TI_rf(GI_rf.Group==g-1);
  TSg_rf=TS_rf(GI_rf.Group==g-1);

  TIg_tg=TI_tg(GI_rf.Group==g-1);
  TSg_tg=TS_tg(GI_rf.Group==g-1);

  STs{g}=nan(size(GIg_rf,1)+1,5);
  EMs{g}=nan(size(GIg_rf,1)+1,5);

  Ttg=[];
  Trf=[];
  for s=1:size(GIg_rf,1)

% Aggregate target time series to match the time resolution of station records
    Fm_tg=nan(length(TIg_rf{s}),1);
    for t=1:length(TIg_rf{s})
      if strcmp(Tcv_rf,'f') % Forward
        Fm_tg(t)=nanmean(TSg_tg{s}(TIg_tg{s}>=TIg_rf{s}(t)-Trs_rf & TIg_tg{s}<TIg_rf{s}(t)));
      elseif strcmp(Tcv_rf,'c') % Centered
        Fm_tg(t)=nanmean(TSg_tg{s}(TIg_tg{s}>=TIg_rf{s}(t)-Trs_rf/2 & TIg_tg{s}<TIg_rf{s}(t)+Trs_rf/2));
      elseif strcmp(Tcv_rf,'b') % Backward
        Fm_tg(t)=nanmean(TSg_tg{s}(TIg_tg{s}>=TIg_rf{s}(t) & TIg_tg{s}<TIg_rf{s}(t)+Trs_rf));
      end
    end
    TSm_rf=TSg_rf{s};
    k=~isnan(Fm_tg) & ~isnan(TSm_rf);
    Fm_tg=cf*Fm_tg(k); % Convert target time series unit to reference's
    TSm_rf=TSm_rf(k);
    TIm_tg=TIg_rf{s}(k,1); % Matched Time Info of target
    Ttg=[Ttg;Fm_tg]; % The groupped records
    Trf=[Trf;TSm_rf];

    if ~isempty(opth) % Output the target time series
      save([opth outn '_' GIg_rf.ID{s}],'TIm_tg','Fm_tg');
    end

% Statistics and Error metrics of each station
    N=length(find(k)); % Sample size
    if N>=20/Trs_rf
      m_tg=mean(Fm_tg); % Mean of target time series
      m_rf=mean(TSm_rf); % Mean of reference time series
      v_tg=var(Fm_tg,1); % Variance of target time series
      v_rf=var(TSm_rf,1); % Variance of reference time series
      STs{g}(s,:)=[N m_tg m_rf v_tg v_rf];

      RMS=sqrt(mean((Fm_tg-TSm_rf).^2)); % Root mean square error
      CRMS=std(Fm_tg-TSm_rf,1); % Centered root mean square error
      CC=corr(Fm_tg,TSm_rf); % Correlation coefficient
      NSE=1-sum((Fm_tg-TSm_rf).^2)/(N*v_rf); % Nash Sutcliff efficiency
      KGE=1-sqrt((CC-1)^2+(m_tg/m_rf-1)^2+(sqrt(v_tg)/sqrt(v_rf)*m_rf/m_tg-1)^2); % KGE
      EMs{g}(s,:)=[RMS CRMS CC NSE KGE];
    end
  end

% Statistics and Error metrics of group
  N=length(Ttg); % Sample size
  if N>=20*size(GIg_rf,1)/Trs_rf && size(GIg_rf,1)>0
    m_tg=mean(Ttg); % Mean of target time series
    m_rf=mean(Trf); % Mean of reference time series
    v_tg=var(Ttg,1); % Variance of target time series
    v_rf=var(Trf,1); % Variance of reference time series
    STs{g}(s+1,:)=[N m_tg m_rf v_tg v_rf];

    RMS=sqrt(mean((Ttg-Trf).^2)); % Root mean square error
    CRMS=std(Ttg-Trf,1); % Centered root mean square error
    CC=corr(Ttg,Trf); % Correlation coefficient
    NSE=1-sum((Ttg-Trf).^2)/(N*v_rf); % Nash Sutcliff efficiency
    KGE=1-sqrt((CC-1)^2+(m_tg/m_rf-1)^2+(sqrt(v_tg)/sqrt(v_rf)*m_rf/m_tg-1)^2); % KGE
    EMs{g}(s+1,:)=[RMS CRMS CC NSE KGE];
  end

  RN=[GIg_rf.ID;mat2cell(['All_G' num2str(g-1,'%0i')],1,6)];
  RN=RN(~isnan(STs{g}(:,1)));
  STs{g}(isnan(STs{g}(:,1)),:)=[];
  EMs{g}(isnan(EMs{g}(:,1)),:)=[];

  STs{g}=table(STs{g}(:,1),STs{g}(:,2),STs{g}(:,3),STs{g}(:,4),STs{g}(:,5),'VariableNames',{'N','m_tg','m_rf','v_tg','v_rf'},'RowNames',RN);
  EMs{g}=table(EMs{g}(:,1),EMs{g}(:,2),EMs{g}(:,3),EMs{g}(:,4),EMs{g}(:,5),'VariableNames',{'RMS','CRMS','CC','NSE','KGE'},'RowNames',RN);
end
end
