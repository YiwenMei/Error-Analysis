% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 08/04/2018

%% Functionality
% This code performs error analysis between target data and reference record.
% Target data and reference record are inputted as time series. 

%% Input
% Stg/Srf: target/reference time series for all location as 1-by-N cell array
%          where each of the N cell stores the time series of a location;
% Gtg/Grf: geographic information of the target/reference time series as an N-by-3
%          table with the column name in order as "ID", "X" and "Y" where "ID"
%          is the name of the N target locations, "X" & "Y" are the x and y coordinates
%          of the locations;
% Ttg/Trf: time information of the target/reference time series as 1-by-N cell
%          array where each of the N cells stores the date number of the time
%          series;
% Rtg/Rrf: time resolution of target/reference time series in number of hour;
% Ctg/Crf: time convention of target/reference time series (can be "f", "b" or
%          "c" representing a forward, backward and centered time convention);
%   cf   : factor used to convert the unit of target data to reference;
%   md   : max distance between the matched target and reference location pair;
%   Thr  : thresholds used to calculate the contigency statistics for different
%          groups (set
%          it to "[]" if no need to output the contigency statistics);
%  outn  : variable name for the target time series outputted (set it to "[]"
%          if no need to output the time series);
%  opth  : output path to store the matched time series with time information
%          (set it to "[]" if no need to output the time series).

%% Output
% STs: statistics of target and reference time series;
% EMs: error metrics derived between target and reference time series;
% CSs: contigency statistics given threshold value "thr".

%% Additional note
% STs includes
%  1)sample size of matched target and reference time series pair,
%  2)mean of target time series,        3)mean of reference time series,
%  4)variance of target time series and 5)variance of reference time series;
% EMs includes
%  1)root mean square error,  2)centered root mean square error,
%  3)correlation coefficient, 4)Nash-Sutcliffe efficiency and
%  5)Kling-Gupta efficiency;
% CSs includes rate of 1)hit, 2)missing, 3)false alarm, and 4)correct negative.

function [STs,EMs,CSs]=comp_TS(Stg,Gtg,Ttg,Rtg,Ctg,Srf,Grf,Trf,Rrf,Crf,cf,md,Thr,outn,opth)
Rrf=Rrf/24; % Convert from hour to day
Rtg=Rtg/24;

%% Target inputted as time series
if strcmp(Ctg,'f') % Make the target record to the centered time convention
  Ttg=cellfun(@(x) x-Rtg/2,Ttg,'UniformOutput',false);
elseif strcmp(Ctg,'b')
  Ttg=cellfun(@(x) x+Rtg/2,Ttg,'UniformOutput',false);
end

% Match target time series to the closest stations
[id,d]=knnsearch([Gtg.X Gtg.Y],[Grf.X Grf.Y],'K',1);
id(d>md)=[];
Stg=Stg(id);
Ttg=Ttg(id);

clear id d md Gtg Rtg Ctg t

%% Error metrics and statistics
Gid=unique(Grf.Group); % Groupping of records
STs=cell(length(Gid),1);
CSs=cell(length(Gid),1);
EMs=cell(length(Gid),1);

for g=1:length(Gid)
  Grf_g=Grf(Grf.Group==g-1,:);
  Trf_g=Trf(Grf.Group==g-1);
  Srf_g=Srf(Grf.Group==g-1);

  Ttg_g=Ttg(Grf.Group==g-1);
  Stg_g=Stg(Grf.Group==g-1);

  STs{g}=nan(size(Grf_g,1)+1,5);
  CSs{g}=nan(size(Grf_g,1)+1,4);
  EMs{g}=nan(size(Grf_g,1)+1,5);

  Stg_a=[];
  Srf_a=[];
  for s=1:size(Grf_g,1)

% Aggregate target time series to match the time resolution of station records
    Stg_m=nan(length(Trf_g{s}),1);
    for t=1:length(Trf_g{s})
      if strcmp(Crf,'f') % Forward
        Stg_m(t)=nanmean(Stg_g{s}(Ttg_g{s}>=Trf_g{s}(t)-Rrf & Ttg_g{s}<Trf_g{s}(t)));
      elseif strcmp(Crf,'c') % Centered
        Stg_m(t)=nanmean(Stg_g{s}(Ttg_g{s}>=Trf_g{s}(t)-Rrf/2 & Ttg_g{s}<Trf_g{s}(t)+Rrf/2));
      elseif strcmp(Crf,'b') % Backward
        Stg_m(t)=nanmean(Stg_g{s}(Ttg_g{s}>=Trf_g{s}(t) & Ttg_g{s}<Trf_g{s}(t)+Rrf));
      end
    end
    Srf_m=Srf_g{s};
    k=~isnan(Stg_m) & ~isnan(Srf_m);
    Stg_m=cf*Stg_m(k); % Convert target time series unit to reference's
    Srf_m=Srf_m(k);
    Ttg_m=Trf_g{s}(k,1); % Matched Time Info of target
    Stg_a=[Stg_a;Stg_m]; % The groupped records
    Srf_a=[Srf_a;Srf_m];

    if ~isempty(opth) % Output the target time series
      save([opth outn '_' Grf_g.ID{s}],'Ttg_m','Stg_m');
    end

% Statistics and Error metrics of each station
    N=length(find(k)); % Sample size
    if N>=20/Rrf
      m_tg=mean(Stg_m); % Mean of target time series
      m_rf=mean(Srf_m); % Mean of reference time series
      v_tg=var(Stg_m,1); % Variance of target time series
      v_rf=var(Srf_m,1); % Variance of reference time series
      STs{g}(s,:)=[N m_tg m_rf v_tg v_rf];

      RMS=sqrt(mean((Stg_m-Srf_m).^2)); % Root mean square error
      CRMS=std(Stg_m-Srf_m,1); % Centered root mean square error
      CC=corr(Stg_m,Srf_m); % Correlation coefficient
      NSE=1-sum((Stg_m-Srf_m).^2)/(N*v_rf); % Nash Sutcliff efficiency
      KGE=1-sqrt((CC-1)^2+(m_tg/m_rf-1)^2+(sqrt(v_tg)/sqrt(v_rf)*m_rf/m_tg-1)^2); % KGE
      EMs{g}(s,:)=[RMS CRMS CC NSE KGE];

      if ~isempty(Thr) % Contigency statistics
        r_h=length(find(Stg_m>Thr(g) & Srf_m>Thr(g)))/N;
        r_m=length(find(Stg_m<=Thr(g) & Srf_m>Thr(g)))/N;
        r_f=length(find(Stg_m>Thr(g) & Srf_m<=Thr(g)))/N;
        r_n=length(find(Stg_m<=Thr(g) & Srf_m<=Thr(g)))/N;
        CSs{g}(s,:)=[r_h r_m r_f r_n]; 
      end
    end
  end

% Statistics and Error metrics of group
  N=length(Stg_a); % Sample size
  if N>=20*size(Grf_g,1)/Rrf && size(Grf_g,1)>0
    m_tg=mean(Stg_a); % Mean of target time series
    m_rf=mean(Srf_a); % Mean of reference time series
    v_tg=var(Stg_a,1); % Variance of target time series
    v_rf=var(Srf_a,1); % Variance of reference time series
    STs{g}(s+1,:)=[N m_tg m_rf v_tg v_rf];

    RMS=sqrt(mean((Stg_a-Srf_a).^2)); % Root mean square error
    CRMS=std(Stg_a-Srf_a,1); % Centered root mean square error
    CC=corr(Stg_a,Srf_a); % Correlation coefficient
    NSE=1-sum((Stg_a-Srf_a).^2)/(N*v_rf); % Nash Sutcliff efficiency
    KGE=1-sqrt((CC-1)^2+(m_tg/m_rf-1)^2+(sqrt(v_tg)/sqrt(v_rf)*m_rf/m_tg-1)^2); % KGE
    EMs{g}(s+1,:)=[RMS CRMS CC NSE KGE];

    if ~isempty(Thr) % Contigency statistics
      r_h=length(find(Stg_a>Thr(g) & Srf_a>Thr(g)))/N;
      r_m=length(find(Stg_a<=Thr(g) & Srf_a>Thr(g)))/N;
      r_f=length(find(Stg_a>Thr(g) & Srf_a<=Thr(g)))/N;
      r_n=length(find(Stg_a<=Thr(g) & Srf_a<=Thr(g)))/N;
      CSs{g}(s+1,:)=[r_h r_m r_f r_n]; 
    end
  end

  RN=[Grf_g.ID;mat2cell(['All_G' num2str(g-1,'%0i')],1,6)];
  RN=RN(~isnan(STs{g}(:,1)));
  STs{g}(isnan(STs{g}(:,1)),:)=[];
  STs{g}=table(STs{g}(:,1),STs{g}(:,2),STs{g}(:,3),STs{g}(:,4),STs{g}(:,5),'VariableNames',{'N','m_tg','m_rf','v_tg','v_rf'},'RowNames',RN);
  EMs{g}(isnan(EMs{g}(:,1)),:)=[];
  EMs{g}=table(EMs{g}(:,1),EMs{g}(:,2),EMs{g}(:,3),EMs{g}(:,4),EMs{g}(:,5),'VariableNames',{'RMS','CRMS','CC','NSE','KGE'},'RowNames',RN);
  if ~isempty(Thr) % Contigency statistics
    CSs{g}(isnan(CSs{g}(:,1)),:)=[];
    CSs{g}=table(CSs{g}(:,1),CSs{g}(:,2),CSs{g}(:,3),CSs{g}(:,4),'VariableNames',{'r_h','r_m','r_f','r_n'},'RowNames',RN);
  end
end
end
