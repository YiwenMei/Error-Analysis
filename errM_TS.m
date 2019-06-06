function [sts,ems,css]=errM_TS(Z_rf,Z_tg,Z_thr)
%% Check the inputs
switch nargin
    case 1; error('Not enough inputs')
    case 2; Z_thr=[];
    case 3
    otherwise; error('Too many input arguments');
end

%% Produce time series
k=~isnan(Z_rf) & ~isnan(Z_tg); % Valid index
if isempty(find(k, 1))
  sts=nan(1,5);
  ems=nan(1,5);
  css=nan(1,5);

else
  N=length(find(k)); % Sample size

% Statistics
  m_tg=mean(Z_tg(k)); % Mean of target map
  m_rf=mean(Z_rf(k)); % Mean of reference map
  v_tg=var(Z_tg(k),1); % Variance of target map
  v_rf=var(Z_rf(k),1); % Variance of reference map
  sts=[N m_tg m_rf v_tg v_rf];

% Error metrics
  RMS=sqrt(mean((Z_tg(k)-Z_rf(k)).^2)); % Root mean square error
  CRMS=std(Z_tg(k)-Z_rf(k),1); % Centered root mean square error
  CC=corr(Z_tg(k),Z_rf(k)); % Correlation coefficient
  NSE=1-sum((Z_tg(k)-Z_rf(k)).^2)/(N*v_rf); % Nash Sutcliff efficiency
  KGE=1-sqrt((CC-1)^2+(m_tg/m_rf-1)^2+(sqrt(v_tg)/sqrt(v_rf)*m_rf/m_tg-1)^2); % Kling-Gupta efficiency
  ems=[RMS CRMS CC NSE KGE];

% Contigency statistics
  if ~isempty(Z_thr)
    r_h=length(find(Z_tg(k)>Z_thr & Z_rf(k)>Z_thr))/N; % Hit rate
    r_m=length(find(Z_tg(k)<=Z_thr & Z_rf(k)>Z_thr))/N; % Missing rate
    r_f=length(find(Z_tg(k)>Z_thr & Z_rf(k)<=Z_thr))/N; % False alarm rate
    r_n=length(find(Z_tg(k)<=Z_thr & Z_rf(k)<=Z_thr))/N; % Correct negative rate
    css=[r_h r_m r_f r_n];
  else
    css=nan(1,5);
  end
end
end
