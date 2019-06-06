function [dPr,Zu,Zl]=RemOut_2D(Z,thr,pfg)
% Check inputs
switch nargin
    case {1:2}; error('Not enough arguments');
    case 3
    otherwise; error('Too many number of arguments');
end

% Statistics of the image
Z=Z(~isnan(Z));
Zm=mean(Z);
Zs=std(Z);

% Calculate the change of removal
Nr=nan(size(thr));
for i=1:length(thr)
  Zu=Zm+thr(i)*Zs;
  Zl=Zm-thr(i)*Zs;
  Nr(i)=length(find(Z>Zu | Z<Zl)); % total number of removal
end
dNr=diff(Nr)./diff(thr); % change of removal with threshold
Nr=movmean(Nr,2); % mean removal
Nr(1)=[];
dPr=dNr./Nr; % change of removal in percent with threshold

% Find the image thresholds
[~,i]=min(dPr); % Minimum of % change of removal
thrm=max(thr(i),3); % At least 3 times of the SD
Zu=Zm+thrm*Zs;
Zl=Zm-thrm*Zs;

% Plotting
if pfg==1
  thrm=movmean(thr,2);
  thrm(1)=[];
  plot(thrm,dPr);
end
end
