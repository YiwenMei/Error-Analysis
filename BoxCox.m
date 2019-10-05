% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/5/2019

%% Functionality
% This function performs the one(two)-parameter Box-Cox transformation of data.
%  It has a parameter calibration mode for estimating the parameter(s) and two
%  calculation modes for calculating the Box-Cox transformed data and its reversed
%  transformation.

%% Input
%  y  : input data for the Box-Cox transformation;
% Mod : mode of operation (can be 'PEst' for parameter calibration, 'A2BC' for
%       calculation of Box-Cox-transformed data, and 'BC2A' for the reversed
%       transformation);
% Prms: parameters or parameter ranges for the Box-Cox transformation.

%% Output
%  yTrs : transformed (or reversed tranformed) data;
% Prms_h: estimated Box-Cox parameter(s).

function [yTrs,Prms_h]=BoxCox(y,Mod,Prms)
%% Check the inputs
narginchk(3,3);
ips=inputParser;
ips.FunctionName=mfilename;
fprintf('%s received 3 required inputs\n',mfilename);

addRequired(ips,'y',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'y',1));
addRequired(ips,'Mod',@(x) any(strcmp(x,{'PEst','A2BC','BC2A'})));
addRequired(ips,'Prms',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'Prms',3));

parse(ips,y,Mod,Prms);
clear ips

%% Parameter estimation or transformed value calculation
if size(Prms,2)==1
  Prms=[Prms zeros(size(Prms))];
end

switch Mod
  case 'PEst' % Calibration mode
    Prms_i=mean(Prms,1);
    Prms_h=fminsearchbnd(@(Prms)BoxCoxOpt(y,Prms),Prms_i,Prms(2,:),Prms(1,:));

  case {'A2BC','BC2A'} % Calculation mode
    Prms_h=Prms;

  otherwise
    error('Mode can only PEst, A2BC or BC2A');
end

%% Box-Cox transformation
yTrs=BoxCoxTrs(y,Mod,Prms_h);
end

function [yTrs,yS]=BoxCoxTrs(y,Mod,Prms)
switch Mod
  case {'PEst','A2BC'}
%% Shifted variable values
    yS=y+Prms(2); % Apply the shift parameter

%% Box-Cox transformation
    if Prms(1)~=0 % Shape parameter not equals to 0
      yTrs=(yS.^Prms(1)-1)/Prms(1);
    else
      yTrs=log(yS);
    end

  case 'BC2A'
%% Reversed Box-Cox transformation
    if Prms(1)~=0 % Shape parameter not equals to 0
      yS=(Prms(1)*y+1).^1/Prms(1);
    else
      yS=exp(y);
    end

%% Shifted variable values
    yTrs=yS-Prms(2); % Apply the shift parameter
end
end

function LL=BoxCoxOpt(y,Mod,Prms)
%% Optimization
[yTrs,yS]=BoxCoxTrs(y,Mod,Prms);
LL=length(yTrs)/2*log(var(yTrs,1))-(Prms(1)-1)*sum(log(yS)); % -log likelihood function
end
