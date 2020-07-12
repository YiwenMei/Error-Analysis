% Yiwen Mei (yiwenm@umich.edu)
% SEAS, University of Michigan
% Last update: 5/3/2020

%% Functionality
% This code counts number of time steps with relative absolute error (RAE=|M-O|/O)
%  less than a given threshold.

%% Input
% OStg: space-time or profiled space time object for the target image stack;
% OSrf: space-time or profiled space time object for the reference image stack;
% Thr : thresholds used to quantify a acceptable RAE magnitude for spin up equilibrium;

% pflg: parallel flag (false/true - squential/parallel, default is false);
% Dthr: Number of day allow to no satisfy the equilibrium criteria for spin up.

%% Output
% Mk: mask showing pixels that do not passed the spin up criteria;
% ND: number of day that a pixel does not passed the spin up criteria.

%% Additional note
% Require V2DTCls.m, Wind2DTCls, or V3DTCls.m.

function [Mk,ND]=comp_spin(OStg,OSrf,Thr,varargin)
%% Check the inputs
narginchk(3,5);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'OStg',@(x) validateattributes(x,{'V2DTCls','Wind2DTCls','V3DTCls'},{'nonempty'},...
    mfilename,'OStg'));
addRequired(ips,'OSrf',@(x) validateattributes(x,{'V2DTCls','Wind2DTCls','V3DTCls'},{'nonempty'},...
    mfilename,'OSrf'));
addRequired(ips,'Thr',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Thr'));

addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));
addOptional(ips,'Dthr',0,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Dthr'));

parse(ips,OStg,OSrf,Thr,varargin{:});
pflg=ips.Results.pflg;
Dthr=ips.Results.Dthr;
clear ips varargin

%% Statistics and error metrics
Tcn=OSrf.TimeCls('begin');
[~,~,sz]=OSrf.GridCls;
ND=zeros(sz);
switch pflg
  case true
    parfor t=1:length(Tcn)
      Stg=OStg.readCls(t);
      Srf=OSrf.readCls(t);

      RAE=100*abs((Stg-Srf)./Srf);
      mk=RAE>Thr;
      ND=ND+mk;
    end

  case false
    for t=1:length(Tcn)
      Stg=OStg.readCls(t);
      Srf=OSrf.readCls(t);

      RAE=100*abs((Stg-Srf)./Srf);
      mk=RAE>Thr;
      ND=ND+mk;
    end
end
Mk=ND/length(Tcn)>Dthr;
end
