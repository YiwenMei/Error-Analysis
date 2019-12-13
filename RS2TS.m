% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/11/2019

%% Functionality
% This code is used to build the time series class (TSCls.m) objects for reference
%  locations. It contains two steps:
%   1)extract time series for the reference locations from a stack of spatial
%     raster images, and
%   2)associate the extracted time series with reference time series of the locations.

% The reference locations can be provided as points representing the ground gauge
%  locations or a mask of a basin. In the first case, reference data are supplied
%  as time series of those locations. In the latter one, data are supplied as
%  a space-time class (V2DTCls.m) object.

%% Input
% OStg: V2DTCls.m object for the target data of the study domain;
% rfT : type of reference data (can be either 'points' or 'shape');
% LOI : geographic information of the reference location(s)
%       - if rfT is 'points', LOI is a table with five fields (ID, X, Y, Cat,
%          and ofs) representing the name, X coordinate, Y coordinate, group
%          ID, and time offset (in number of hour to UTC) of the locations;
%       - if rfT is 'shape', LOI is a space map class (V2DCls.m) object for the
%          basin mask with 1/0 represents locations within/without the basin
%          (note that the mask mush be within the study domain);
% Drf : data of the reference location(s)
%       - if rfT is 'points', Drf is a cell array stores the time series data
%          for all the reference locations (within each cell is a table with
%          two columns for the date numbers and record values for one location);
%       - if rfT is 'shape', Drf is a V2DTCls.m object for the reference data
%         of the study domain;
% opth: directory to store the TSCls.m objects;

% pflg: parallel flag (false - default, use no parallel computing; true - otherwise);
%  cf : conversion factor for the target data.

%% Output
% OTfn: cell array stores the full name list of all the TSCls.m objects for the
%        reference location(s);

% TSCls.m objects for the location(s) are located in opth.

%% Additional note
% - Naming of the TSCls.m object:
%    The reference location ID and variable type are used to differentiate the
%    file names. For example, if there is a station in Grf called '1001' and
%    the variable type specified by OStg.vtp is 'Tair', the output file for the
%    object will be px_1001_Tair.mat in opth.
% - Require TSCls.m and V2DTCls.m. If rfT is 'shape', V2DCls.m is also required.

function OTfn=RS2TS(OStg,rfT,LOI,Drf,opth,varargin)
%% Check the inputs
narginchk(5,7);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'OStg',@(x) validateattributes(x,{'V2DTCls','Wind2DTCls'},{'nonempty'},...
    mfilename,'OStg'));
addRequired(ips,'rfT',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'rfT'));
addRequired(ips,'LOI',@(x) validateattributes(x,{'V2DCls','table'},{'nonempty'},...
    mfilename,'LOI'));
addRequired(ips,'Drf',@(x) validateattributes(x,{'V2DCls','cell'},{'nonempty'},...
    mfilename,'Drf'));
addRequired(ips,'opth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'opth'));

addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},...
    mfilename,'pflg'));
addOptional(ips,'cf',1,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'cf'));

parse(ips,OStg,rfT,LOI,Drf,opth,varargin{:});
pflg=ips.Results.pflg;
cf=ips.Results.cf;
clear ips varargin

%% Calculation with different inputs
switch rfT
  case 'shape'
    OTfn=RS2TS_sub1(OStg,LOI,cf,Drf,opth,pflg);

  case 'points'
    OTfn=RS2TS_sub2(OStg,LOI,cf,Drf,opth,pflg);
end
end

function OTfn=RS2TS_sub2(OStg,LOI,cf,Drf,opth,pflg)
%% ID of locations of interest Grids of the target image
[X,Y]=OStg.GridCls;
X=reshape(X,numel(X),1);
Y=reshape(Y,numel(Y),1);
x=LOI.X;
y=LOI.Y;

id=knnsearch([X Y],[x y],'K',1);

%% Read and extract time series from target image stack
Stg=[];
switch pflg
  case true
    parfor t=1:length(OStg.Fnm)
      Z=OStg.readCls(t)*cf;
      Stg=[Stg;Z(id)'];
    end

  case false
    for t=1:length(OStg.Fnm)
      Z=OStg.readCls(t)*cf;
      Stg=[Stg;Z(id)'];
    end
end

Ttg=OStg.TimeCls(OStg.TmR); % time line of target image stack in UTC

%% Construct the time series class
OTfn={};

X=X(id);
Y=Y(id);

for s=1:size(Stg,2)
  TS1=Stg(:,s);
  TL1=Ttg(~isnan(TS1));
  TS1=TS1(~isnan(TS1));

  TR2=Drf{s}.Properties.CustomProperties.T_Res;
  TC2=Drf{s}.Properties.CustomProperties.T_Conv;
  TL2=Drf{s}.Dnum;
  TS2=table2array(Drf{s}(:,2));
  TL2=TL2(~isnan(TS2));
  TS2=TS2(~isnan(TS2));

  Gtg.ID=LOI.ID{s};
  Gtg.X=X(s);
  Gtg.Y=Y(s);
  unt=Drf{s}.Properties.VariableUnits{2};
  ofn=fullfile(opth,sprintf('px_%s_%s.mat',LOI.ID{s},OStg.vtp));
  OTS=TSCls(OStg.vtp,cf*OStg.Ulm,cf*OStg.Llm,Gtg,LOI.Cat(s),OStg.ofs,TS1,TL1,OStg.TmR,...
      OStg.TmC,LOI.ofs(s),TS2,TL2,TR2,TC2,unt,{ofn});

  save(ofn,'OTS');
  OTfn=[OTfn;{ofn}];
end
end

function OTfn=RS2TS_sub1(OStg,OSmk,OSrf,cf,opth,pflg)
%% ID of mask on target and reference image stacks
if ~isempty(find(min([OStg.GIf(3,:);OSrf.GIf(3,:)],[],1)<OSmk.GIf(3,:), 1))
  error('Resolutions of mask must be finer than the target and reference images');
else
  [x,y]=OSmk.GridCls;
  mk=OSmk.readCls;
  x(mk==0)=[];
  y(mk==0)=[];
  x=reshape(x,numel(x),1);
  y=reshape(y,numel(y),1);
end

% Target image
[X1,Y1]=OStg.GridCls;
X1=reshape(X1,numel(X1),1);
Y1=reshape(Y1,numel(Y1),1);
id1=knnsearch([X1 Y1],[x y],'K',1);

% Reference image
[X2,Y2]=OSrf.GridCls;
X2=reshape(X2,numel(X2),1);
Y2=reshape(Y2,numel(Y2),1);
id2=knnsearch([X2 Y2],[x y],'K',1);

%% Read and extract the target image data
TL1=OStg.TimeCls(OStg.TmR); % time line
TL2=OSrf.TimeCls(OSrf.TmR);

TS1=[];
TS2=[];
switch pflg
  case true
    parfor t=1:length(OStg.Fnm)
      Z=OStg.readCls(t)*cf;
      TS1=[TS1;Z(id1)'];
      Z=OSrf.readCls(t);
      TS2=[TS2;Z(id2)'];
    end

  case false
    for t=1:length(OStg.Fnm)
      Z=OStg.readCls(t)*cf;
      TS1=[TS1;Z(id1)'];
      Z=OSrf.readCls(t);
      TS2=[TS2;Z(id2)'];
    end
end

%% Construct the time series class
[~,bfn,~]=OSmk.Fnm;
Gtg.ID=bfn;
Gtg.X=nanmean(x);
Gtg.Y=nanmean(y);

TS1=nansum(TS1,2);
N=length(find(TS1~=OStg.ndv));
TS1=TS1./N;
TS1(N/length(id1)<.5)=NaN;
TL1=TL1(~isnan(TS1));
TS1=TS1(~isnan(TS1));

TS2=nansum(TS2,2);
N=length(find(TS2~=OSrf.ndv));
TS2=TS2./N;
TS2(N/length(id2)<.5)=NaN;
TL2=TL2(~isnan(TS2));
TS2=TS2(~isnan(TS2));

OTfn=fullfile(opth,sprintf('ba_%s_%s.mat',bfn,OStg.vtp));
OTS=TSCls(OSrf.vtp,cf*OStg.Ulm,cf*OStg.Llm,Gtg,OSmk.gid,OStg.ofs,TS1,TL1,OStg.TmR,...
    OStg.TmC,OSrf.ofs,TS2,TL2,OSrf.TmR,OSrf.TmC,OSrf.unt,{OTfn});

save(OTfn,'OTS');
OTfn={OTfn};
end
