% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 1/2/2021

%% Functionality
% This code is used to build the time series class (TSCls.m) objects for reference
%  locations. It contains two steps:
%   1)extract time series for the reference locations from the target spatial
%     raster image stack, and
%   2)associate the extracted time series with the reference time series.

% Reference geographic information can be provided as points representing the
%  ground gauge locations or masks of basins. In the first case, reference data
%  are supplied as time series of those locations. In the latter one, data are
%  supplied as a space-time class (V2DTCls.m) object or time series stored in
%  cell array.

%% Input
% OStg: V2DTCls.m object for the target data of the study domain;
% rfT : type of reference data (can be 'points', 'shape-imgs' or 'shape-TS');
% LOI : geographic information of the reference location(s)
%       - if rfT is 'points', LOI is a table with five fields (ID, X, Y, Cat,
%          and ofs) representing the name, X coordinate, Y coordinate, group
%          ID, and time offset (in number of hour to UTC) of the locations;
%       - if rfT is 'shape-X', LOI is the shapefile names of the study basins;
% Drf : data of the reference location(s)
%       - if rfT is 'points', Drf is a cell array stores the time series data
%          for all the reference locations (within each cell is a table with
%          two columns for the date numbers and record values for one location);
%       - if rfT is 'shape-imgs', Drf is a V2DTCls.m object for the reference data
%         of the study domain;
%       - if rfT is 'shape-TS', Drf is a cell array stores the basin-average
%         time series of the basin shapes;
% opth: directory to store the TSCls.m objects.
% xTfn: Extracted time series for location of interests.

% pflg: parallel flag (false - default, use no parallel computing; true - otherwise);
%  cf : conversion factor for the target data (deflaut is 1).

%% Output
% OTfn: cell array stores the full name list of all the TSCls.m objects for the
%        reference location(s);

% Extracted time series for the location of interests are stored as xTfn;
% Paired time series in the name of pTS.lid.asn.mat as a TSCls.m objects for
%  the location(s) are located in opth (see additional note for more details).

%% Additional note
% - Naming of the TSCls.m object:
%    The reference location ID and variable type are used to differentiate the
%    file names. For example, if there is a station in LOI called '1001' and
%    the variable type specified by OStg.vtp is 'Tair', the output file for the
%    object will be pTS.1001.Tair.mat in opth.
% - Require TSCls.m, V2DTCls.m, and V3DTCls.

function OTfn=RS2TS(OStg,rfT,LOI,Drf,opth,xTfn,varargin)
%% Check the inputs
narginchk(6,8);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'OStg',@(x) validateattributes(x,{'V2DTCls','Wind2DTCls','V3DTCls'},{'nonempty'},...
    mfilename,'OStg'));
addRequired(ips,'rfT',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'rfT'));
switch rfT
  case 'points'
    addRequired(ips,'LOI',@(x) validateattributes(x,{'table'},{'nonempty'},mfilename,'LOI'));
    addRequired(ips,'Drf',@(x) validateattributes(x,{'cell'},{'nonempty'},mfilename,'Drf'));
  case 'shape-imgs'
    addRequired(ips,'LOI',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'LOI'));
    addRequired(ips,'Drf',@(x) validateattributes(x,{'V2DTCls','Wind2DTCls','V3DTCls'},{'nonempty'},...
        mfilename,'Drf'));
  case 'shape-TS'
    addRequired(ips,'LOI',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'LOI'));
    addRequired(ips,'Drf',@(x) validateattributes(x,{'cell'},{'nonempty'},mfilename,'Drf'));
end
addRequired(ips,'opth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'opth'));
addRequired(ips,'xTfn',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'xTfn'));

addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));
addOptional(ips,'cf',1,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'cf'));

parse(ips,OStg,rfT,LOI,Drf,opth,xTfn,varargin{:});
pflg=ips.Results.pflg;
cf=ips.Results.cf;
clear ips varargin

%% Calculation with different inputs
switch rfT
  case 'points'
    OTfn=RS2TS_sub1(OStg,LOI,cf,Drf,opth,xTfn,pflg);

  case 'shape-imgs'
    OTfn=RS2TS_sub2(OStg,LOI,cf,Drf,opth,xTfn,pflg);

  case 'shape-TS'
    OTfn=RS2TS_sub3(OStg,LOI,cf,Drf,opth,xTfn,pflg);
end
end

function OTfn=RS2TS_sub1(OStg,LOI,cf,Drf,opth,xTfn,pflg)
%% ID of locations of interest Grids of the target image
[X,Y]=OStg.GridCls;
X=reshape(X,numel(X),1);
Y=reshape(Y,numel(Y),1);
x=LOI.Xcor;
y=LOI.Ycor;
id=knnsearch([X Y],[x y],'K',1);
clear x y

%% Read and extract time series from target image stack
if exist(xTfn,'file')~=2 % Extracted time series not existed
  TC2=Drf{1}.Properties.CustomProperties.T_Conv;
  Ttg=OStg.TimeCls(TC2); % time line of target image stack in UTC

  if isa(OStg,'V3DTCls')
    NL=length(OStg.VHt);
    fprintf('The target image is profiled\n');
  elseif isa(OStg,'V2DTCls')
    NL=1;
    fprintf('The target image is single-level\n');
  end
  fprintf('Start to extract time series from target image stack\n');

  Stg=nan(length(Ttg),length(id),NL);
  switch pflg
    case true
      parfor t=1:length(OStg.Fnm)
        Z=OStg.readCls(t)*cf;
        Stg(t,:,:)=Z(id+length(X)*(0:NL-1));
      end

    case false
      for t=1:length(OStg.Fnm)
        Z=OStg.readCls(t)*cf;
        Stg(t,:,:)=Z(id+length(X)*(0:NL-1));
      end
  end
  save(xTfn,'Stg','Ttg');
  Wmsg=lastwarn;
  if ~isempty(Wmsg)
    save(xTfn,'Stg','Ttg','-v7.3');
  end
  fprintf('Extracted time series saved to %s\n',xTfn);

else
  load(xTfn,'Stg','Ttg');
  fprintf('Matched time series file %s loaded\n',xTfn);
end

%% Construct the time series class
OTfn={};
X=X(id);
Y=Y(id);
for s=1:size(Stg,2)
  TS1=permute(Stg(:,s,:),[1 3 2]);
  TL1=Ttg(any(~isnan(TS1),2));
  TS1=TS1(any(~isnan(TS1),2),:);

  TR2=Drf{s}.Properties.CustomProperties.T_Res;
  TC2=Drf{s}.Properties.CustomProperties.T_Conv;
  TL2=Drf{s}.Dnum;
  TS2=table2array(Drf{s}(:,2:end));
  TL2=TL2(any(~isnan(TS2),2));
  TS2=TS2(any(~isnan(TS2),2),:);

  Gtg.Sid=LOI.Sid{s};
  Gtg.Xcor=X(s);
  Gtg.Ycor=Y(s);
  unt=Drf{s}.Properties.VariableUnits{2};
  if isa(OStg,'V3DTCls')
    vh2=Drf{s}.Properties.CustomProperties.VerHgt;
    OTS=TSCls(OStg.vtp,cf*OStg.Ulm,cf*OStg.Llm,Gtg,LOI.Gid(s),OStg.ofs,TS1,TL1,OStg.TmR,OStg.TmC,...
        LOI.Tofs(s),TS2,TL2,TR2,TC2,OStg.VHt,vh2(2:end),unt);
  elseif isa(OStg,'V2DTCls')
    OTS=TSCls(OStg.vtp,cf*OStg.Ulm,cf*OStg.Llm,Gtg,LOI.Gid(s),OStg.ofs,TS1,TL1,OStg.TmR,OStg.TmC,...
        LOI.Tofs(s),TS2,TL2,TR2,TC2,0,0,unt);
  end

  ofn=fullfile(opth,sprintf('pTS.%s.%s.mat',LOI.Sid{s},OStg.vtp));
  save(ofn,'OTS');
  fprintf('Matched "%s" time series pair for station "%s" saved\n',OStg.vtp,LOI.Sid{s});
  OTfn=[OTfn;{ofn}];
end
end

function OTfn=RS2TS_sub2(OStg,sfn,cf,OSrf,opth,xTfn,pflg)
%% Basin shapefile
fprintf('Start to extract data for basins');
Bfl=shaperead(sfn);
Bfl=struct2table(Bfl);
idpth=fileparts(sfn);
for i=1:size(Bfl,1)
  Bid=Bfl.Sid{i};
  bfn=fullfile(idpth,sprintf('id-%s.mat',Bid));

  if exist(bfn,'file')~=2
% Target image
    [X,Y]=OStg.GridCls;
    [id1,id]=inpolygon(reshape(X,numel(X),1), reshape(Y,numel(Y),1), Bfl.X{i}, Bfl.Y{i});
    id1=any([id1 id],2);

% Reference image
    [X,Y]=OSrf.GridCls;
    [id2,id]=inpolygon(reshape(X,numel(X),1), reshape(Y,numel(Y),1), Bfl.X{i}, Bfl.Y{i});
    id2=any([id2 id],2);
    x=mean(X(id2));
    y=mean(Y(id2));

% Output the id file
    Gid=Bfl.Gid{i};
    save(bfn,'id1','id2','x','y','Bid','Gid');
  end
end
clear X Y id id1 id2 Bid Gid x y
fprintf(' --- Done\n');

%% Calculate the basin average value
if exist(xTfn,'file')~=2 % Extracted time series not existed
  Ttg=OStg.TimeCls(OSrf.TmC); % time line of target image stack in UTC
  Trf=OSrf.TimeCls(OSrf.TmC); % Note that time convention converted to the reference one

  fprintf('Start to calculate basin average time series from target image stack\n');
  Stg=nan(length(OStg.Fnm),size(Bfl,1));
  Srf=nan(length(OStg.Fnm),size(Bfl,1));
  switch pflg
    case true
      parfor t=1:length(OStg.Fnm)
        ws_tg=OStg.readCls(t)*cf;
        ws_rf=OSrf.readCls(t);
        [Ts1,Ts2]=RS2TS_sub2_sub(ws_tg,ws_rf,Bfl,idpth);
        Stg(t,:)=Ts1; 
        Srf(t,:)=Ts2; 
      end

    case false
      for t=1:length(OStg.Fnm)
        ws_tg=OStg.readCls(t)*cf;
        ws_rf=OSrf.readCls(t);
        [Ts1,Ts2]=RS2TS_sub2_sub(ws_tg,ws_rf,Bfl,idpth);
        Stg(t,:)=Ts1; 
        Srf(t,:)=Ts2; 
      end
  end
  clear Ts1 Ts2
  save(xTfn,'Stg','Srf','Ttg','Trf');
  fprintf('Basin average "%s" saved to %s\n',OStg.vtp,xTfn);

else
  load(xTfn,'Stg','Srf','Ttg','Trf');
  fprintf('Basin average time series file %s loaded\n',xTfn);
end

%% Construct the time series class
OTfn={};
for s=1:size(Bfl,1)
  TS1=Stg(:,s);
  TL1=Ttg(~isnan(TS1));
  TS1=TS1(~isnan(TS1));
  TS2=Srf(:,s);
  TL2=Trf(~isnan(TS2));
  TS2=TS2(~isnan(TS2));

  bs=matfile(Bfl.ID_name{s});
  Gtg.Sid=bs.Bid;
  Gtg.Xcor=bs.x;
  Gtg.Ycor=bs.y;
  OTS=TSCls(OSrf.vtp, cf*OStg.Ulm, cf*OStg.Llm, Gtg, bs.Gid, 0, TS1, TL1, OStg.TmR, OSrf.TmC,...
      0, TS2, TL2, OSrf.TmR, OSrf.TmC, 0,0,OSrf.unt);

  ofn=fullfile(opth,sprintf('pTS.%s.%s.mat',Gtg.ID,OStg.vtp));
  save(ofn,'OTS');
  fprintf('Basin average "%s" time series pair for basin "%s" saved\n',OStg.vtp,Gtg.ID);
  OTfn=[OTfn;{ofn}];
end
end

function OTfn=RS2TS_sub3(OStg,sfn,cf,Drf,opth,xTfn,pflg)
%% Basin shapefile
fprintf('Start to extract data for basins');
Bfl=shaperead(sfn);
Bfl=struct2table(Bfl);
idpth=fileparts(sfn);
for i=1:size(Bfl,1) % Loop over the basins
  Bid=Bfl.Sid{i};
  bfn=fullfile(idpth,sprintf('id-%s.mat',Bid));

  if exist(bfn,'file')~=2
% Target image
    [X,Y]=OStg.GridCls;
    [id1,id]=inpolygon(reshape(X,numel(X),1), reshape(Y,numel(Y),1), Bfl.X{i}, Bfl.Y{i});
    id1=any([id1 id],2);

    x=mean(X(id1)); % Centroid
    y=mean(Y(id1));

% Output the id file
    Gid=Bfl.Gid{i};
    save(bfn,'id1','x','y','Bid','Gid');
  end
end
clear X Y id id1 Bid Gid x y
fprintf(' --- Done\n');

%% Calculate the basin average value
if exist(xTfn,'file')~=2 % Extracted time series not existed
  TC2=Drf{1}.Properties.CustomProperties.T_Conv;
  Ttg=OStg.TimeCls(TC2); % time line of target image stack in UTC

  fprintf('Start to calculate basin average time series from target image stack\n');
  Stg=nan(length(OStg.Fnm),size(Bfl,1));
  switch pflg
    case true
      parfor t=1:length(OStg.Fnm)
        ws_tg=OStg.readCls(t)*cf;
        Ts=RS2TS_sub2_sub(ws_tg,[],Bfl,idpth);
        Stg(t,:)=Ts; 
      end

    case false
      for t=1:length(OStg.Fnm)
        ws_tg=OStg.readCls(t)*cf;
        Ts=RS2TS_sub2_sub(ws_tg,[],Bfl,idpth);
        Stg(t,:)=Ts; 
      end
  end
  clear Ts
  save(xTfn,'Stg','Ttg');
  fprintf('Basin average "%s" saved to %s\n',OStg.vtp,xTfn);

else
  load(xTfn,'Stg','Ttg');
  fprintf('Basin average time series file %s loaded\n',xTfn);
end

%% Construct the time series class
OTfn={};
for s=1:size(Bfl,1)
  TS1=Stg(:,s);
  TL1=Ttg(~isnan(TS1));
  TS1=TS1(~isnan(TS1));

  TR2=Drf{s}.Properties.CustomProperties.T_Res;
  TC2=Drf{s}.Properties.CustomProperties.T_Conv;
  TL2=Drf{s}.Dnum;
  TS2=table2array(Drf{s}(:,2:end));
  TL2=TL2(any(~isnan(TS2),2));
  TS2=TS2(any(~isnan(TS2),2),:);

  bfn=fullfile(idpth,sprintf('id-%s.mat',Bfl.Sid{s}));
  bs=matfile(bfn);
  Gtg.Sid=bs.Bid;
  Gtg.Xcor=bs.x;
  Gtg.Ycor=bs.y;
  unt=Drf{s}.Properties.VariableUnits{2};
  if isa(OStg,'V3DTCls')
    vh2=Drf{s}.Properties.CustomProperties.VerHgt;
    OTS=TSCls(OStg.vtp, cf*OStg.Ulm, cf*OStg.Llm, Gtg, bs.Gid, 0, TS1, TL1, OStg.TmR, OStg.TmC,...
        0, TS2, TL2, TR2, TC2, OStg.VHt, vh2(2:end), unt);
  elseif isa(OStg,'V2DTCls')
    OTS=TSCls(OStg.vtp, cf*OStg.Ulm, cf*OStg.Llm, Gtg, bs.Gid, 0, TS1, TL1, OStg.TmR, OStg.TmC,...
        0, TS2, TL2, TR2, TC2, 0,0,unt);
  end

  ofn=fullfile(opth,sprintf('pTS.%s.%s.mat',Gtg.Sid,OStg.vtp));
  save(ofn,'OTS');
  fprintf('Basin average "%s" time series pair for basin "%s" saved\n',OStg.vtp,Gtg.Sid);
  OTfn=[OTfn;{ofn}];
end
end

function [Ts1,Ts2]=RS2TS_sub2_sub(Z1,Z2,bfl,idpth)
Ts1=nan(1,size(bfl,1));
Ts2=nan(1,size(bfl,1));
for i=1:size(bfl,1)
  bfn=fullfile(idpth,sprintf('id-%s.mat',bfl.Sid{i}));
  id1=matfile(bfn);
  if ~isempty(Z2)
    id2=id1.id2;
  end
  id1=id1.id1;

  z=Z1(id1);
  N=sum(~isnan(z));
  if N/length(z)>.6 % Threshold of available data
    z=nansum(z)/N;
  else
    z=NaN;
  end
  Ts1(i)=z;

  if ~isempty(Z2)
    z=Z2(id2);
    N=sum(~isnan(z));
    if N/length(z)>.6
      z=nansum(z)/N;
    else
      z=NaN;
    end
    Ts2(i)=z;
  end
end
end
