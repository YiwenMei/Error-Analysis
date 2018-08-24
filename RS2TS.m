% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 08/13/2018

%% Functionality
% This code extract time series for locations of interest over the target image
% stack of every time step.

%% Input
% Ftg : list of full name of target images for all time steps;
% Ttg : an array of date number of the target images;
% Ntg : field name of target image if the target images are in .nc, .nc4, .hdf,
%       .hdf5 and .mat (for .mat file, f_tg is the Matlab variable name);
% ndv : no data value of target image;
% Gtg : x and y boundary and resolution of the target image (2-by-3 array where
%       row 1 stores the left boundary, right boundary and resolution of x and
%       row 2 stores the top boundary, bottom boundary and resolution of y);
% Grf : geographic information of the reference data as an N-by-4 table with
%       the column name in order as "ID", "X", "Y" and "Group" where "ID" is
%       the name of the N target locations, "X" & "Y" are the x and y coordinates
%       of the locations and "Group" is the grouping of the locations;
%  md : max distance between the matched target and reference location pair (to
%       extract time series of grid cells that have at least one station fall
%       within, set it to the half length of the diagonal of the grid cell);
% opth: directory to store the extracted data time series with date number (set
%       it to "[]" if no need to output the time series);
% outn: name for the outputted time series (set it to "[]" if no need to output
%       the time series).

%% Output
% Stg: time series of interested locations extracted from the target image
%        stack stored as an 1-by-N cell array where each of the N cell stores
%        the time series of a location;
% Gtg: geographic information of the target time series as an N-by-3 table
%        with the column name in order as "ID", "X" and "Y" where "ID" is the
%        names of the locations and "X" & "Y" are their x and y coordinates;
% Ttg: time information of the target time series stored as an 1-by-N cell
%        array where each of the cell stores a time series of a location.

function [Stg,Gtg,Ttg]=RS2TS(Ftg,Ttg,Ntg,ndv,Gtg,Grf,md,opth,outn)
% Target domain
x=Gtg(1,1)+Gtg(1,3)/2:Gtg(1,3):Gtg(1,2)-Gtg(1,3)/2;
y=Gtg(2,1)-Gtg(2,3)/2:-Gtg(2,3):Gtg(2,2)+Gtg(2,3)/2;
[X,Y]=meshgrid(x,y);
X=reshape(X,numel(X),1);
Y=reshape(Y,numel(Y),1);

% Grid ID of interested locations
[id,d]=knnsearch([X Y],[Grf.X Grf.Y],'K',1);
id(d>md)=[];
Gtg=table(Grf.ID,X(id),Y(id),Grf.Group,'VariableNames',{'ID','X','Y','Group'});
if ~isempty(opth) % Output the target time series
  save([opth 'GI_px' outn],'Gtg');
end

% Read the target image data
Stg1=nan(size(Ftg,1),length(id));
for t=1:size(Ftg,1)
  [~,~,fex]=fileparts(Ftg(t,:));
  if strncmp(fex,'.tif',4) % compatable for .tiff
    Z=reshape(double(imread(Ftg(t,:))),1,numel(X));
  elseif strncmp(fex,'.nc4',3) % compatable for .nc
    Z=reshape(double(ncread(Ftg(t,:),Ntg)),1,numel(X));
  elseif strncmp(fex,'.hdf',4) % compatable for .hdf5
    Z=reshape(double(hdfread(Ftg(t,:),Ntg)),1,numel(X));
  elseif strcmp(fex,'.asc') || strcmp(fex,'.txt')
    Z=reshape(double(dlmread(Ftg(t,:),'',6,0)),1,numel(X));
  else
    load(Ftg(t,:),Ntg);
    Z=reshape(Ntg,1,numel(X));
  end

% Extract time series of interested locations
  Stg1(t,:)=Z(id);
end
Stg1(Stg1==ndv)=NaN;

% Output the target time series data with date number
if ~isempty(opth)
  for s=1:size(Gtg,1)
    Stg=Stg1(:,s);

    save([opth 'px' outn '_' Gtg.ID{s}],'Ttg','Stg');
  end
end

Stg=mat2cell(Stg1,size(Stg1,1),ones(size(Stg1,2),1));
Ttg=repmat(Ttg,1,size(Gtg,1));
Ttg=mat2cell(Ttg,size(Ttg,1),ones(size(Ttg,2),1));
end
