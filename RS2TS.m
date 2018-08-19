% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 08/13/2018

%% Functionality
% This code extract time series for locations of interest over the target image
% stack of every time step.

%% Input
% Fn_tg: list of full name of target images for all time steps;
% TI_tg: an array of date number of the target images;
% f_tg : field name of target image if the target images are in .nc, .nc4, .hdf,
%        .hdf5 and .mat (for .mat file, f_tg is the Matlab variable name);
%  ndv : no data value of target image;
% GI_tg: x and y boundary and resolution of the target image (2-by-3 array where
%        row 1 stores the left boundary, right boundary and resolution of x and
%        row 2 stores the top boundary, bottom boundary and resolution of y);
% GI_rf: geographic information of the reference data as an N-by-4 table with
%        the column name in order as "ID", "X", "Y" and "Group" where "ID" is
%        the name of the N target locations, "X" & "Y" are the x and y coordinates
%        of the locations and "Group" is the grouping of the locations;
%  md  : max distance between the matched target and reference location pair
%        (to extract time series of grid cells that have at least one station
%        fall within, set it to the half length of the diagonal of the grid cell);
% opth : directory to store the extracted data time series with date number (set
%        it to "[]" if no need to output the time series);
% outn : name for the outputted time series (set it to "[]" if no need to output
%        the time series).

%% Output
% TS_tg: time series of interested locations extracted from the target image
%        stack stored as an 1-by-N cell array where each of the N cell stores
%        the time series of a location;
% GI_tg: geographic information of the target time series as an N-by-3 table
%        with the column name in order as "ID", "X" and "Y" where "ID" is the
%        names of the locations and "X" & "Y" are their x and y coordinates;
% TI_tg: time information of the target time series stored as an 1-by-N cell
%        array where each of the cell stores a time series of a location.

function [TS_tg,GI_tg,TI_tg]=RS2TS(Fn_tg,TI_tg,f_tg,ndv,GI_tg,GI_rf,md,opth,outn)
% Target domain
x=GI_tg(1,1)+GI_tg(1,3)/2:GI_tg(1,3):GI_tg(1,2)-GI_tg(1,3)/2;
y=GI_tg(2,1)-GI_tg(2,3)/2:-GI_tg(2,3):GI_tg(2,2)+GI_tg(2,3)/2;
[X,Y]=meshgrid(x,y);
X=reshape(X,numel(X),1);
Y=reshape(Y,numel(Y),1);

% Grid ID of interested locations
[id,d]=knnsearch([X Y],[GI_rf.X GI_rf.Y],'K',1);
id(d>md)=[];
GI_tg=table(GI_rf.ID,X(id),Y(id),GI_rf.Group,'VariableNames',{'ID','X','Y','Group'});
if ~isempty(opth) % Output the target time series
  save([opth 'GI_px' outn],'GI_tg');
end

% Read the target image data
TS_tg1=nan(size(Fn_tg,1),length(id));
for t=1:size(Fn_tg,1)
  [~,~,fex]=fileparts(Fn_tg(t,:));
  if strncmp(fex,'.tif',4) % compatable for .tiff
    Z=reshape(double(imread(Fn_tg(t,:))),1,numel(X));
  elseif strncmp(fex,'.nc4',3) % compatable for .nc
    Z=reshape(double(ncread(Fn_tg(t,:),f_tg)),1,numel(X));
  elseif strncmp(fex,'.hdf',4) % compatable for .hdf5
    Z=reshape(double(hdfread(Fn_tg(t,:),f_tg)),1,numel(X));
  elseif strcmp(fex,'.asc') || strcmp(fex,'.txt')
    Z=reshape(double(dlmread(Fn_tg(t,:),'',6,0)),1,numel(X));
  else
    load(Fn_tg(t,:),f_tg);
    Z=f_tg;
  end

% Extract time series of interested locations
  TS_tg1(t,:)=Z(id);
  TS_tg1(TS_tg1==ndv)=NaN;
end

% Output the target time series data with date number
if ~isempty(opth)
  for s=1:size(GI_tg,1)
    TS_tg=TS_tg1(:,s);

    save([opth 'px' outn '_' GI_tg.ID{s}],'TI_tg','TS_tg');
  end
end

TS_tg=mat2cell(TS_tg1,size(TS_tg1,1),ones(size(TS_tg1,2),1));
TI_tg=repmat(TI_tg,1,size(GI_tg,1));
TI_tg=mat2cell(TI_tg,size(TI_tg,1),ones(size(TI_tg,2),1));
end
