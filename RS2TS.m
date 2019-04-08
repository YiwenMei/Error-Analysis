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
%       .hdf5 and .mat (for .mat file, Ntg is the Matlab variable name);
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
% opth: directory to store the "Stg", "Gtg" and "Ttg" (see Output section for
%       the meanings); set it to "[]" if no need to output the time series;
% outn: string included in the outputted files for "Stg", "Gtg" and "Ttg" (see
%       additional note for more details); set it to "[]" if no need to output
%       the time series.

%% Output
% Stg: time series of interested locations extracted from the target image stack
%      stored as an 1-by-N cell array where each of the N cell stores the time
%      series of a location;
% Gtg: geographic information of the target time series as an N-by-3 table with
%      the column name in order as "ID", "X" and "Y" where "ID" is the names
%      of the locations and "X" & "Y" are their x and y coordinates;
% Ttg: time information of the target time series stored as an 1-by-N cell array
%      where each of the cell stores a time series of a location.

%% Additional note
% Naming of outputs:
% - Gtg is stored in a .mat file located in opth with the name "GI_pxoutn.mat".
%   For example, if opth is 'C:\Users\user\Desktop\Work\HiMat\Data\MTS\CHP' and
%   outn is 'CHIRP', then the full name of file stores Gtg will be
%   C:\Users\user\Desktop\Work\HiMat\Data\MTS\CHP\GI_pxCHIRP.mat
% - Stg and Ttg are store together in .mat files for all locations. The station
%   IDs included in Grf are used to differentiate the file names. For example,
%   if there is a station in "Grf" called '1001' and with the same opth and outn
%   in the previous example, then the name of the file stores Stg and Ttg of
%   station 1001 will be
%   C:\Users\user\Desktop\Work\HiMat\Data\MTS\CHP\pxCHIRP_1001.mat

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

% Output the target time series geolocation information
if ~isempty(opth)
  save([opth 'GI_px' outn],'Gtg');
end

% Read the target image data
Stg1=nan(size(Ftg,1),length(id));
parfor t=1:size(Ftg,1)
  Z=read2Dvar(Ftg(t,:),Ntg,ndv);

% Extract time series of interested locations
  Stg1(t,:)=Z(id);
end
Stg1(Stg1==ndv)=NaN;

% Output the target time series with date number
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
