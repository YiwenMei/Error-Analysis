function OFl=VerInterp(ofl,mflg,varargin)
%% Check the inputs
narginchk(2,4);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'ofl',@(x) validateattributes(x,{'cell'},{'nonempty'},mfilename,'ofl'));
addRequired(ips,'mflg',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'mflg'));

addOptional(ips,'Hflg','Reference',@(x) validateattributes(x,{'char','double'},{'nonempty'},...
    mfilename,'Hflg'));
addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));

parse(ips,ofl,mflg,varargin{:});
Hflg=ips.Results.Hflg;
pflg=ips.Results.pflg;
clear ips varargin

%% Kernel of vertical interpolate
OFl=cell(size(ofl));
switch pflg
  case true
    parfor n=1:size(ofl,1)
      OFl(n,:)=VerInterp_sub(ofl(n,:),Hflg,mflg);
    end

  case false
    for n=1:size(ofl,1)
      [~,nm,~]=fileparts(ofl{n,1});
      fprintf('%i. %s - ',n,nm);
      OFl(n,:)=VerInterp_sub(ofl(n,:),Hflg,mflg);
      fprintf('Done\n');
    end
end
end

function ofl=VerInterp_sub(Ofl,Hflg,mflg)
ofl=cell(size(Ofl));
for i=1:size(Ofl,2)
% Load the TSCls.m object
  OTs=matfile(Ofl{i});
  OTs=OTs.OTS;

% Interpolate to the reference depth
  switch Hflg
    case 'Reference'
      OTs.TS1=VIfun(OTs.TS1,OTs.Vh1,OTs.Vh2,mflg);
      OTs.Vh1=OTs.Vh2;
      NL=length(OTs.Vh2);

    case 'Target'
      OTs.TS2=VIfun(OTs.TS2,OTs.Vh2,OTs.Vh1,mflg);
      OTs.Vh2=OTs.Vh1;
      NL=length(OTs.Vh1);
      
    otherwise
      fldn=fieldnames(OTs);
      k=strcmp(fldn,'Vh');
      Hfdn=fldn(k);
      k=strcmp(fldn,'TS');
      fldn=fldn(k);
      for s=1:2
        Stg_rf=VIfun(extractfield(OTs,fldn{s}),extractfield(OTs,Hfdn{s}),Hflg,mflg);
        OTS=setfield(OTs,fldn{s},Stg_rf);
        OTS=setfield(OTs,Hfdn{s},Hflg);
      end
      NL=length(Hflg);
  end
    
% Save the output
  ofl{i}=cell(NL,1);
  for l=1:NL
    if any(~isnan(OTs.TS2(:,l)))
      OTS=OTs;

      TL1=OTS.TL1;
      TS1=OTS.TS1(:,l);
      TL1(isnan(TS1))=[];
      TS1(isnan(TS1))=[];
      OTS.TL1=TL1;
      OTS.TS1=TS1;

      TL2=OTS.TL2;
      TS2=OTS.TS2(:,l);
      TL2(isnan(TS2))=[];
      TS2(isnan(TS2))=[];
      OTS.TL2=TL2;
      OTS.TS2=TS2;

      [pth,ofn,~]=fileparts(Ofl{i});
      ofn=fullfile(pth,sprintf('%s_p%02i.mat',ofn,l));
      save(ofn,'OTS');
      ofl{i}{l}=ofn;
    end
  end
end
end

function Stg_rf=VIfun(Stg,Htg,Hrf,mflg)
switch mflg
  case 'Linear'
    K=diff(Stg,1,2)./repmat(diff(Htg,1,2),size(Stg,1),1);
    K=[K(:,1) K K(:,end)];
    J=[0 Htg];
    dh=repmat(Hrf,length(J),1)-repmat(J',1,length(Hrf));
    [~,j]=min(abs(dh),[],1);
    J=j-1;
    J(J<1)=1;
    Stg_rf=K(:,j).*repmat(dh(j+(0:5:5*(length(j)-1))),size(Stg,1),1)+Stg(:,J);

  case 'Block'
    Stg_rf=nan(size(Stg,1),length(Hrf));
    for l=1:length(Hrf)
      Db=2*(sum(Hrf(l:-2:1))-sum(Hrf(l-1:-2:1))); % Depth of bottom
      dD=Db-Hrf(l); % Half-thickness
      k=Htg>=Hrf(l)-dD & Htg<Db;
      Stg_rf(:,l)=nanmean(Stg(:,k),2);
    end
end
end
