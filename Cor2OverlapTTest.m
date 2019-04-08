function [H,Ha,pv]=Cor2OverlapTTest(Sref,S1,S2,alp)
%% Compute the correlation
[cc,pv]=corr([Sref S1 S2],'rows','complete');
cr1=cc(1,2); % sample 1 to ref
cr2=cc(1,3); % sample 2 to ref
c12=cc(2,3); % sample 1 to sample 2
N=length(find(~isnan(Sref) & ~isnan(S1) & ~isnan(S2))); % Sample size

%% Compute the test statistic
zr1=log((1+cr1)/(1-cr1))/2; % Fisher transform
zr2=log((1+cr2)/(1-cr2))/2;
c_2m=(cr1^2+cr2^2)/2;
f=min([(1-c12)/(2*(1-c_2m)) 1]);
h=1+c_2m*(1-f)/(1-c_2m);
ts=(zr1-zr2)*sqrt((N-3)/(2*(1-c12)*h));

%% Hypothesis test
y=1-abs((tcdf(ts,N-3)-tcdf(-ts,N-3)));
H=y<alp;
pts=tinv(1-alp/2,N-3);
Ha=abs(ts)>pts;
end
