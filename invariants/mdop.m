function [taus]=mdop(y,taumax)

%function taus=mdop(y,taumax)
%
%My implementation of the MDOP (maximising derivatives on projection)
%method of Nichkawde (PRE 87, 022905 (2013)) to find an optimal non-uniform
%embedding
%
%


y=y(:);
if nargin<2,
    taumax=[];
end;
ly=length(y);

if isempty(taumax),
    taumax=floor(ly/30);
end;

mpay=inf; %positive marginal payoff in the objective funtion
d=1;
taus=0;
w=taumax+1;
lastphi=inf;
figure; hold on;
while mpay>0 %& length(taus)==length(unique(taus)),
    %embed with current best strategy
    x=embed(y,taus);
    %find nearest neighbours - saving space at the begining for new lags
    ind=nearestk(x(:,w:end),ceil(taumax/2));
    ind=ind+taumax; 
    %demoninator
    denom=myrms(x(:,w:end)-x(:,ind));
    denom=ones(w,1)*denom;
    %numerator
    numer=[];
    for i=-taumax:0,
        numer=[abs(x(1,i+(w:end))-x(1,i+ind)); numer];
    end;
    %average
    logbeta = mean(log( (numer./denom)' ));
    plot(logbeta);
    [thisphi,tau]=max(logbeta);
    semilogy(tau,thisphi,'rp');
    taus=[taus tau];
    mpay= -(thisphi-lastphi);
    lastphi=thisphi;
    drawnow;
    title(['cycle ',int2str(length(taus))]);
end;
hold off;
    
if length(taus)~=length(unique(taus)),
    disp('MDOP: Duplicate lags');
end;
    

function d=myrms(x)

[mx,nx]=size(x);
if mx==1,
    d=sqrt(x.^2);
else,
    d=sqrt(sum(x.^2));
end;

