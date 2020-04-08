function [m,dcm,eps0,cim,bins]=judd(y,de,tau,nbins,nt,dt);

% function [m,dc,eps0,ci,bins] = judd(y,de,tau,nbins,nt,dt);
%
% compute correlation dimension (dc) using Judd's algorithm
%
% de range of embedding dimensions (default 2:20)
% tau embedding lag (default 1)
% nbins number of bins of interpoint distances (default 200)
% nbins may be a vector or a scalar if it is a vector, then it is assumed
% to be the bin locations.
% nt is the number of temporal neighbours to excluse (default 0)
% dt is the topological dimension of the set (default, 2)
%
% dc is a cell array of the same size as m (the list of embedding
% dimensions). For embedding in m(i) dimensions dc{i}(1,:) are the eps0
% values, dc{i}(2,:) is the correlation dimension estimates (for the
% corresponding eps0) and dc{i}(3,:) is the estimated fitting error.
%
% Michael Small
% ensmall@polyu.edu.hk
% 25/4/02

nout=nargout;

if nargin<6,
    dt=2;
end;
if nargin<5,
    nt=0;
end;
if nargin<4,
    nbins=200;
end;
if nargin<3,
    tau=1;
end;
if nargin<2,
    de=2:20;
end;

if dt<1,
    dt=1; %dt can't be less than one
end;

nde=length(de);

%parameters
maxn=5000; %maximum number of points to use
pretty=0;  %pictures?
noccup=20; % minimum number of occupied bins to fit to.
maxci=0.9; % upperbound on correlation integral
errorbound=0.01; %maximum fitting error to blindly accept

%data
y=y(:);
n=length(y);

%rescale to mean=0 & std=1
y=y-mean(y);
y=y./std(y);

%init
m=[];
d=[];
k=[];
s=[];
b=[];
cim=[];

%get bins : distributed logarithmically
if max(size(nbins))==1,
    binl=1+log(min(diff(unique(y))));   %smallest diff
    binh=log(max(de)*(max(y)-min(y)));%seems to work
    binstep=(binh-binl)./(nbins-1);
    bins=binl:binstep:binh;
    bins=exp(bins);
else
    bins=sort(nbins);
    nbins=length(bins);
end;

%disp
disp(['Judd''s Algorithm (n=',int2str(n),'; tau=',int2str(tau),'; nbins=',int2str(nbins),'; nt=',int2str(nt),'; dt=',int2str(dt),')']);

%get distributions of interpoint distances
if n>2*maxn, %why sample with replacement when you could without?
    %distribution of interpoint distances
    %compute distrib. from maxn ref. pairs of points
    np=interpoint(y,de,tau,bins(1:(end-1)),maxn.^2,nt);
    %number of interpoint distances
    ntot=maxn.^2;
   disp(['Using ',int2str(maxn),'^2 reference points (ntot=',int2str(ntot),')']);
else,
    %distribution of interpoint distances
    %compute distrib. using all points
    np=interpoint(y,de,tau,bins(1:(end-1)),nt);
    %number of interpoint distances
    nx=length(y)-(max(de)-1)*tau;
    ntot=nx*(nx-(1+2*nt));
    disp(['Using all points (ntot=',int2str(ntot),')']);
end;

%loop on de
for mi=1:nde,
    
    disp(['Fitting for m=',int2str(de(mi))]);
    
    %compute correlation integral
    ci=cumsum(np(:,mi)'./ntot);
    ind=find(ci<maxci);
    
    dc=[];
    eps0=[];
    errs=[];
    while(sum(diff(ci(ind))>0)>noccup), %keep going so long as noccup
        %bins are occupied
        %fit to find D and a
        opt=optimset('TolX',1e-6,'TolFun',1e-6,'display','notify',...
            'MaxFunEvals',10000,...
            'MaxIter',10000);	%   'LevenbergMarquardt','on',
        xi=[de(mi) 1]; %initial guess
        xi=fminsearch('judd_da',xi,opt,bins(ind),ci(ind)); % do it once (to est. d)
        xi=[xi(1) xi(2) zeros(1,dt-1)];
        [xi,gfit,exitflag]=fminsearch('judd_da',xi,opt,bins(ind),ci(ind)); % do it again
        
        %compute normalised error
        gfit=gfit./sum(ci(ind).^2);
        
        %should we give up?
        if ~exitflag,
            %fitting didn't converge .. so quit
            disp('WARNING: Fitting failed to converge.');
            break;
        end;
        
        %was it good enough?
        if gfit<errorbound & xi(1)>0,
            dc=[dc xi(1)];
            eps0=[eps0 bins(ind(end))];
            errs=[errs gfit];
        end;
        
        %display is necessary
        if pretty,
            figure(gcf);
            clf;
            subplot(211);
            loglog(bins(ind),ci(ind),'k:');hold on;
            loglog(bins(ind),ci(ind),'r');
            tfit=judd_fit(xi(1),xi(2:end),bins(ind));
            loglog(bins(ind),tfit,'g-');
            axis([bins(1) bins(end) max(min(ci(ci>0)),min(tfit)) 1]);
            grid on;
            xlabel('log(\epsilon)');
            ylabel('log(P_\mu(\epsilon))');
            title(['Correlation Integral (m=',int2str(de(mi)), ...
                ' and tau=',int2str(tau),')']);
            subplot(212);
            plot(bins(ind),tfit(ind),'g-');hold on;
            plot(bins(ind),ci(ind),'r');
            plot(bins(ind),abs(tfit(ind)-ci(ind)),'b');
            title(['\epsilon_0=',num2str(bins(ind(end))),', d_c=', ...
                num2str(xi(1)),', gfit=',num2str(gfit)]);
            xlabel('\epsilon');
            ylabel('P_\mu(\epsilon)');
            drawnow;
        end;
        
        ind(end)=[];
    end;
    
    %normalisation factor for the bins
    bbox=(max(y)-min(y))*sqrt(de(mi));%diag of bounding box in R^m
    
    %remember the gki and ss
    cim=[cim;ci];
    dcm{mi}=[eps0./bbox;dc;errs];
end;


m=de;

%output display?
if pretty,
    figure(gcf);
    clf;hold on;
    for i=1:nde,
        if min(size(dcm{i}))>1,
            semilogx(dcm{i}(1,:),dcm{i}(2,:));
        end;
    end;
    xlabel('log(\epsilon_0)');
    ylabel('d_c(\epsilon_0)');
end;



