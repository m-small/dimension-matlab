function rhomax=findrho(y,de,tau,target,Aopt);

%function rho=findrho(y,de,tau,ptarget,Aopt);
%
%find the optimal value of rho for the pls surrogate generation algorithm. 
%Uses the premise that the optimal value is that which produces surrogates 
%that have the greatest number of short sequences identical to the data.
%Now the number of short identical sequences, is related to the transition
%probability (i.e. the probability that we don't just follow the current
%state), if the transition probability p is 0.5 the the probability of a
%sequence of length 2 (the shortest allowable) p(1-p) is maximal. Lower the
%probability for longer sequences.
%
%if de is a vector then it is a vector of lags (min 0).
%
%computation is pseudo-analytic (certainly not stochastic)
%
%Aopt is optional
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk

%parameters


if nargin<5,
    Aopt=[];
end;
if nargin<4
    target=0.5;
end;
if nargin<3,
  tau=[];
end;
if nargin<2,
  de=[];
end;
if isempty(de),
  de=3;
end;
if isempty(tau),
  tau=1;
end;

y=y(:);
ny=length(y);
if length(de)>1,
    x=embed(y,de);
else,
    x=embed(y,de,tau);
end;
[de,n]=size(x);
if isempty(Aopt),
    Aopt=ones(1,de);
end;

%memory limitation 
maxsize=3000;
if n>maxsize,
    x=x(:,end+(1-maxsize:0));
    n=maxsize;
    disp('WARNING: too much data to handle in findrho3 ... truncating');
end;

%compute the switch probabilities for each point
%first, the L2-norm^2
dd=zeros(n,n);
for i=1:de, %loop on de and compute the distance.^2
    dd=dd+Aopt(i)*(ones(n,1)*x(i,:)-x(i,:)'*ones(1,n)).^2;
end;
dd=sqrt(dd);


tol=0.0000001;
rhou=std(y);
rhol=rhou/2;
pl=1;pm=1;
        pp=exp(-0.5.*dd/rhou);
        pu=sum(pp);
        pu=(pu-diag(pp)')./pu;
        pu=mean(pu);   
        pp=exp(-0.5.*dd/rhol);
        pl=sum(pp);
        pl=(pl-diag(pp)')./pl;
        pl=mean(pl);

while (rhou-rhol)>tol,
    if pl>target %both bounds too big
        pu=pl;
        rhou=rhol;
        rhol=rhol/2;
        %recompute pl
        pp=exp(-0.5.*dd/rhol);
        pl=sum(pp);
        pl=(pl-diag(pp)')./pl;
        pl=mean(pl);
        rhom=mean([rhol,rhou]);
    elseif pu<target, %or too small
        pl=pu;
        rhol=rhou;
        rhou=rhou*2;        
        %recompute pu
        pp=exp(-0.5.*dd/rhou);
        pu=sum(pp);
        pu=(pu-diag(pp)')./pu;
        pu=mean(pu);
        rhom=mean([rhol,rhou]);
    else, %dead on, so tighten them
        rhom=mean([rhol,rhou]);
        %compute pm
        pp=exp(-0.5.*dd/rhom);
        pm=sum(pp);
        pm=(pm-diag(pp)')./pm;
        pm=mean(pm);
        if pm<target,
            rhol=rhom;
            pl=pm;
        else,
            rhou=rhom;
            pu=pm;
        end;
    end;
    disp([num2str(rhol),'<',num2str(rhom),'<',num2str(rhou)]);
end;
disp(['Final prob=',num2str(pm)]);
rhomax=rhom;
