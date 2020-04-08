function [m,d,k,s,gki]=gka(y,de,tau,nbins,nt,pretty);

% function [m,d,k,s,gki]=gka(y,de,tau,nbins,nt);
%
% compute correlation dimension (d) and entropy (k) and noise level
% (s) using the GKA
%
% de range of embedding dimensions (default 2:20)
% tau embedding lag (default 1)
% nbins number of bins of interpoint distances (default 200)
% nt is the number of temporal neighbours to excluse (default 0)
%
% For more info, read README
%
% Michael Small
% ensmall@polyu.edu.hk
% 28/2/02

nout=nargout;

if nargin<6,
    pretty=[];
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

if isempty(nt), nt=0; end;
if isempty(nbins), nbins=200; end;

nt=max(nt,1); %nt>=1
nde=length(de);

%parameters
maxn=5000; %maximum number of points to use
hcmin=0.05;%0.25; %absolute minimum bandwidth
hcmax=3;   %absolute maximum bandwidth
if isempty(pretty),
    pretty=0;  %pictures?
end;

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
gkim=[];
ssm=[];

%get bins : distributed logarithmically
binl=log(min(diff(unique(y))))-1;   %smallest diff 
%if isempty(binl),binl=eps;end;  %just in case the data is crap
binh=log(max(de)*(max(y)-min(y)))+1;%seems to work
%if isnan(binh),binh=log(max(de)*max(y))+1;end; %just in case the data is crap
binstep=(binh-binl)./(nbins-1);
bins=binl:binstep:binh;
bins=exp(bins);

%disp
disp(['GKA (n=',int2str(n),'; tau=',int2str(tau),'; nbins=',int2str(nbins),'; nt=',int2str(nt),')']);
disp(['hcmin=',num2str(hcmin),' & hcmax=',num2str(hcmax)]);
disp('Computing histogram');

%estimate distribution of interpoint distances 
if n>2*maxn, %why sample with replacement when you could without?
  disp(['Using ',int2str(maxn),' reference points'])
  %distribution of interpoint distances
  %compute distrib. from maxn ref. points
  np=interpoint(y,de,tau,bins(1:(end-1)),maxn^2,nt);  
  %number of interpoint distances      
  ntot=maxn^2;                    
else,
  disp('Using all points');
  %distribution of interpoint distances
  %compute distrib. using all points
  np=interpoint(y,de,tau,bins(1:(end-1)),0,nt);
  %number of interpoint distances      
  ntot=n-(de(end)-1)*tau;
  if nt>0,
    ntot=(ntot-2*nt+2)*(ntot-2*nt+1)+2*(nt-1)*(ntot-nt+1)-(nt-1)*nt;
  else
    ntot=ntot^2;
  end;
end;

disp('First pass:');
disp(sprintf(' m\t D\t s\t B')); 

%set the first guess here, next first guess can then be last final
xi=[de(1) 0.1 1];

%loop on de
for mi=1:length(de),
  m=de(mi);

    %bandwidth is computed directly from bins
    bands=mean(embed([0 bins],2,1));
    
    %compute Gaussian kernel correlation integral
    npt=np(:,mi);%./ntot;
    for i=1:nbins,
        gki(i)=sum(npt'.*exp(-(0.5*bins./bands(i)).^2))./ntot;
	    ss(i)=sum((npt'.*exp(-(0.5*bins./bands(i)).^2)./ntot).^2) ./(nbins-1) ...
	       - (gki(i)./(nbins-1)).^2;
    end;
    ss=sqrt(abs(ss));  %standard deviations of bin sizes
    ss=ss./sum(ss);
    
    %fit to find D and s
    hc=0.8;     % the upper cut off of the linear scaling region: This 
                % param is CRUCIAL. If things are going wrong, then
		        % this is probably a good place to start looking
    ind=find(bands<hc);
    opt=optimset('TolX',1e-6,'TolFun',1e-6,'display','notify');
    xi=fminsearch('gka_dsb',xi,opt,m,bands(ind),gki(ind),ss(ind)); % do it once (to est. s)
    hc=3*xi(2);                     % set upper cut off at three times noise
    hc=min(max(hc,hcmin),hcmax);    % set hcmin<=hc<=hcmax
    ind=find(bands<hc);
    xi=fminsearch('gka_dsb',xi,opt,m,bands(ind),gki(ind),ss(ind)); % do it again
    d=[d xi(1)];
    s=[s xi(2)];
    b=[b xi(3)];
    
    %display is necessary
    if pretty,
      figure(gcf);
      clf;
      subplot(211);
      loglog(bands,gki,'k:');hold on;
      loglog(bands(ind),gki(ind),'r');
      loglog(bands(ind),gki(ind)-ss(ind),'r:');
      loglog(bands(ind),gki(ind)+ss(ind),'r:');
      gfit=gkifit(xi(1),xi(2),xi(3),m,bands);
      loglog(bands,gfit,'g-');
      axis([bands(1) bands(end) max(min(gki(gki>0)),min(gfit)) 1]);
      grid on;
      xlabel('log(\epsilon)');
      ylabel('log(T_m(\epsilon))');
      title(['Gaussian Kernel Correlation Integral (m=',int2str(m), ...
		    ' and tau=',int2str(tau),')']);
      subplot(212);
      errorbar(bands(ind),gki(ind),ss(ind),'r');hold on;
      plot(bands(ind),gfit(ind),'g-');
      plot(bands(ind),abs(gfit(ind)-gki(ind)),'b');
      xlabel('\epsilon');
      ylabel('T_m(\epsilon)');
      drawnow;
    end;
    
    %display
    disp(sprintf(' %d\t %0.3f\t %0.3f\t %0.3f',m,xi));

    %remember the gki and ss
    gkim=[gkim;gki];
    ssm=[ssm;ss];
end;

%now fit to find K
if nde>1, %need more than one embedding dimension
  
  disp('Computing phi:');
  %find K and then phi
  %note: b(m)=phi.*exp(-m*K*tau), so b(m+1)/b(m)=exp(-K*tau) and 
  % K = log(b(m)/b(m+1))/tau
  %and this is what we do as an initial guess
  %
  %First, make a guess for K and phi
  k=b(1:(nde-1))./b(2:nde);
  k(k<=0)=eps; %ensure non-negativity
  k=log(k)/tau;  %these are the first guess(es) at K
  phi=b(1:(nde-1)) .* exp(de(1:nde-1).*k.*tau);
  phi=mean(phi); %initial guess of phi
  k1=mean(k);    %initial guess of k 
  %Second, do nonlinear fit
  opt=optimset('TolX',1e-6,'TolFun',1e-6,'display','notify');%...
	   %   'LevenbergMarquardt','on');
  xi=[k1 phi];
  xi=fminsearch('gka_kphi',xi,opt,b,de,tau);
  k1=xi(1);
  phi=xi(2);
  disp(['  K=',num2str(k1),' and phi=',num2str(phi)]);
  
  %Just to iron out any wrinkles,
  %fit d,s and k all over again 
  disp('Final fit:');
  disp(sprintf(' m  \t D  \t K  \t s  \t e (n) \t   hc'));
  for i=1:nde,
    hc=3*s(i);                      % set upper cut off at three times noise
    ind=find(bands<hc); 
    hc=min(max(hc,hcmin),hcmax);    % set hcmin<=hc<=hcmax
    ind=find(bands<hc);
    opt=optimset('TolX',1e-6,'TolFun',1e-6,'display','notify');
    xi=[d(i) k1 s(i)];
    gki=gkim(i,:);
    xi=fminsearch('gka_dks',xi,opt,phi,de(i),tau,bands(ind),gki(ind),ss(ind)); 
    rms=gka_dks(xi,phi,de(i),tau,bands(ind),gki(ind),ss(ind)); 
    disp(sprintf(' %d\t %0.3f\t %0.3f\t %0.4f\t %0.2g (%d)\t %0.3f',de(i),xi,rms,length(ind),hc));
    d(i)=xi(1);
    k(i)=xi(2);
    s(i)=xi(3);
    
    %display is necessary
    if pretty,
      figure(gcf);
      clf;
      subplot(211);
      loglog(bands,gki,'k:');hold on;
      loglog(bands(ind),gki(ind),'r');
      loglog(bands(ind),gki(ind)-ss(ind),'r:');
      loglog(bands(ind),gki(ind)+ss(ind),'r:');
      gfit=gkifit(xi(1),xi(3),phi*exp(-de(i)*xi(2)*tau),de(i),bands);
      loglog(bands,gfit,'g-');
      axis([bands(1) bands(end) max(min(gki(gki>0)),min(gfit)) 1]);
      grid on;
      xlabel('log(\epsilon)');
      ylabel('log(T_m(\epsilon))');
      title(['Gaussian Kernel Correlation Integral (m=',int2str(de(i)), ...
		    ' and tau=',int2str(tau),')']);
      subplot(212);
      errorbar(bands(ind),gki(ind),ss(ind),'r');hold on;
      plot(bands(ind),gfit(ind),'g-');
      plot(bands(ind),abs(gfit(ind)-gki(ind)),'b');
      xlabel('\epsilon');
      ylabel('T_m(\epsilon)');
      drawnow;
    end;

  end;
  
else,
  
  k=b;
  disp('WARNING: cannot compute entropy from a single embedding dimension');
  
end;

m=de;

%output display?
if pretty,
  figure(gcf);
  clf;
  subplot(311);
  plot(m,d,'r');
  xlabel('embedding dimension, m');ylabel('correlation dimension d');
  subplot(312);
  plot(m,k);
  xlabel('embedding dimension, m');ylabel('entropy, k');
  subplot(313);
  plot(m,s);
  xlabel('embedding dimension, m');ylabel('noise level, s');
end;

  
%check for additional output arguments
if nout>4,
  gki=[bands; gki; ss;];
end;

