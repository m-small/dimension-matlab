function [yp,np,ord]=localoptimiser(X,yi,de,nx,npmax,ordmax);

%function [yp,np,ord]=localoptimiser(X,yi,de,nx,npmax,ordmax);
%
%
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk
%

%find nearest neighbours of y(i)
rr=rms((X-yi*ones(1,nx))')'; %RMS distance
rr(end)=[];                  %always exclude the last point (it doesn't
                             %have an image). 
[rdist,rind]=sort(rr);       %X(:,rind(i)) is the i-th closest point to
                             %yi, it is rdist=rr(rind) distance away
                             
orders=0:ordmax;
%for oi=1:length(orders), %loop through model order o
oi=1;
  del=[];
  o=orders(oi);
  neighs=2:npmax;
  for ni=1:length(neighs); %loop through number of neighbours n
    
    %build mode of this size
    n=neighs(ni);
    lam=X(:,rind(1:n)+1);
    if n>1,
       lam=mean(lam')';
    end;
    
    yp=X(:,rind(1:n)+1);
    mss=yp-lam*ones(1,n);    %prediction errors
    mss=mean(sum(mss.^2)); %mean sum of squares
    
    %and compute MDL
%    Q=2*eye(de);
%    del=abs(lam)*0.01; % first guess precisions
%    del=find_deltas(Q,del);keyboard;
    del=sqrt(2)*ones(de,1);

    mdl(ni)=description_length(mss,lam,del,n);
%    aka(ni)=n*log(mss)+2;
%    sch(ni)=n*log(mss)+log(n);
    
  end;
%end;

%find best MDL model
[minmdl,ni]=min(mdl./neighs);
np=neighs(ni);
%plot(mdl./neighs);hold on;
%plot(np,mdl(ni)./np,'rp');hold off;drawnow;
yp=mean(X(:,rind(1:np)+1)')';
disp(['Optimal # neigbours = ',int2str(np)]);
ord=nan;


