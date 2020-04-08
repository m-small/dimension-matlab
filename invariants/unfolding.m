function [d,nfnn]=unfolding(y,th,de,tau);
  
%function [de,nfnn]=unfolding(y,th,de,tau)
%
%estimate the minimum unfolding dimension by calculating when the
%proportion of false nearest neighbours if first below th.
%
% The number of false nearest neighbours are calculated for the
% time series y embedded in dimension de with lag tau. 
%
%for each pair of values (de,tau) the data y is embeded and the
%nearest neighbour to each point (excluding the immediate
%neighbourhood of n points) is determined. If the ratio of the
%distance of the next points and these points is greater  
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk
%

if nargin<4,
  tau=firstzero(y);
  disp(['tau = ',int2str(tau)]);
end;

if nargin<3,
  de=[1:10];
  disp(['de = ',int2str(de(1)),':',int2str(de(end))]);
end;

if nargin<2,
  th=0.01;
  disp(['th = ',num2str(th)]);
end;

dsp=5; %separation cutoff
y=y(:);

n=2*tau;%exclude nearpoints

nfnn=[];
locamin=0;

for d=de,
  
  %embed the data
  X=embed(y,d,tau);
  [d,nx]=size(X);

  falseneighbours=0;
  total=0;
    %find the nearest neighbours of each point
    ind=nearest(X(:,1:(nx-1)),tau); %whooh hooo!
    

    %distance between each point and its nearest neighbour
    d0=rms(X(:,(1:(nx-1)))'-X(:,ind)');
    %... and after one time step
    d1=rms(X(:,2:nx)'-X(:,ind+1)');

    %exclude any coincident points
    d1(d0==0)=[];
    d0(d0==0)=[];
    
    %calculate the proportion fnn
    prop=sum((d1./d0)>dsp)/length(d0);


  disp(['de='int2str(d),', n(fnn)=',num2str(prop*100),'%']);
							
  %is data sufficiently unfolded?
  if (prop<th) 
    nfnn=prop;
    return;
  end;

  %or maybe a local min
  if (length(nfnn)>1)
    if (min(prop,nfnn(end))>nfnn(end-1)),
      localmin=1;
      localmini=length(nfnn)-1;
    end;
  end;
  
  nfnn=[nfnn prop];
  

end;

if (localmin),
  %we didn't go subthreshold, but we did have a local min.
  nfnn=nfnn(localmini);
  d=de(i);
else
  %otherwise, just do the best we can
  [nfnn,i]=min(nfnn);
  d=de(i);
end;
