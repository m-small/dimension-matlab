function [Kkdl]=lfnn(y,de,tau,nb,beta)
  
%function [dl,nfnn]=fnn(y,de,tau,nb,beta)
%
%determine the number of local false nearest neighbours for the time
%series y embedded in local embedding dimension dl<de with lag tau. 
%
%for each pair of values (de,tau) the data y is embeded and the
%nearest neighbour to each point (excluding the immediate
%neighbourhood of n points) is determined. If the ratio of the
%distance of the next (kth) points and these points is greater than
%th then they are counted as false nearest neighbours.
%
% default:
% th=5
% kth=1
%
% p(i,j) is the proportion of false nearest neighbours for de(i)
% and tau(j).
%
%Michael Small
%2/28/2012
%michael.small@uwa.edu.au

ra=mean(abs(y-mean(y))); %attractor size

if nargin<5,
    beta=0.1;%fraction of attractor size
        disp(['beta = ',int2str(nb)]);

end;
if nargin<4,
    nb=40; %nbunmber of neighbours to choose
    disp(['nb = ',int2str(nb)]);
end;
if nargin<3,
  tau=firstzero(y);
  disp(['tau = ',int2str(tau)]);
end;

if nargin<2,
  de=10;
  disp(['de = ',int2str(de)]);
end;

Kkdl=[];
p=[];
  px=[];
  %embed the data
  X=embed(y,de,tau);
  [dx,nx]=size(X);
  hand=waitbar(0,'lfnn: working');
  for i=1:nx,
      dists=rms(X'-ones(nx,1)*X(:,i)');
      dists(i)=inf;
      [dummy,ind]=sort(dists); %ind are the nb closest neighbours in R^de
      ind=ind(1:nb);
      
      %PCA needs u and s^2 
      [u,s,v]=svd(X(:,ind),'econ'); %'econ' since we don't need v
      for dl=1:de,
            %in PCA dl-dim space
            normit=u(:,1:dl)'*(X(:,ind)-X(:,i)*ones(1,nb)); % u or u' - that is the question
            if dl==1,
                dists=abs(normit);
            else
                dists=rms(normit);
            end;

            [dummy,ind1]=min(dists); %the nearest neighbour in dl projection
            indnn=ind(ind1);
            %now we go back to de-space
            k=nx-max(i,indnn)-1;%how far from the end are we
            biggy=abs(y(i+(de-1)*tau+1:i+(de-1)*tau+k)-y(indnn+(de-1)*tau+1:indnn+(de-1)*tau+k));
            delta=find(biggy>beta*ra);
            if isempty(delta),
                 Kkdl(i,dl)=nan;
            else
                Kkdl(i,dl)=delta(1);
            end;
      end;
      waitbar(i/nx,hand);
  end;  
  close(hand);
  
