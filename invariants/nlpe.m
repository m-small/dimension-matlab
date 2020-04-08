function e=nlpe(y,de,tau);

%function e=nlpe(y,de,tau);
%function e=nlpe(y,v);
%
%compute the normalised "drop-one-out" constant interpolation nonlinear
%prediction error for embedding dimension de and lag tau or for embedding
%strategy v (v>0)
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk
%

if nargin<3,
    tau=1;
end;
if nargin<2
    de=3;
end;

%normalised.
y=y-mean(y(:));
y=y/std(y(:));

if min(size(y))>1,
    x=y;
    y=y(1,2:end);
    x=x(:,1:(end-1));
elseif max(size(de))>1,
    v=de(de>0);
    [x,y]=embed(y,v-1);
else,
    [x,y]=embed(y,[-1 0:tau:((de-1)*tau)]);
end;

[de,n]=size(x);



dd=zeros(n,n);
for i=1:de, %loop on de and compute the distance.^2
    dd=dd+(ones(n,1)*x(i,:)-x(i,:)'*ones(1,n)).^2;
end;
% dd is the distnace .^2 with inf on the diag
warning off MATLAB:divideByZero
dd=dd+1./(1-eye(n,n));
warning on MATLAB:divideByZero
%near is the index of the nearest neighbour of each point
[dist,near]=min(dd);
%the prediction error is
e=y(near)-y;
e=mean(e.^2);
    


