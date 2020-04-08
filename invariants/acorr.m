function [p,t]=acorr(y,rmax);

% function [p,t]=acorr(y,r);
%
% Calculates an approximation to the autocorrelation function for a timeseries
% y.
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk

if nargin<2,
  rmax=length(y);
end;

X=mean(y);
N=length(y);
wait=waitbar(0,'acorr : working...');
for r=0:rmax
  waitbar(r/rmax,wait);
	p(rmax+1+r)=1/(N)*dot(y(1:N-r)-X,y(r+1:N)-X);
	p(rmax+1-r)=p(rmax+1+r);
end;
p=p/(p(rmax+1));
t=-rmax:1:rmax;
close(wait);

