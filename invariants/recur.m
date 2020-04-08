function [rx,tol]=recur(y,tol);

%function rx=recur(y,tol);
%
%return the sparse recurrence matrix rx (such that spy(rx) is the
%recurrence plot) for the embedded data y, such that all points within a
%distance tol are marked
%
%if tol is not given, it is set such that each point has atleast one
%recurrence
%
% MS
% 18/5/07

if nargin<2,
    tol=[];
end;
[dy,ny]=size(y);
rx=sparse(ny,ny);

if isempty(tol);
    [d,i]=nearneigh(y,1);
    tol=max(d)+eps;
    clear d i
end;

hand=waitbar(0,'recur, running...');
for i=1:ny,
    d=rms((y-y(:,i)*ones(1,ny))');
    ind=find(d<tol);
    ind(ind==i)=[];
    rx(i,ind)=1;
    waitbar(i/ny,hand);
end;
close(hand);
