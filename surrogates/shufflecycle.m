function [zz,splits,pz] = ...
    shufflecycle(z,p,pt,t,tt,where,shift)

% [z,splits]=shufflecycle(z,where,shift)
%    or shufflecycle(z,p,pt,t,tt,where,shift)          shuffle the cycles of z
%
% seperates the cycles of z at peak inspiration and (randomly) reorders
% them - aka Theiler.
% where specifies where to chop up the timeseries, either;
%  where='peak'
%  where='trough' or
%  where='centre' (mean, going up)
% 
% shift is used to avoid a possible problem which occurs when sampling 
% frequency and underlying system frequency are incommesurate (i.e. always, 
% p=1). If shift=1 then adjacent (shuffled) cycles are shifted vertically to 
% fit. If shift=0 then no vertical shifting is performed (this avoids
% incommesurate issues but may introduce aparent discontinuities).
%
% default='peak', shift=1
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk

% Copyright (c) 1997 by Michael Small.
%
% Please see the copyright notice included in this distribution
% for full details.
%
%
% File   shufflephase.m
%   $Id$
%
% Created by Michael Small (<watchman@>) on Wed May 14 1997
%
% $Log$

na=nargin;
if na<7,
  shift=[];
end;
if na<6
  where=[];
end;
if na<3,
  pt=[];
end;
if na<2,
  p=[];
end;
if na<4
  where=p;
  shift=pt;
end;
if isempty(where),
  where='peak';
end;
if isempty(shift),
 shift=1;
end;
  
[lz,b]=size(z);
if b~=1; z=z'; [lz,b]=size(z); end;

if na<5 & ~strcmp(where,'zeros');
  [p,pt,t,tt]=peaktrough([z; z(1:min(length(z),500))]);
                     % need to add a bit to make it the same length
  ind=find(pt<lz & tt<lz);
  p=p(ind);pt=pt(ind);t=t(ind);tt=tt(ind);
end;



if strcmp(where,'peak'),
  lt=diff(pt);
  tm=pt;
elseif strcmp(where,'trough')
  lt=diff(tt);
  tm=tt;
elseif strcmp(where,'centre');
  zz=z(1:lz-1)-mean(z); nz=z(2:lz)-mean(z);
  tm=find(zz<0 & nz>0);
  lt=diff(tm);
  clear zz nz
else
  disp('invalid splitting method');
end;

splits=tm;
n=length(tm)-1;
start=z(1:tm(1));
finish=z(tm(n+1)+1:lz);

if ~isempty(start),
  zz=start;
else
  zz=mean(p);
end;
pz=zz;
while n>0
  ind=floor(rand(1,1)*n)+1;
  vert=zz(length(zz));
  pts=(tm(ind)+1):(tm(ind)+lt(ind));
  pz=[pz; z(pts);];
  zz=[zz; z(pts)-(z(tm(ind))-vert)*shift;];
  tm(ind)=[];
  lt(ind)=[];
  n=n-1;
end;

vert=zz(length(zz))*shift;
if ~isempty(finish),
  zz=[zz; finish-(finish(1)-vert)*shift;];
end;
z=zz;

% End of shufflephase.m
