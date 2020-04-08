function s = ...
    alg2v2(y,n)

% z=alg2v2(y,n)  alg2 surrogates ala Schreiber and Schmitz
%
% Iteratively generated algorithim 2 surrogates, the iterative bit is an
% attempt to make the autocorrelation and power spectrum more like the
% data. y is data, z is surrogate and n is the number of iterations of
% alg.
%

% Copyright (c) 1997 by Michael Small.
%
% Please see the copyright notice included in this distribution
% for full details.
%
%
% File   alg2v2.m
%   $Id$
%
% Created by Michael Small (<watchman@>) on Wed May 21 1997
%
% $Log$

if nargin<2
  n=100;
end;
fy=fft(y);
afy=abs(fy);
[ry,iy]=sort(y);
s=shuffle(y);
for i=1:n,
  %fourier transform
  fs=fft(s);
  %fiddle magnitude
  fs=fs./abs(fs).*afy;
  %invert
  s=ifft(fs);
  %re-rank order
  [rs,is]=sort(s);
  s(is)=ry;
end;


% End of alg2v2.m
