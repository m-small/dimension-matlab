%
%   smallshuffle.m
%
%   shuffle surrogate but only move data on average n places
%   i.e. shuffle data but preserve variable volatility.

function y = smallshuffle(y, n);

%   y: data
%   n: amplitude of Gaussian noise added to a rank order
%
%   Copyright (c) 2005 By Michael Small and Tomo Nakamura
%
% NAME
%   smallshuffle.m
%
% HISTORY
%   Michael Small - 15 July 2004 (Thu): Created
%   Tomo Nakamura - 1 Feb 2005 (Tue): Revised
%   Tomo Nakamura - 3 Feb 2005 (Thu): Revised
%-----------------------------------------------------------------

y = y(:);
ly = length(y);
num_yi = (1:ly)';

%	re-ordering
yy = num_yi + n * randn(ly, 1);
[dn, di] = sort(yy);
y = y(di);

