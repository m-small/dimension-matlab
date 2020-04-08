% function np=interbinref(x,bins,nref,nt)
%
% calculate the number of interpoint distances of the embedded time
% series in the bins specified by bins. 
% I.e.
%  np(k) is *approximately* the number of pairs of points x(:,i)
%  and x(:,j) (with |j-i|>nt) such that
%                bins(k-1)<||x(:,i)-x(:,j)||<bins(k)           
%  [assuming bins(0)=0]
%
% Unlike interbin, this version does not do the calculation for all
% pairs, but only for nref randomly chosen points x(:,i).
%
% nt is the Theiler window (default nt=0)
% nref is the number of (randomly chosen) reference vectors x(:,i)
% bins are the k bin sizes (must be ascending)
% x is a d-by-n matrix: the d-dimensional embedding of n-points
%
% see interbin
%
%
% For more info, read README
%
% Michael Small
% ensmall@polyu.edu.hk
% 26/2/02
