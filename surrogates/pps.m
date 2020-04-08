
% function [yp,yi]=pps(y,de,tau,n,rad,y1);
%                 =pps(y,emb,[],n,rad,[]);
%
% n-point PPS (pseudo-periodic surrogate) of y, embedding with de and
% tau. Randomised with radius rad.
% yp is the surrogate, yi is the indices of the embedded version of 
% y selected for the surrogate
% emb is a (vector embedding strategy.
% emb=[0:tau:(de-1)*tau] is equivalent to providing de and tau
%
% y1 is the (optional) index of the first datum in the surrogate.
%
% Implmenented in PPS.C
%
% MS 6/1/00, 13/7/04
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk




