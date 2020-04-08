% function np=interpoint(y,de,tau,bins,nref,nt)
%
% calculate the number of interpoint distances of the time series y
% embedded in dimension de (can be a vector) with lag tau (must be
% a scalar) in the bins specified by bins. 
% I.e.
%  np(k) is *approximately* the number of pairs of points x(:,i)
%  and x(:,j) (with |j-i|>nt) such that
%                bins(k-1)<||x(:,i)-x(:,j)||<bins(k)           
%  [assuming bins(0)=0]
%
% If nref>0 then do the calculation only for nref randomly chosen 
% *pairs* of points (x(:,i),x(:,j)).
% If nref==0 then do the calculation for all the embedded data
%
% de vector of embedding dimensions
% embedding lag
% nt is the Theiler window (default nt=0)
% nref is the number of (randomly chosen) pairs (x(:,i),x(:,j))
% bins are the k bin sizes (must be ascending)
% y is a 1-by-n vector: the (unembedded) time series
%
% NOTE: because all the computations are done at once, we save some
% time, but the length of the data available for lower dimensional
% embeddings is reduced (as a matter of convenience) to be the same
% as the highest dimensional embedding.
% NOTE: If nt=0 then exclude nothing (include i=j), if nt=1 exclude
% *only* i=j. Hence, if nt=0 and nref=0 the lowest bin will always have
% at least n-de+1 occupants!
%
% see also interbin and interbinref
%
%
% For more info, read README
%
% Michael Small
% ensmall@polyu.edu.hk
% 30/4/02
