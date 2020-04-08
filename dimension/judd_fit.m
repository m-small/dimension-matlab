function y=judd_fit(d,a,eps);

%function y=judd_fit(d,a,b,eps);
%
% Evaluate the fit to the CI at bands eps
%
% For more info, read README
%
% Michael Small
% ensmall@polyu.edu.hk
% 28/2/02

pt=length(a);
ne=length(eps);

et=ones(1,ne);
for i=2:pt;
  et=[et; eps.^(i-1)];
end;

y=(eps.^d).*(a*et);

