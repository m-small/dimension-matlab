function fit=judd_da(xi,eps,ci);

%function fit=gka_dks(xi,eps,ci);
%
%return RMS error of the fit of 
%  eps^d*p(eps) and
%    p(eps)=a(1)+a(2)*eps+a(3)*eps^2+ ... +a(n)*eps^(n-1)
%where
% xi(1)=d    : correlation dimension
% xi(2)=a(1) : polynomial coeeficient terms
% xi(2)=a(2) :     "          "         "
% xi(4)=a(3) :     "          "         "
%  ...
% xi(n)=a(n-1):    "          "         "
% eps        : epsilon (bandwidth/viewing scale)
% ci(h)=ci   : Correlation integral
%
%no check of args is done (to speed up calculation)
%
% For more info, read README 
%
% Michael Small
% ensmall@polyu.edu.hk
% 28/2/02

d=xi(1);
a=xi(2:end);
a=a(:)';
pt=length(a);
ne=length(eps);

et=ones(1,ne);
for i=2:pt;
  et=[et; eps.^(i-1)];
end;

fit=(eps.^d).*(a*et)-ci;
fit=sum(fit.^2); %Euclidean norm

