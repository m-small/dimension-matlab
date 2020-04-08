function fit=gka_kphi(xi,b,m,tau);

%function fit=gka_kphi(xi,b,m,tau);
%
%return RMS error of the fit of 
% b(m)=phi*exp(-m*K*tau)
%
%where
% xi(1)=phi  : constant
% xi(2)=K    : entropy
% b          : b(m)=phi*exp(-m*K*tau)
% m=m        : embedding dimension
% tau=tau    : embedding lag
%
%no check of args is done (to speed up calculation)
%
% For more info, read README
%
% Michael Small
% ensmall@polyu.edu.hk
% 28/2/02

k=xi(1);
phi=xi(2);


fit=phi.*exp(-m*k*tau);
fit=b-fit;
fit=sum(fit.^2); %Euclidean norm


