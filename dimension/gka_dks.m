function fit=gka_dks(xi,phi,m,tau,bands,gki,ss);

%function fit=gka_dks(xi,phi,m,tau,bands,gki,ss);
%
%return RMS error of the fit of 
%  phi*exp(-m*k*tau)*(h^2/(h^2+s^2))^(m/2)*((h^2+s^2)/m)^(D/2) to gki(h)
%where
% xi(1)=D    : correlation dimension
% xi(2)=s    : noise level
% xi(3)=b    : coefficient (entropy term times constant)
% m=m        : embedding dimension
% tau=tau    : embedding lag
% phi        : constant
% h=bands    : bandwidths
% gki(h)=gki : Gaussian kernel correlation integral
%
%no check of args is done (to speed up calculation)
%
% For more info, read README 
%
% Michael Small
% ensmall@polyu.edu.hk
% 28/2/02

d=xi(1);
k=xi(2);
s=xi(3);

h2=bands.^2;

fit=(phi*exp(-m*k*tau)) .* [(h2./(h2+s.^2)).^(m/2)] .* [((h2+s.^2)./m).^(d/2)];
fit=gki-fit;
fit=(fit.^2)*ss'; %Euclidean norm


