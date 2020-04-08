function y=gkifit(d,s,b,m,bands);

%function y=gkifit(d,s,b,m,bands);
%
% Evaluate the fit to the GKI at bands
%
% For more info, read README
%
% Michael Small
% ensmall@polyu.edu.hk
% 28/2/02

h2=bands.^2;

y =b * [(h2./(h2+s.^2)).^(m/2)] .* [((h2+s.^2)./m).^(d/2)];


