% xs=iAAFTsurr_HRV(x,zones,Fs, specflag);
%
% Inputs:
% - x = input signal
% - f1 and f2 = This algorithm preserves the fourier phases between those
%               frequencies
% - Fs = Sample frequency of the observed signal. 
% -specflag: 1-> exact amplitude spectrum, otherise exact signal
%            distribution
%
% Outputs:
% Xs = iAAFT surrogate generated by randomizing the phases only in the 
%      selected frequency band.
% c = number of itetarions before convergence	

% Included in this ditribution with permission of the author:
%
% Diego Guar�n
% dlguarin@gmail.com
% 26/08/2011 

function [Xs,c] = BPRsurrogates (X,f1,f2,Fs,specflag)



if (nargin<1)
	Xs = [];
	return;
end
%checking the inputs
n=length(X);
if (nargin<5)
    specflag = 1;
end
if (nargin<4)
    Fs=1;
    specflag = 1;
end
if (nargin<3)
    Fs=1;
    f2=0.5;
    f1=0;
	specflag = 1;
end
if (nargin<2)
    Fs=1;
    f2=0.5;
    f1=0;
	specflag = 1;
end
if f1>f2
    f1=0;
    f2=0.5;
end

if (isempty(X))
	return;
end

max_it = 1000;
N = length(X);

% Initial Conditions
%rn=generate_iAAFT(X,specflag);
rn=surrogates(X,'alg3',1);
%rn_c=rn;
FX=fft(X);
Xsorted = sort(X);	% Desired signal distribution
Yamp = abs(FX);	% Desired amplitude spectrum
Yang = angle(FX); % Desired phases

if rem(N,2)
    Yang=Yang(1:(end+1)/2);
else
    Yang=Yang(1:(end/2)+1);
end



% Randomization of the desired phases
f=(1:length(Yang))*Fs/length(X);

zone= f>f1 & f<f2;

c = 1;
err = 1;

sn_comp=[];

while (c<max_it && err>1e-50)
	% Match Amplitude Spec
    phase=angle(fft(rn));
    if rem(N,2)
        phase=phase(1:(end+1)/2);
    else
        phase=phase(1:(end/2)+1);
    end
    nozone=1:length(phase);
    nozone(zone)=[];
    phase(nozone)=Yang(nozone);
    
    if rem(N,2)
        phase=[phase ; -phase(end:-1:2)];
    else
        phase=[phase ; -phase(end-1:-1:2)];
    end
    
    sn = real(ifft(Yamp.*exp(sqrt(-1).*phase)));
	% Scale to Original Signal Distribution
	[sns INDs] = sort(sn);
	rn(INDs) = Xsorted;

    %seeking for the stop criterium
    if ~isempty(sn_comp)
        % actual sn vs the previous sn
        err=sqrt((1/n)*sum((sn-sn_comp).^2));
    end
    if c>10
        %a copy of the actual sn
        sn_comp=sn;
    end

	c = c+1;
end

if c==max_it
    error('Convergence problem!!!')
end

if (specflag==1)
	Xs = sn;	% Exact Amplitude Spectrum
else
	Xs = rn;	% Exact Signal Distribution
end