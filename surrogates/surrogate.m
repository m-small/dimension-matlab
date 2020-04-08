function xp = surrogate(x,alg,rot)

% xp = surrogate(x,alg,rot)
% Surrogate data generator
% handles multichannel data (channels are columns)
% alg = algorithm 1 or 2 (default alg=1)
% rot = method of randomizing phases (default rot=1)
%     +ve = preserve correlation of channels 
%     -ve = destroy correlation of channels 
%       1 = add random phase
%       2 = replace phase
% if alg = -1 or -2, then apply enpoint correction
%
%
% NB: if using one channel then input vector MUST BE a column vector. MAS.
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk

% History
% Surrogate algorithm 2 ala James Theiler
% Written by Tanya Schmah
% Modified by Alistair Mees for Matlab 4
% Modified by Kevin Judd into a function
% Rewritten by Kevin Judd Aug95
% Modified by Michael Small

if nargin<3
  rot= 1;
end
if nargin<2
  alg= 1;
end
if ~any(abs(alg)==[1 2]),error('Only know algorithms 1 and 2');end;

%figure
%subplot(211)
%plot(x)
%drawnow

[n,d] = size(x);
y= zeros(n,d);

if abs(alg)==2
  % normal rank
  j= (1:n)';
  r= randn(size(x));
  [sr,ri]= sort(r);
  [sx,xi]= sort(x);
  [sxi,xii]= sort(xi);
  for k=1:d
    y(:,k)= sr(xii(:,k));
  end
else
  y= x;
end
m= mean(y);
y= y-m(ones(n,1),:);

% endpoint correction
if alg<0
  y0= y(1,:);
  y1= y(n,:);
  u= (0:n-1)'/(n-1);
  c= u*y1 + (1-u)*y0;
else
  c= 0;
end;

% fft
fy = fft(y-c);

% randomize phases (Must be same change each channel)
if rot>0
  phz= rand(n,1);
  phz= phz(:,ones(1,d));
else
  phz= rand(n,d);
end
if abs(rot)==1
  rot= exp(1) .^ (2*pi*sqrt(-1)*phz);
  fyp= fy .* rot;
else
  rot= exp(1) .^ (2*pi*sqrt(-1)*phz);
  fyp= abs(fy) .* rot;
end

%keyboard

% ifft
yp= real(ifft(fyp)) + c + m(ones(n,1),:);

if abs(alg)==2
  [syp,ypi] = sort(yp);
  [sypi,ypii] = sort(ypi);
  for k=1:d
    xp(:,k) = sx(ypii(:,k));
  end;
else
  xp= yp;
end;

%subplot(212)
%plot(xp)

return

%pause
%subplot(211)
%xf=filtfilt(ones(1,10)/10,1,x);  plot(xf);
%subplot(212)
%xpf=filtfilt(ones(1,10)/10,1,xp);plot(xpf);



