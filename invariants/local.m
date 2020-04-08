function [y,np]=local(y,de,tau,n,yi,sig);

%function yp=local(y,de,tau,n,yi,noise);
%           =local(X,[],[],n,yi,noise);
%
%An n point local model simulation of y (embedded in de
%dimensions with a lag of tau, or X a de-by-m matrix starting at yi 
%(either a point or a
%datum number. At each time step, the local image is
%obtained by a de point interpolation of the image of the nearest
%points. The number of neighbouring points and the order of that
%interpolation (constant, linear, quadratic...) is determined using MDL.
%noise is the standard deviation of gaussian dynamic
%noise (default 0).
%
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk

%algorithim parameters
npmax=200;
ordmax=0;

%sort out the parameters
na=nargin;
if na<6,
  sig=0;
  if na<5,
    yi=[];
    if na<4,
      n=[];
      if na<3,
	tau=[];
	if na<2,
	  de=[];
	end;
      end;
    end;
  end;
end;


%sort out embedding parameters
if isempty(tau),
  tau=1;
end;
if isempty(de),
  de=3;
end;

%sort out y ... and embed it (if necessary)
[dy,ny]=size(y);
if dy~=1 & ny ~=1, %y is a matrix don't embed
  if ny<dy,
    y=y';
    [dy,ny]=size(y);
  end;
  de=dy;
  X=y;
else, % y is a vector
  y=y(:);
  ny=max(dy,ny);
  X=embed(y,de,tau); %embed
end;
[de,nx]=size(X);

%final<ly, fix n...
if isempty(n),
  n=ny;
end;

%and figure out the initial conidtion
if isempty(yi),
  yi=ceil(rand(1,1)*ny);
end;
[di,ni]=size(yi);
if max(di,ni)==1,
  %yi is (hopefully) an integer index.
  if yi~=round(yi),
    disp('WARNING: Rounding yi');
    yi=round(yi);
  end;
  if yi<1 | yi>nx,
    disp('WARNING: yi out of range... guessing');
    yi=ceil(rand(1,1)*ny);
  end;
  yi=X(:,yi);
else,
  %yi is a coordinate point
  if di==1,
    yi=yi';
    [di,ni]=size(yi);
  end;
  if di~=de
    disp('WARNING: yi wrong size... guessing');
    yi=ceil(rand(1,1)*ny);
    yi=X(:,yi);
  end;
end;

%finally... get down to the algorithim.
y=yi;
for i=1:n,

  %find the optimal neighbourhood size np and model order ord.
  [yp,np(i),ord] = localoptimiser(X,yi,de,nx,npmax,ordmax);
  
  %iterate to y(i+1)
  yi = yp + sig*randn(de,1);
  y =[y yi];

end;


  
  
