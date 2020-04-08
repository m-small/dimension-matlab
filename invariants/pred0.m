function y=pred0(y,de,tau,n,yi,sig,np);

%function yp=pred0(y,de,tau,n,yi,noise,np);
%
%An n point local CONSTANT model simulation of y (embedded in de
%dimensions with a lag of tau) starting at yi (either a point or a
%datum number. At each time step, the local linear image is
%obtained by a de point interpolation of the image of the nearest
%np points (default np=de+1);. noise is the standard deviation of 
%gaussian dynamic noise (default 0).
%
%
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk
%

%sort out the parameters
na=nargin;
if na<7,
    np=[];
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
end;
if isempty(sig),
    sig=0;
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
if isempty(np),
    np=de+1;
end;
y=yi;
for i=1:n,
    %find de nearest neighbours of y(i)
    rr=rms((X-yi*ones(1,nx))')';
    rr(end)=[]; %always exclude the last point (it doesn't have an image)
    neigh=[];ndist=[];
    for j=1:np, %get the np closest points (do it this way to save time)
        [minx,mini]=min(rr);
        rr(mini)=inf;
        neigh(j)=mini;
        ndist(j)=minx;
    end;
    %determine weights from ndist
    lam=ones(np,1);
    %lam=1./(ndist+eps);
    lam=lam./sum(lam); %make it sum to one
    lam=lam(:);
    %iterate to find y(i+1)
    yi=X(:,neigh+1)*lam + sig*randn(de,1);
    y=[y yi];
end;

  
  
