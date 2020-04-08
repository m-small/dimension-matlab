function [dl,off,lnx,lne,aopt]=window(y,maxd,approx);

%function [dl,off,lnx,lne,aopt]=window(y,maxd,approx);
%
%compute DL(d) for 1:maxd to determine best embedding window
%approx is optional and specifies the order of the model approximation
% (0 linear ala Kantz, 1 constant, 2 constant optimal selection, 
% 3 constant optimised with a simplex search)
%
% DL=d+N/2*(1+ln(2*pi))+(N-d)/2*ln(E(e^2))+d/2*ln(E(x^2))
%
% dl is the full expression
% off is the first two terms
% lne is the third 
% lnx is the fourth
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk
%

linear=0;

if nargin<2,
    maxd=100;
end;
if nargin<3,
    approx=~linear;
end;
y=y(:);

n=length(y);
d=1:maxd;
lnx=d/2.*log(cumsum( (y(d)-mean(y)).^2)'./d );
off=n/2*(1+log(2*pi))+d;
aopt=[];

h=waitbar(0,'window working...');
for i=1:maxd,
    if approx==0,
        %Kantz's local linear thingy (too many params)
        err(i)=prederr(y,i,1);
    else,
        x=embed(y,i,1);
        %find the best "projection"
        if approx>2,
            %do full nonlinear simplex search
            aopt=fminsearch(inline('rms(x(1,2:end)-x(1,nearest(x(:,1:(end-1)),0,avect)))','avect','x'),ones(1,i),[],x);
        elseif approx==2,
            %add iff it improves
            aopt=[aopt 1];
            if i>1,
                %test new extended embedding ... 
                xi=nearest(x(:,1:(end-1)),0,aopt);
                err(i)=rms(x(1,2:end)-x(1,xi+1));
                %is it any better (in terms of RMS(err))?
                if err(i)>err(i-1),
                    aopt(end)=0;
                end;
           %     olddl=newdl;
           %     newdl=log(err(i).^2).*(n-d(i))./2+off(i)+lnx(i);
           %     if olddl<newdl
           %         aopt(end)=0;
           %     end;
            else
           %     xi=nearest(x(:,1:(end-1)),0,aopt);
           %     err(i)=rms(x(1,2:end)-x(1,xi+1));
           %     newdl=log(err(i).^2).*(n-d(i))./2+off(i)+lnx(i);
            end;
        else,
            %include all (constant Euclidean norm)
            aopt=ones(1,i);  
        end;
        %simple, nearest neighbour analogue        
        xi=nearest(x(:,1:(end-1)),0,aopt);
        err(i)=rms(x(1,2:end)-x(1,xi+1));
    end;
    waitbar(i/maxd,h,['window working... ',int2str(i),'/',int2str(maxd)]);
end;

if approx==0,
    err=err.*std(y);
end;
lne=(n-d)/2.*log(err.^2);
dl=off+lnx+lne;

close(h);

plot(dl,'r');
hold on;
plot(off,'g-.');
plot(lnx,'k-');
plot(lne,'b-');
plot(dl,'ro');
[m,i]=min(dl);
plot(i,m,'rp');
ax=axis;
[m,i]=min(lne);
plot([i i]',ax(3:4),'b:');
hold off;
