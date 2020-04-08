%you've unpacked it and compiled all the mex files.
%now test it with the following

xi=randn(2,1);
for i=1:1000,
  xi=ikeda(xi);
end;
x=[];
for i=1:3000,
  xi=ikeda(xi)+randn(2,1)*0.01;
  x=[x xi];
end;

z=x(1,:)+randn(1,3000)*0.02;

plot(z(1:(end-1)),z(2:end),'.');

%compute complexity
cmp=complexity(z,2:10);

%compute dimensions
figure;

[m,d,k,s]=gka(z,2:20,1,250,3);
