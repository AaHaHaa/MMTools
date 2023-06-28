clearvars; close all;

s = 4;
m=rand(s,s,1)+1i*rand(s,s,1); m = (m+pagectranspose(m))/2;

m = gpuArray(m);

tic;
expD = pageexpm_skew_Hermitian(m);
a=expD(:,:,1)-expm(m(:,:,1));max(abs(a(:)))
toc

addpath('expm/');
tic;
b=myexpm_(m,[],[],1,0,1);
c=b(:,:,1)-expm(m(:,:,1));max(abs(c(:)))
toc