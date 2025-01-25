function [fig,ax] = radialPcolor(r,A2)
%RADIALPCOLOR Summary of this function goes here
%   Detailed explanation goes here

r_max = max(r); % m

Nr = length(r);
Nx = Nr*2;
x = linspace(-r_max,r_max,Nx); % m
[xx,yy] = meshgrid(x,x);

rr = sqrt(xx.^2+yy.^2);

A2 = reshape(interp1(r.',A2.',rr(:),'linear',0),Nx,Nx);

fig = figure;
pcolor(x*1e6,x*1e6,A2);
shading interp; colormap(jet);
xlabel('x (\mum)');
ax = gca;

end