%% Simulation of uniformly distributed points on a spherical shell 
% e.g. to determine rate of random colocalisation of nuclear bodies in Arabidopsis nucleoplasm

Ro = 1;   %outer radius of nucleus ~5 microns
Ri = 0.9999;  %
d = 0.04;  % overlap resolution is ~200 nm

N = 10000;  % number of sample points

%% Exact solution for ball
L=linspace(0,2*Ro,N);
%for the complete ball R_i=0, the pdf is exactly:
pdf_ball = 3*(L/Ro).^2-9/4*(L/Ro).^3+3/16*(L/Ro).^5;
%for which the mean distance is 
Lmean_ball = 36/35*Ro;
plot(L,pdf_ball)

%% Exact solution for sphere 
pdf_sphere=3/4*Ro*sqrt((L/Ro).^2-1/4*(L/Ro).^4);
hold on
plot(L,pdf_sphere)
plot(L,pdf_sphere.^(1/5))
%% Numerical solution for spherical shell

V = 4/3*pi*(Ro^3-Ri^3);
A = pi*(Ro^2-Ri^2);


Cv = N/V
Ca = N/A

%sample N uniform points on a spherical shell using polar coordinates

u = 2*rand(N,1)-1;     %axial coordinate
phi = 2*pi*rand(N,1);  %azimuthal coordinate
r = (Ri^3+(Ro^3-Ri^3).*rand(N,1)).^(1/3);   %If PDF is uniform over volume, CDF has to ~shell volume, so shell radius ~CDF^(1/3)
x = r.*cos(phi).*(1-u.^2).^0.5;
y = r.*sin(phi).*(1-u.^2).^0.5;
z = r.*u;

%calculate pairwise distances along z-projection
%delta = sqrt((x-x(1)).^2+(y-y(1)).^2+(z-z(1)).^2);
%delta = pdist([x;y]);
delta = pdist([x,y,z]);

%calculate integral up to overlapping spot size d
randomcoloc = sum(sum(delta<d))/N^2

plot(delta(1:300))

%generate probability density function of pairwise distance
h = histogram(delta,'Normalization','probability');
probdelta = h.Values;
bin = h.BinWidth;

%mean distance
Lmean = sum(bin*(1:length(probdelta)).*probdelta)/sum(probdelta)

%plot point clouds
f2 = figure;
plot3(x,y,z);pbaspect([1,1,1]);

f3 = figure;
scatter(x,y);pbaspect([1,1,1]);


