


% Define 3D data
[r, g, b] = meshgrid(0:1:255);
gr = floor(0.299*r + 0.587*g + 0.114*b);

% define the slice planes
%[xi, yi] = meshgrid(0:1:255);

xi = [];
yi = [];
zi = 20;

% slice(r,g,b,gr,xi,yi,zi);
% colormap(gray(256));
% shading interp

volumeViewer
