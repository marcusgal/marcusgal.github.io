function [ ] = plot_circ_in_square( coords)
%PLOT_CIRC_IN_SQUARE Plots circles in a square (2 dimensions only)
%  coords - The coordinates for the circle centres, in the range [0,1]
%    Coordinates will be adjusted.

numcircs = floor(size(coords,2) / 2);

% evaluate points
mindist = circ_in_square( numcircs, 2, [1;1], coords );

radius = (abs(mindist) / 2);

xcoords = (radius + coords(1:numcircs))                / (1 + 2*radius);
ycoords = (radius + coords((numcircs + 1):2*numcircs)) / (1 + 2*radius);

radius = radius /(1+2*radius);

% Centre points
scatter(xcoords, ycoords, '.');
hold on;
xlim([0,1]);
ylim([0,1]);

% Draw cicles
numpoints = 64; % how many points in a circle
for i = 1:numcircs
    % from some code on drawing circles
    theta=linspace(0,2*pi,numpoints);
    rho=ones(1,numpoints)*radius;
    [circX,circY] = pol2cart(theta,rho);
    circX=circX+xcoords(i);
    circY=circY+ycoords(i);
    plot(circX,circY,'b-');
end
hold off;