function [center, Area] = centroidHull(x,y)% x & y being vectors with all the coordinates
k = convhull(x,y);
[xc, yc] = centroid(polyshape(x(k),y(k)));
Area = area(polyshape(x(k),y(k))); % plus area.GM, Oct 2019
center = [xc yc];
end