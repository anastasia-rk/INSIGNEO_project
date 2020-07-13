function[knots,sigma] = setup_gaussian_support(grid_limits,nx,ny)
% function creates a grid of cardinally spaced nodes for B-spline support
% Input:
% grid_limits (vector 4x1) - the size of the area of interest 
% (for image    - % pixels in % formal [x_start y_start x_end y_end])
% nx (scalar)   - number of splines along x-axis
% ny (scalar)   - number of splines along y-axis
% Output:
% knots (2 x nx*ny) - ordered pairs of coordinates for each 2d Gaussian function 
%%
dx = (grid_limits(3)- grid_limits(1))/(nx+1);
dy = (grid_limits(4)- grid_limits(2))/(ny+1);
% Set the width of the basis function to the larger interval
if dx > dy
    sigma = dx/2;
else 
    sigma = dy/2;
end
% for radial basis functions knots are located at centers of bf
for i=1:nx
    x_knot(i) = grid_limits(1) + dx*(i);
end
for i=1:ny
    y_knot(i) = grid_limits(2) + dy*(i);
end
% combine centers into knot array
k = 1;
for i=1:nx
    for j=1:ny
        knots(1,k) = x_knot(i);
        knots(2,k) = y_knot(j);
        k = k + 1;
    end
end



