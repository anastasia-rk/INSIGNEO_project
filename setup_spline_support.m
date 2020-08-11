function[knots] = setup_spline_support(grid_limits,nx,ny,order)
% function creates a grid of cardinally spaced nodes for B-spline support
% Input:
% grid_limits (vector 4x1) - the size of the area of interest 
% (for image    - % pixels in % formal [x_start y_start x_end y_end])
% nx (scalar)   - number of splines along x-axis
% ny (scalar)   - number of splines along y-axis
% Output:
% knots (2 x 2*nx*ny) - ordered pairs of coordinates for each 2d spline 
%% Number of intervals between nodes in each direction - specifically for cubic b-splines7
order = order - 1;
knot_nx = nx + order;
knot_ny = ny + order;
% length between nodes in each direction
dx = (grid_limits(3)- grid_limits(1))/(knot_nx);
dy = (grid_limits(4)- grid_limits(2))/(knot_ny);
% positions of nodes
for i=1:knot_nx+1
    x_knot(i) = grid_limits(1) + dx*(i-1);
end
for i=1:knot_ny+1
    y_knot(i) = grid_limits(2) + dy*(i-1);
end
k = 1;
% for each scubic spline support knots are shifted by one quarter of the length
for i=1:knot_nx-order
    for j=1:knot_ny-order
        knots(1,k)      = x_knot(i); 	% left-most point, x coordinate
        knots(1,k+1) 	= x_knot(i+4); 	% right-most point, x coordinate
        knots(2,k)      = y_knot(j); 	% left-most point, y coordinate
        knots(2,k+1) 	= y_knot(j+4);  % right-most point, x coordinate
        k = k + 2;
    end
end
