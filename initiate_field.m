function[Theta] = initiate_field(body,wound,limits,knots)
dx = (limits(3) - limits(1))/100;
dy = (limits(4) - limits(2))/100;
coordinate_x = [limits(1):dx:limits(3)];
coordinate_y = [limits(2):dy:limits(4)];
N = length(coordinate_x);
M = length(coordinate_y);
 
[X_grid, Y_grid] = meshgrid(coordinate_x,coordinate_y);

% z = ax + by + c;
c = (body*limits(3) - wound*limits(1))/(limits(1)+limits(3));
a = (wound - c)/limits(3);
b = 0;
% The surface approximation is Z = A*theta
% generate vector z
for i=1:N
    for j=1:M
        Z((i-1)*M + j,1) = a*coordinate_x(i) + b*coordinate_y(j) + c;
        Z_plot(j,i) =  Z((i-1)*M + j,1);
    end
end

% figure; 
% plot3(X_grid,Y_grid,Z_plot);

% generate matrix A
ll = size(knots,2) - 1;
for i=1:N
    for j=1:M
        k = 0;
        for index = 1:2:ll 
            support_x = knots(1,index:index+1);
            support_y = knots(2,index:index+1);
            k = k + 1;
            A((i-1)*M + j,k) = tensorproductbspline(4,support_x(1,1),support_x(1,2),support_y(1,1),support_y(1,2),coordinate_x(i),coordinate_y(j));           
        end
    end
end
% Least square estimate
 Theta = A\Z;
end