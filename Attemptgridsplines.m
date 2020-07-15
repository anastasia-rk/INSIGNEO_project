%input variables 
grid_limits = [0 0 100 100]
% order of the splines
order = 4
% number of splines along x and y axis 
no.splinesx = 5;
no.splinesy = 5;
%start and end points of the grid of splines
ax = 10; bx = 70;
ay = 10; by = 70;


dx = (grid_limits(3)-grid_limits(1))/100;
dy = (grid_limits(4)-grid_limits(2))/100;
xval = grid_limits(1):dx:grid_limits(3);
yval = grid_limits(2):dy:grid_limits(4);
length_splinex = (bx-ax)/no.splinesx;
length_spliney= (by-ay)/no.splinesy;
j = ax:length_splinex:bx;
i = ay:length_spliney:by;
ll = no.splinesx*no.splinesy;
theta = 10*1:ll;
    for k = 1:length(i)-1
        for n = 1:length(j)-1
surfb(4,j(1,n),j(1,n+1),i(1,k),i(1,k+1),xval,yval,n*k)
hold on
        end 
    end 
   

    
   


 




