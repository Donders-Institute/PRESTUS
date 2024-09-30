function [x,y] = getArc(center, radius, angleStart, angleEnd, npoints)
arguments
    center (1,2) double
    radius (1,1) double
    angleStart (1,1) double
    angleEnd (1,1) double
    npoints (1,1) int8 = 200
end
theta = linspace(angleStart, angleEnd, npoints);
x = radius * cos(theta) + center(1); 
y = radius * sin(theta) + center(2); 
end