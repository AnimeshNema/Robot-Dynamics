%theta = [0, -pi/2, 0, 0, 0, pi];  
%alpha = [-pi/2, 0, -pi/2, pi/2, -pi/2, 0];
%d = [290, 0, 0, 302, 0, 72];
%a = [0, 270, 70, 0, 0, 0];

function [T] = dhparam2matrix(theta, d, a, alpha)
T = [cos(theta) (-sin(theta)*cos(alpha)) (sin(theta)*sin(alpha)) (a*cos(theta)); sin(theta) cos(theta)*cos(alpha) (-cos(theta)*sin(alpha)) a*sin(theta); 0 (sin(alpha)) cos(alpha) d; 0 0 0 1];
end




