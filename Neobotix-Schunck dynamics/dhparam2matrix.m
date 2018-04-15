function M = dhparam2matrix(theta, disp_z, alpha, disp_x)
old = digits(3);
M = cell(5,1);
for i = 1:5
    M{i} = [cos(theta(i)) (-sin(theta(i))*cos(alpha(i))) (sin(theta(i))*sin(alpha(i))) (disp_x(i)*cos(theta(i)));
       sin(theta(i)) cos(theta(i))*cos(alpha(i)) (-cos(theta(i))*sin(alpha(i))) disp_x(i)*sin(theta(i));
       0 (sin(alpha(i))) cos(alpha(i)) disp_z(i);
       0 0 0 1];
end
end