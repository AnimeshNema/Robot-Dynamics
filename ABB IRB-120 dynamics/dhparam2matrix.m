
function T_M = dhparam2matrix(theta, alpha, d, a)
old = digits(2);
for i= 1:6
   T_M(:,:,i)= vpa([cos(theta(i)) (-sin(theta(i))*cos(alpha(i))) (sin(theta(i))*sin(alpha(i))) (a(i)*cos(theta(i)));
       sin(theta(i)) cos(theta(i))*cos(alpha(i)) (-cos(theta(i))*sin(alpha(i))) a(i)*sin(theta(i));
       0 (sin(alpha(i))) cos(alpha(i)) d(i);
       0 0 0 1]);
   
end
end

