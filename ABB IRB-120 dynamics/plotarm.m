function T = plotarm(Q1,Q2,Q3,Q4,Q5,Q6);

syms q1 q2 q3 q4 q5 q6
theta = deg2rad([q1, q2-90, q3, q4, q5, q6+180]);  
alpha = deg2rad([-90, 0, -90, 90, -90, 0]);
d = [290, 0, 0, 302, 0, 72];
a = [0, 270, 70, 0, 0, 0];
T_M = dhparam2matrix(theta, alpha, d, a);

X=[0];
Y=[0];
Z=[0];

for i=1:6
    T_M(:,:,i) = vpa(subs(T_M(:,:,i), [q1, q2, q3, q4, q5, q6],[Q1,Q2,Q3,Q4,Q5,Q6]));
end

T_0_1 = T_M(:,:,1);
T_0_2 = T_M(:,:,1)*T_M(:,:,2);
T_0_3 = T_M(:,:,1)*T_M(:,:,2)*T_M(:,:,3);
T_0_4 = T_M(:,:,1)*T_M(:,:,2)*T_M(:,:,3)*T_M(:,:,4);
T_0_5 = T_M(:,:,1)*T_M(:,:,2)*T_M(:,:,3)*T_M(:,:,4)*T_M(:,:,5);
T_0_6 = T_M(:,:,1)*T_M(:,:,2)*T_M(:,:,3)*T_M(:,:,4)*T_M(:,:,5)*T_M(:,:,6);

X = round([X, T_0_1(1,4) , T_0_2(1,4),  T_0_3(1,4), T_0_4(1,4), T_0_5(1,4), T_0_6(1,4)]);
Y = round([Y, T_0_1(2,4) , T_0_2(2,4),  T_0_3(2,4), T_0_4(2,4), T_0_5(2,4), T_0_6(2,4)]);
Z = round([Z, T_0_1(3,4) , T_0_2(3,4),  T_0_3(3,4), T_0_4(3,4), T_0_5(3,4), T_0_6(3,4)]);
%for i=1:6
 %   P = intermediate(i,0,T_M);
  %  X = round([X, P(1,4)]);
   % Y = round([Y, P(2,4)]);
    %Z = round([Z, P(3,4)]);
%end

figure
scatter3(X,Y,Z);
hold on 
plot3(X,Y,Z);
xlabel('x')
ylabel('y')
zlabel('z')
hold off



