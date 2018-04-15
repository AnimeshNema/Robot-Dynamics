clc,clear;

syms q1 q2 q3 q4 q5 q6
theta = deg2rad([q1, q2-90, q3, q4, q5, q6+180]);  
alpha = deg2rad([-90, 0, -90, 90, -90, 0]);
d = [290, 0, 0, 302, 0, 72];
a = [0, 270, 70, 0, 0, 0];

T_M = dhparam2matrix(theta, alpha, d, a)

T_total = T_M(:,:,1)*T_M(:,:,2)*T_M(:,:,3)*T_M(:,:,4)*T_M(:,:,5)*T_M(:,:,6)
T_simplified= simplify(T_total)

T_0_6_Q5 = subs(T_simplified, [q1,q2,q3,q4,q5,q6], [0,0,0,0,0,0]);
T_0_6_Q5 = round(T_0_6_Q5)

T_0_6_Q6 = subs(T_simplified, [q1,q2,q3,q4,q5,q6], [-45,30,-30,-30,-45,180]);
T_0_6_Q6 = round(T_0_6_Q6)

%T1 = plotarm(-50,40,-30,-30,-50,180)
T_0_6_Q8_1 = subs(T_simplified, [q1,q2,q3,q4,q5,q6], [-50,40,-30,-30,-50,180]);
T_0_6_Q8_1 = round(T_0_6_Q8_1)
%T2 = plotarm(-11,-4,26,-9,-27,-16)
T_0_6_Q8_2 = subs(T_simplified, [q1,q2,q3,q4,q5,q6], [-11,-4,26,-9,-27,-16]);
T_0_6_Q8_2 = round(T_0_6_Q8_2)
%T3 = plotarm(-5,50,-36,-24,-36,0)
T_0_6_Q8_3= subs(T_simplified, [q1,q2,q3,q4,q5,q6], [-5,50,-36,-24,-36,0]);
T_0_6_Q8_3 = round(T_0_6_Q8_3)
%T4 = plotarm(-45,35,-35,-35,-40,0)
T_0_6_Q8_4 = subs(T_simplified, [q1,q2,q3,q4,q5,q6], [-45,35,-35,-35,-40,0]);
T_0_6_Q8_4 = round(T_0_6_Q8_4)
%T5 = plotarm(-10,-5,-20,-20,-10,0);
T_0_6_Q8_5 = subs(T_simplified, [q1,q2,q3,q4,q5,q6], [-10,-5,-20,-20,-10,0]);
T_0_6_Q8_5 = round(T_0_6_Q8_5)




