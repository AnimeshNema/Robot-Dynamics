clc; clear all;
syms L1 L2 L3 m1 m2 m3 mL q1 q2 q3 q1_dot q2_dot q3_dot q1_acc q2_acc q3_acc g
%g = 9.8;
old = digits(4);
%% ANS-1 - POSITION KINEMATICS
%Rotation is about x-axis
% Rotation from frame 0 to frame 1
R_0_1 = [1, 0, 0;  0, cos(q1), -sin(q1); 0, sin(q1), cos(q1)];
% Rotation from frame 1 to frame 2
R_1_2 = [1, 0, 0;  0, cos(q2), -sin(q2); 0, sin(q2), cos(q2)];
% Rotation from frame 2 to frame 3
R_2_3 = [1, 0, 0;  0, cos(q3), -sin(q3); 0, sin(q3), cos(q3)];

%displacement along y-axis 
%(projection of point with respect to the previous frame)
P_0_1 = [0 ; L1*cos(q1) ; L1*sin(q1)];
P_1_2 = [0 ; L2*cos(q2) ; L2*sin(q2)];
P_2_3 = [0 ; L3*cos(q3) ; L3*sin(q3)];

%Transformation matrix 
T_0_1 = [R_0_1 , P_0_1 ; 0,0,0,1];
T_1_2 = [R_1_2 , P_1_2 ; 0,0,0,1];
T_2_3 = [R_2_3 , P_2_3 ; 0,0,0,1];
T_0_2 = simplify(T_0_1* T_1_2);

%Homogeneous Transfrmation matrix from base to tip frame
T_0_3 = simplify(T_0_1* T_1_2* T_2_3)

%% ANS-2 VELOCITY KINEMATICS

%position vector for T_0_3 matrix.
P_0_3 = T_0_3(1:3,4);

T_0_2 = T_0_1*T_1_2;

% The rotation component of x-axis (Since, Rotation is about x-axis)
R_0_X1= T_0_1(1:3,1);
R_0_X2= T_0_2(1:3,1);
R_0_X3= T_0_3(1:3,1);

% Differentiating the position with respect to the joint angles
d_P_0_3_1 = diff(P_0_3 , q1);
d_P_0_3_2 = diff(P_0_3 , q2);
d_P_0_3_3 = diff(P_0_3 , q3);

% JACOBIAN FOR FORWARD VELOCITY KINEMATICS 
J = simplify([d_P_0_3_1 , d_P_0_3_2, d_P_0_3_3 ; R_0_X1, R_0_X2, R_0_X3])

%Eliminating linear velocity of x-axis and angular velocity of y & z axis)
Reduced_J = [J(2,:);J(3,:);J(4,:)];

% VELOCITY KINEMATICS EQUATION - (x_dot = J(q) * q_dot)
q_dot = [q1_dot ; q2_dot ; q3_dot]; 
x_dot = simplify((J)*q_dot)

% ANS 3 - NUMERIC SOLUTION
% L1= 60cm, L2= 40cm, L3= 15cm, q1= 30deg, q2= 20deg, q3= -25deg

%Ans 3a
T_0_3_ans3a = vpa(simplify(subs(T_0_3, [q1,q2,q3,L1,L2,L3], [(pi/180)*30,(pi/180)*20,-(pi/180)*25,0.6,0.4,0.15])))
Tip_Position = [T_0_3_ans3a(1:3,4)]
disp('units: m')
rot_angle = vpa(converttoEuler(T_0_3_ans3a))
disp('units: degrees')
%Ans 3b
J_ans3b = vpa(simplify(subs(J, [q1,q2,q3,L1,L2,L3], [(pi/180)*30,(pi/180)*20,-(pi/180)*25,0.6,0.4,0.15])))

%Ans 3c
x_dot_ans3c = vpa(simplify(subs(x_dot, [q1,q2,q3,L1,L2,L3,q1_dot,q2_dot,q3_dot], [(pi/180)*30,(pi/180)*20,-(pi/180)*25,0.6,0.4,0.15,0.785398,0.785398,0.785398])))
disp('units : m/sec')

%% ANS-4 JOINT TORQUE CALCULATION
% Torque = Transpose(J) * Force_at_tip
%Tau = Torque
%mL = mass at tip
%Ftip = Force at tip
Ftip = [0;0;-mL*9.8;0;0;0];

%Ans 4a Symbolic torque matrix
Tau = transpose(J_ans3b)*Ftip

%Ans 4b Numeric torque values
Tau_ans4b = vpa(simplify(subs(Tau, [q1,q2,q3,L1,L2,L3,q1_dot,q2_dot,q3_dot,mL], [(pi/180)*30,(pi/180)*20,-(pi/180)*25,0.6,0.4,0.15,0.785398,0.785398,0.785398,1.5])))
disp('units : N-m ') 

%% ANS-5 Lagrangian equation
syms Q1(t) Q2(t) Q3(t) t  
% Mass 1,2 and 3 are at center of each link and mass mL at end of link3
%Coordinates for mass1
y1 = (L1/2)*cos(Q1);
z1 = (L1/2)*sin(Q1);
%Coordinates for mass2
y2 = (L2/2)*cos(Q1 + Q2) + L1*cos(Q1);
z2 = (L2/2)*sin(Q1 + Q2) + L1*sin(Q1);
%Coordinates for mass3
y3 = L2*cos(Q1 + Q2) + L1*cos(Q1) + (L3/2)*cos(Q1 + Q2 + Q3);
z3 = L2*sin(Q1 + Q2) + L1*sin(Q1) + (L3/2)*sin(Q1 + Q2 + Q3);
%Coordinates for massL
y4 = L2*cos(Q1 + Q2) + L1*cos(Q1) + L3*cos(Q1 + Q2 + Q3);
z4 = L2*sin(Q1 + Q2) + L1*sin(Q1) + L3*sin(Q1 + Q2 + Q3);

%finding velocity (derivate yn,zn with respect to time)
y1_d = diff(y1,t);
z1_d = diff(z1,t);
y2_d = diff(y2,t);
z2_d = diff(z2,t);
y3_d = diff(y3,t);
z3_d = diff(z3,t);
y4_d = diff(y4,t);
z4_d = diff(z4,t);

%square of velocity
v1sq = ((y1_d)^2 + (z1_d)^2);
v2sq = ((y2_d)^2 + (z2_d)^2);
v3sq = ((y3_d)^2 + (z3_d)^2);
v4sq = ((y4_d)^2 + (z4_d)^2);

% Answer 5a
%Total kinetic Energy = (1/2)*m1*(v1*v1)+ (1/2)*m2*(v2*v2)+ (1/2)*m3*(v3*v3)+ (1/2)*mL*(v4*v4)
K = simplify((1/2)*m1*v1sq + (1/2)*m2*v2sq + (1/2)*m3*v3sq + (1/2)*mL*v4sq)

%Total Potential Energy = m1*g*z1 + m2*g*z2 + m3*g*z3 + mL*g*z4
P = simplify(m1*g*z1 + m2*g*z2 + m3*g*z3 + mL*g*z4)

% Answer 5b
% Lagrangian = total kinetic energy - total potential energy (L = K - P)
L = vpa(simplify(K - P));
Lagrangian = vpa(simplify(subs(L,[L1,L2,L3,m1,m2,m3,mL,g], [0.6,0.4,0.15,3,2,0.5,2,9.8])))

% Answer 5c
%Euler Lagrangian Equation: (d/dt)*(dL/d(q_dot)) - dL/dq = Tau
L_ = subs(L,[Q1(t),Q2(t),Q3(t),diff(Q1(t), t),diff(Q2(t), t),diff(Q3(t), t)] ,[q1, q2, q3, q1_dot, q2_dot, q3_dot]);

%dL/dq
dL_dq1 = diff(L_, q1);
dL_dq2 = diff(L_, q2);
dL_dq3 = diff(L_, q3);

%(dL/d(q_dot))
dL_dqdot1 = diff(L_, q1_dot);
dL_dqdot2 = diff(L_, q2_dot);
dL_dqdot3 = diff(L_, q3_dot);

dL_1dq = subs(dL_dq1, [q1, q2, q3, q1_dot, q2_dot, q3_dot], [Q1(t),Q2(t),Q3(t),diff(Q1(t), t),diff(Q2(t), t),diff(Q3(t), t)]);
dL_2dq = subs(dL_dq2, [q1, q2, q3, q1_dot, q2_dot, q3_dot], [Q1(t),Q2(t),Q3(t),diff(Q1(t), t),diff(Q2(t), t),diff(Q3(t), t)]);
dL_3dq = subs(dL_dq3, [q1, q2, q3, q1_dot, q2_dot, q3_dot], [Q1(t),Q2(t),Q3(t),diff(Q1(t), t),diff(Q2(t), t),diff(Q3(t), t)]);

dL_dq1dot = subs(dL_dqdot1, [q1, q2, q3, q1_dot, q2_dot, q3_dot], [Q1(t),Q2(t),Q3(t),diff(Q1(t), t),diff(Q2(t), t),diff(Q3(t), t)]); 
dL_dq2dot = subs(dL_dqdot2, [q1, q2, q3, q1_dot, q2_dot, q3_dot], [Q1(t),Q2(t),Q3(t),diff(Q1(t), t),diff(Q2(t), t),diff(Q3(t), t)]);
dL_dq3dot = subs(dL_dqdot3, [q1, q2, q3, q1_dot, q2_dot, q3_dot], [Q1(t),Q2(t),Q3(t),diff(Q1(t), t),diff(Q2(t), t),diff(Q3(t), t)]);

%(d/dt)*(dL/d(q_dot))
d_dt_dL_dq1dot = diff(dL_dq1dot, t);
d_dt_dL_dq2dot = diff(dL_dq2dot, t);
d_dt_dL_dq3dot = diff(dL_dq3dot, t);

Tau1_ = simplify(d_dt_dL_dq1dot - dL_1dq);
Tau2_ = simplify(d_dt_dL_dq2dot - dL_2dq);
Tau3_ = simplify(d_dt_dL_dq3dot - dL_3dq);
%Tau_  = [Tau1_, Tau2_, Tau3_];

%Torque Values
Tau1 = simplify(subs(Tau1_,[Q1(t),Q2(t),Q3(t),diff(Q1(t), t),diff(Q2(t), t),diff(Q3(t), t), diff(Q1(t), t, t), diff(Q2(t), t, t), diff(Q3(t), t, t)] ,[q1, q2, q3, q1_dot, q2_dot, q3_dot, q1_acc, q2_acc, q3_acc]));
Tau2 = simplify(subs(Tau2_,[Q1(t),Q2(t),Q3(t),diff(Q1(t), t),diff(Q2(t), t),diff(Q3(t), t), diff(Q1(t), t, t), diff(Q2(t), t, t), diff(Q3(t), t, t)] ,[q1, q2, q3, q1_dot, q2_dot, q3_dot, q1_acc, q2_acc, q3_acc]));
Tau3 = simplify(subs(Tau3_,[Q1(t),Q2(t),Q3(t),diff(Q1(t), t),diff(Q2(t), t),diff(Q3(t), t), diff(Q1(t), t, t), diff(Q2(t), t, t), diff(Q3(t), t, t)] ,[q1, q2, q3, q1_dot, q2_dot, q3_dot, q1_acc, q2_acc, q3_acc]));

%Tau Matrix 
Tau = [Tau1; Tau2; Tau3];
Tau = Tau(t);
%Numerical Tau Matrix
Tau_numerical = vpa(simplify(subs(Tau, [L1,L2,L3,m1,m2,m3,mL,g], [0.6,0.4,0.15,3,2,0.5,2,9.8])))

%Intermediate calculations to calculate M, C and G matrix
M11 = simplify((Tau1 - subs(Tau1, q1_acc,0)) / q1_acc);
M12 = simplify((Tau1 - subs(Tau1, q2_acc,0)) / q2_acc);
M13 = simplify((Tau1 - subs(Tau1, q3_acc,0)) / q3_acc);
M21 = simplify((Tau2 - subs(Tau2, q1_acc,0)) / q1_acc);
M22 = simplify((Tau2 - subs(Tau2, q2_acc,0)) / q2_acc);
M23 = simplify((Tau2 - subs(Tau2, q3_acc,0)) / q3_acc);
M31 = simplify((Tau3 - subs(Tau3, q1_acc,0)) / q1_acc);
M32 = simplify((Tau3 - subs(Tau3, q2_acc,0)) / q2_acc);
M33 = simplify((Tau3 - subs(Tau3, q3_acc,0)) / q3_acc);


%C11 = simplify((Tau1 - subs(Tau1, q1_dot,0)) / q1_dot);
%C12 = simplify((Tau1 - subs(Tau1, q2_dot,0)) / q2_dot);
%C13 = simplify((Tau1 - subs(Tau1, q3_dot,0)) / q3_dot);
%C21 = simplify((Tau2 - subs(Tau2, q1_dot,0)) / q1_dot);
%C22 = simplify((Tau2 - subs(Tau2, q2_dot,0)) / q2_dot);
%C23 = simplify((Tau2 - subs(Tau2, q3_dot,0)) / q3_dot);
%C31 = simplify((Tau3 - subs(Tau3, q1_dot,0)) / q1_dot);
%C32 = simplify((Tau3 - subs(Tau3, q2_dot,0)) / q2_dot);
%C33 = simplify((Tau3 - subs(Tau3, q3_dot,0)) / q3_dot);

G1 = simplify(subs(Tau1, [q1_acc, q2_acc, q3_acc, q1_dot, q2_dot, q3_dot],[0,0,0,0,0,0]));
G2 = simplify(subs(Tau2, [q1_acc, q2_acc, q3_acc, q1_dot, q2_dot, q3_dot],[0,0,0,0,0,0]));
G3 = simplify(subs(Tau3, [q1_acc, q2_acc, q3_acc, q1_dot, q2_dot, q3_dot],[0,0,0,0,0,0]));

%Inertia Matrix
M = simplify([M11,M12,M13; M21,M22,M23; M31,M32,M33])
M = M(t);
%Numerical Inertia Matrix
M_numerical = vpa(simplify(subs(M,[L1,L2,L3,m1,m2,m3,mL,g],[0.6,0.4,0.15,3,2,0.5,2,9.8])))

%Gravity Matrix
G = simplify([G1; G2; G3])
G = G(t);
%Numerical Gravity Matrix
G_numerical = vpa(simplify(subs(G,[L1,L2,L3,m1,m2,m3,mL,g],[0.6,0.4,0.15,3,2,0.5,2,9.8])))

% Coriollis force matrix
q_acc = [q1_acc; q2_acc; q3_acc];
Cq = [Tau - M*q_acc - G];
C11 = simplify(Cq(1)/q1_dot);
C12 = simplify(Cq(1)/q2_dot);
C13 = simplify(Cq(1)/q3_dot);
C21 = simplify(Cq(2)/q1_dot);
C22 = simplify(Cq(2)/q2_dot);
C23 = simplify(Cq(2)/q3_dot);
C31 = simplify(Cq(3)/q1_dot);
C32 = simplify(Cq(3)/q2_dot);
C33 = simplify(Cq(3)/q3_dot);

C = simplify([C11,C12,C13; C21,C22,C23; C31,C32,C33])
%Numerical Coriollis force Matrix
C_numerical = vpa(simplify(subs(C,[L1,L2,L3,m1,m2,m3,mL,g],[0.6,0.4,0.15,3,2,0.5,2,9.8])))

% Matrix-Vector Relationship: M(q)q_acc + C(q,q_dot)q_dot + G(q) = Tau
Matrix_Vector_Relationship = M*q_acc + C*q_dot + G == Tau

%Numerical Relation
Numerical_Matrix_Vector_Relationship = M_numerical*q_acc + C_numerical*q_dot + G_numerical == Tau_numerical

%% ANS 6 Newton Euler Dynamics
%initial state
w0 = 0;
acc0_e = 0;
acc0_c = 0;
acc0 = 0;

%given link lenghts
r_0_1 = [0;0.6;0];
r_1_2 = [0;0.4;0];
r_2_3 = [0;0.15;0];
rc_0_1 =[0;0.3;0];
rc_1_2 =[0;.2;0];
rc_2_3 =[0;0.075;0];
rc_4_3 =[0;-0.075;0];
rc_1_0 =[0;-0.3;0];
rc_2_1 =[0;-.2;0];

%forward recursion

%Angular velocity 
w1 = transpose(R_0_1)*[0;0;0] + [1;0;0]*q1_dot;
w2 = transpose(R_1_2)*w1 + [1;0;0]*q2_dot;
w3 = transpose(R_2_3)*w2 + [1:0;0]*q3_dot;

qb1 = (q1_dot*[1;0;0]);
qb2 = (q2_dot*[1;0;0]);
qb3 = (q3_dot*[1;0;0]);
%Angular accelaration
acc1 = transpose(R_0_1)*[0;0;0] + [1;0;0]*q1_acc + cross(w1,qb1);
acc2 = transpose(R_1_2)*acc1 + [1;0;0]*q2_acc + cross(w2,qb2);
acc3 = transpose(R_2_3)*acc2 + [1;0;0]*q3_acc + cross(w3,qb3);

%linear accelaration
acc1_e = transpose(R_0_1)*[0;0;0] + cross(acc1,r_0_1) + cross(w1,(cross(w1,r_0_1)));
acc2_e = transpose(R_1_2)*acc1_e + cross(acc2,r_1_2) + cross(w2,(cross(w2,r_1_2)));
acc3_e = transpose(R_2_3)*acc2_e + cross(acc3,r_2_3) + cross(w3,(cross(w3,r_2_3)));

acc1_c = transpose(R_0_1)*[0;0;0] + cross(acc1,rc_0_1) + cross(w1,cross(w1,rc_0_1));
acc2_c = transpose(R_1_2)*acc1_e + cross(acc2,rc_1_2) + cross(w2,cross(w2,rc_1_2));
acc3_c = transpose(R_2_3)*acc2_e + cross(acc3,rc_2_3) + cross(w3,cross(w3,rc_2_3));

%backward recursion

%force calculation
f4 = [0;0;2*9.8];
R_4_3 = eye(3);
R_0_3 = R_0_1*R_1_2*R_2_3;

%for converting f4 from global frame to frame 3
R_3_0 = transpose(R_0_3);

f3 = R_3_0*f4 + mL*acc3_c - mL*g;
f2 = transpose(R_2_3)*f3 + mL*acc2_c - m3*9.8;
f1 = transpose(R_1_2)*f2 + mL*acc1_c - m2*9.8;


%Torque Calculation
T4 = 0;
T3 = R_4_3*[0;0;0] - cross(f3,rc_2_3) + cross((R_4_3*f4), rc_4_3) + acc3;
T2 = transpose(R_2_3)*T3 - cross(f2,rc_1_2) + cross((transpose(R_2_3)*f3), rc_2_1) + acc2;
T1 = transpose(R_1_2)*T2 - cross(f1,rc_0_1) + cross((transpose(R_1_2)*f2), rc_1_0) + acc1;

T_newton = [T1(1,:);T2(1,:);T3(1,:)];

%Intermediate calculations to calculate M, C and G matrix
M_11 = simplify((T1(1,:) - subs(T1(1,:), q1_acc,0)) / q1_acc);
M_12 = simplify((T1(1,:) - subs(T1(1,:), q2_acc,0)) / q2_acc);
M_13 = simplify((T1(1,:) - subs(T1(1,:), q3_acc,0)) / q3_acc);
M_21 = simplify((T2(1,:) - subs(T2(1,:), q1_acc,0)) / q1_acc);
M_22 = simplify((T2(1,:) - subs(T2(1,:), q2_acc,0)) / q2_acc);
M_23 = simplify((T2(1,:) - subs(T2(1,:), q3_acc,0)) / q3_acc);
M_31 = simplify((T3(1,:) - subs(T3(1,:), q1_acc,0)) / q1_acc);
M_32 = simplify((T3(1,:) - subs(T3(1,:), q2_acc,0)) / q2_acc);
M_33 = simplify((T3(1,:) - subs(T3(1,:), q3_acc,0)) / q3_acc);

G1_ = simplify(subs(T1(1,:), [q1_acc, q2_acc, q3_acc, q1_dot, q2_dot, q3_dot],[0,0,0,0,0,0]));
G2_ = simplify(subs(T2(1,:), [q1_acc, q2_acc, q3_acc, q1_dot, q2_dot, q3_dot],[0,0,0,0,0,0]));
G3_ = simplify(subs(T3(1,:), [q1_acc, q2_acc, q3_acc, q1_dot, q2_dot, q3_dot],[0,0,0,0,0,0]));

%Inertia Matrix
M_ = simplify([M_11,M_12,M_13; M_21,M_22,M_23; M_31,M_32,M_33])
%M_ = M(t);
%Numerical Inertia Matrix
M_numerical_ = vpa(simplify(subs(M_,[L1,L2,L3,m1,m2,m3,mL,g],[0.6,0.4,0.15,3,2,0.5,2,9.8])))

%Gravity Matrix
G_ = simplify([G1_; G2_; G3_])
%G = G(t);
%Numerical Gravity Matrix
G_numerical_ = vpa(simplify(subs(G_,[L1,L2,L3,m1,m2,m3,mL,g],[0.6,0.4,0.15,3,2,0.5,2,9.8])))

% Coriollis force matrix
q_acc = [q1_acc; q2_acc; q3_acc];
Cq_ = [T_newton - M_*q_acc - G_];
C11_ = simplify(Cq_(1,:)/q1_dot);
C12_ = simplify(Cq_(1,:)/q2_dot);
C13_ = simplify(Cq_(1,:)/q3_dot);
C21_ = simplify(Cq_(2,:)/q1_dot);
C22_ = simplify(Cq_(2,:)/q2_dot);
C23_ = simplify(Cq_(2,:)/q3_dot);
C31_ = simplify(Cq_(3,:)/q1_dot);
C32_ = simplify(Cq_(3,:)/q2_dot);
C33_ = simplify(Cq_(3,:)/q3_dot);

C_ = simplify([C11_,C12_,C13_; C21_,C22_,C23_; C31_,C32_,C33_])
%Numerical Coriollis force Matrix
C_numerical_ = vpa(simplify(subs(C_,[L1,L2,L3,m1,m2,m3,mL,g],[0.6,0.4,0.15,3,2,0.5,2,9.8])))

% Matrix-Vector Relationship: M(q)q_acc + C(q,q_dot)q_dot + G(q) = Tau
Matrix_Vector_Relationship = M_*q_acc + C_*q_dot + G_ == T_newton

%Numerical Relation
Numerical_Matrix_Vector_Relationship = M_numerical_*q_acc + C_numerical_*q_dot + G_numerical_ == T_newton

%ans 6 c
% difference in both methods 
M_diff = simplify(M_numerical_ - M_numerical)
C_diff = simplify(C_numerical_ - C_numerical)
G_diff = simplify(G_numerical_ - G_numerical)

% function to convert angles to Euler
function rot_angle = converttoEuler(Tfinal)
% Convert rotation matrix to Euler angles
singularity_check = sqrt(Tfinal(1,1)^2 + Tfinal(2,1)^2);
if singularity_check < 10^-6
    PHI = atan2(Tfinal(2,3), Tfinal(2,2)); % corresponds to RX in Robot Studio
    THETA = atan2(-Tfinal(3,1), sqrt(Tfinal(1,1)^2+Tfinal(2,1)^2)); % corresponds to RY in Robot Studio
    PSI = 0; % corresponds to RZ in Robot Studio
else
    THETA = atan2(-Tfinal(3,1),sqrt(Tfinal(1,1)^2+Tfinal(2,1)^2)); % corresponds to RY in Robot Studio
    PHI = atan2(Tfinal(3,2),Tfinal(3,3)); % corresponds to RX in Robot Studio
    PSI = atan2(Tfinal(2,1),Tfinal(1,1)); % corresponds to RZ in Robot Studio
end

rot_angle = [PHI,THETA,PSI]*180/pi;
end



