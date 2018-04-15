clc; clear all;
syms q1 q2 q3 b c d e f g alp1 beta1 alp2 beta2 alp3 beta3 a wL wr wo ws thet rf ro rs x_o y_o thet_o x_qR y_qR z_qR pi

%% D-H PARAMETERS
theta = [0, q1-(pi/2), q2, q3, -(pi/2)];
disp_z =[b+d, 0, 0, 0, 0];
alpha = [-(pi/2), 0, 0, 0, -(pi/2)];
disp_x= [-c, e, f, g, 0];

%% ANSWER 2
M = dhparam2matrix(theta, disp_z, alpha, disp_x);
for n = 1:5
    %disp(sprintf('T_%d_%d', n-1 , n))
    vpa(M{n});
end

T_0_5 = M{1}*M{2}*M{3}*M{4}*M{5};

%% THE TRANSFORMATION MATRIX FROM BASE TO TIP FRAME
T_0_5 = simplify(T_0_5)

%% ANSWER 3

% The position vector for T_0_5
X_0_5 = T_0_5(1:3,4);

T_0_1 = M{1};
T_0_2 = M{1}*M{2};
T_0_3 = M{1}*M{2}*M{3};

% The rotation matrix 
R_0_Z1= T_0_1(1:3,3);
R_0_Z2= T_0_2(1:3,3);
R_0_Z3= T_0_3(1:3,3);

% Differentiating the position with respect to the joint angles
d_X_0_5_1 = diff(X_0_5 , q1);
d_X_0_5_2 = diff(X_0_5 , q2);
d_X_0_5_3 = diff(X_0_5 , q3);

%% JACOBIAN FOR FORWARD VELOCITY KINEMATICS 
%% (x_dot = J(q) * q_dot) 
J = [d_X_0_5_1 , d_X_0_5_2, d_X_0_5_3 ; R_0_Z1, R_0_Z2, R_0_Z3]

%Reduced Jacobain(Linear y velocity as well as angular x and z are zero)
Reduced_J_matrix= [J(1,:);J(3,:);J(5,:)]

%% ANSWER 4- TRANSPOSE OF JACOBIAN FOR FORCE-TORQUE RELATIONSHIP
%% Tau = Transpose(J) * Ftip
J_T = transpose(J)

%% ANSWER 5-A 
% Transformation matrix from base to tip frame for the given Configuration 
T_0_5_ans5a = vpa(simplify(subs(T_0_5, [q1,q2,q3,b,c,d,e,f,g], [(pi/180)*20,(pi/180)*90,(pi/180)*30,.424,.300,.380,.328,.323,.0824])))

%% ANSWER 5-B
% Forward Velocity Kinematics for the given configuration
q_dot = [(pi/180)*30; (pi/180)*30; (pi/180)*30];
J_ans5b = vpa(simplify(subs(J, [q1,q2,q3,b,c,d,e,f,g], [(pi/180)*20,(pi/180)*90,(pi/180)*30,.424,.300,.380,.328,.323,.0824])));
X_dot = vpa(simplify(J_ans5b*q_dot))
disp('units : m/sec and rad/sec') 

%% ANSWER 5-C
% Calculating Torque through Torque-force relationship.
F_tip = [30; 0; 0; 0; 0; 0];
%J_T_ans5b = vpa(simplify(subs(J_T, [q1,q2,q3,b,c,d,e,f,g], [20,90,30,.424,.300,.380,.328,.323,.0824])));
Tau = vpa(simplify(transpose(J_ans5b)*F_tip))
disp('units : N-m ') 

%% ANSWER 6 
%% INVERSE VELOCITY KINEMATICS
x_dot = [ 0; 0; -.1; 0; 0; 0];
Q_dot = vpa(simplify(pinv(J_ans5b)*x_dot))
disp('units :rad/sec ') 

%% ANSWER 7 
% The rotation constraint matrix for left fixed wheel
FW1_RCM = [sin(alp1 + beta1), -cos(alp1 + beta1) , -(a/2)*cos(beta1)];
% The sliding constraint matrix for left fixed wheel 
FW1_SCM = [cos(alp1 + beta1), sin(alp1 + beta1) , -(a/2)*sin(beta1)];
% The rotation constraint matrix for right fixed wheel
FW2_RCM = [sin(alp2 + beta2), -cos(alp2 + beta2) , -(a/2)*cos(beta2)];
% The sliding constraint matrix for right fixed wheel 
FW2_SCM = [cos(alp2 + beta2), sin(alp2 + beta2) , -(a/2)*sin(beta2)];
% The rotation constraint matrix for omni wheel
OW3_RCM = [sin(alp3 + beta3), -cos(alp3 + beta3) , -(a/2)*cos(beta3)];
% The sliding constraint matrix for omni wheel 
OW3_SCM = [cos(alp3 + beta3), sin(alp3 + beta3) , (a/2)*sin(beta3)];

%% Rotation matrix from robot to world frame 
R_rob_0 = [cos(thet), sin(thet), 0; -sin(thet), cos(thet), 0; 0, 0, 1];

R_0_rob = inv(R_rob_0);

c_dot= [x_o ; y_o; thet_o];

%% SIX EQUATIONS FOR THE KINEMATIC CONSTRAINTS
% The rolling kinematic constraint equation for left fixed wheel
Eq1 = (FW1_RCM) * ((R_rob_0) * (c_dot)) - (rf*wL) == 0
% The sliding kinematic constraint equation for left fixed wheel
Eq2 = (FW1_SCM) * ((R_rob_0) * (c_dot)) == 0
% The rolling kinematic constraint equation for right fixed wheel
Eq3 = (FW2_RCM) * ((R_rob_0) * (c_dot)) - (rf*wr) == 0
% The sliding kinematic constraint equation for right fixed wheel
Eq4 = (FW2_SCM) * ((R_rob_0) * (c_dot)) == 0
% The rolling kinematic constraint equation for omni wheel
Eq5 = (OW3_RCM) * ((R_rob_0) * (c_dot)) - (ro*wo) == 0
% The sliding kinematic constraint equation for omni wheel
Eq6 = (OW3_SCM) * ((R_rob_0) * (c_dot)) -(rs*ws) == 0

%% ANSWER 8

% Combined Rolling Constraint Matrix
J_M = [ FW1_RCM ; FW2_RCM ; OW3_RCM]
% Combined Sliding Constraint Matrix
C_M = [ FW1_SCM ; FW2_SCM ; OW3_SCM]

r_J2_M = [ rf*wL; rf*wr ; ro*wo];
r_C2_M = [0;0;rs*ws];

JC = [J_M ; C_M]

% THE COMBINED EQUATION FOR THE KINEMATIC CONSTRAINTS
Equation = (JC) * ((R_rob_0) * (c_dot)) - ([r_J2_M ; r_C2_M]) == 0 

% (IGNORING OMNI WHEEL CONSTRAINTS AND ONE OF THE SLIDIING CONSTRAINTS OF FIXED WHEEL)
Reduced_J_M = [ FW1_RCM ; FW2_RCM];
Reduced_C_M = [FW1_SCM];
Reduced_JC = [ Reduced_J_M ; Reduced_C_M];
reduced_r_J2_M = [ rf*wL; rf*wr];
reduced_r_C2_M = [0];

% THE REDUCED EQUATION 
Reduced_Equation = (Reduced_JC) * ((R_rob_0) * (c_dot)) - ([rf*wL; rf*wr; 0]) == 0

%% ANSWER 9 (VELOCITY KINEMATICS OF MOBILE ROBOT BASE)
c_dot = simplify((R_0_rob) * (inv(Reduced_JC)) * ([rf*wL; rf*wr; 0]))
disp('units : m/sec and rad/sec ') 

%% ANSWER 10 (Numeric Velocity Kinematics for Mobile Robot Base) 
c_dot_ans10 = vpa(simplify(subs(c_dot, [alp1,beta1,alp2,beta2,a,wL,wr,thet,rf], [(pi/180)*90,0,-(pi/180)*90,pi,.507,pi,2*pi,(pi/180)*30,0.143])))  
disp('units : m/sec and rad/sec ') 
%% ANSWER 11
% The tranformation matrix from world to robot frame
T_W_R = [ cos(thet), -sin(thet), 0, 2.5 ; sin(thet), cos(thet), 0, 1.5 ; 0, 0, 1, 0 ; 0, 0, 0, 1];
T_W_R_numeric = vpa(simplify(subs(T_W_R, [thet], [(pi/180)*30])));
% The transformation matrix from World to tip frame (T_W_R * T_R_T)
T_W_T = vpa(simplify((T_W_R) * (T_0_5)))

%% ANSWER 12 (Numeric answer for Transformation from world to tip frame).   
T_W_T_ans12 = vpa(simplify((T_W_R_numeric) * (T_0_5_ans5a)))
T_W_T_tip_Position = [T_W_T_ans12(1:3,4)]
disp('units: m')
rot_angle = converttoEuler(T_W_T_ans12)
disp('units: degree')
%% ANSWER 13 (Combined Velocity Kinematics) 
R_special = [cos(thet),-sin(thet),0,0,0,0 ;sin(thet),cos(thet),0,0,0,0 ;0,0,1,0,0,0; 0,0,0,cos(thet),-sin(thet),0; 0,0,0,sin(thet),cos(thet),0; 0,0,0,0,0,1];
q_dot_R = [x_qR; y_qR; z_qR; wL; wr]; 
new_col = [(rf/2);0;0;0;0;(-rf/a)];
new_col2 = [(rf/2);0;0;0;0;(rf/a)];
J_special = [J , new_col , new_col2];
R_J_special = (R_special * J_special); 

X_dot_R = vpa(simplify((R_J_special) * q_dot_R ))
disp('units : m/sec and rad/sec ') 

%% ANSWER 14 (FORCE PROPAGATION)
F_tip_R = [30*cos((pi/180)*30);30*sin((pi/180)*30);0;0;0;0];
J_special_R = vpa(simplify(subs(R_J_special, [q1,q2,q3,b,c,d,e,f,g,wL,wr,rf,a,thet], [(pi/180)*20,(pi/180)*90,(pi/180)*30,.424,.300,.380,.328,.323,.0824,pi,2*pi,.143,.507,(pi/180)*30])));

Tau_R = vpa(simplify( transpose(J_special_R) * F_tip_R))
disp('units : N-m ') 
Tau_Wheels = [Tau_R(4,1); Tau_R(5,1)]
disp('units : N-m ') 

% A function for calculating the Tranformations for answer 3
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














