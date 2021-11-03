clear all
close all

syms L1 L2 L3 L4 L5 positive
syms q1 q2 q3 q4 q5 q6 real

ex = [1 0 0]';ey = [0 1 0]';ez = [0 0 1]';zz = zeros(3,1);
h1=ez;h2=ey;h3=ey;h4=ex;h5=ey;h6=ex;
p01=L1*ez;p12=zz;p23=L2*ez;p34=L3*ez+L4*ex;p45=zz;p56=zz;p6T=L5*ex;

R01 = rot(h1,q1);
R12 = rot(h2,q2);
R23 = rot(h3,q3);
R34 = rot(h4,q4);
R45 = rot(h5,q5);
R56 = rot(h6,q6);

R03 = R01*R12*R23;
R13 = R12*R23;
R35 = R34*R45;
R36 = R35*R56;

R30 = R03';
R31 = R13';
R32 = R23';
R0T = R01*R12*R23*R34*R45*R56;

p14_3 = R31*p12 + R32*p23 + p34;
p24_3 = R32*p23 + p34;
p4T_3 = R34*p45 + R35*p56 + R36*p6T;

J4_3_eric = [[R31*h1; crossmat(R31*h1)*p14_3], [R32*h2; crossmat(R32*h2)*p24_3], [h3; crossmat(h3)*p34],[R34*h4; zz], [R35*h5; zz], [R36*h6; zz]];

J4_3_eric = simplify(expand(J4_3_eric));

JT_0_eric=phi(R03,R30*R0T*p6T)*J4_3_eric; %
JT_0_eric = simplify(JT_0_eric);

detJ4_3 = simplify(det(J4_3_eric));
detJT_0 = simplify(det(JT_0_eric));
cccheck = detJ4_3-detJT_0;

%detJ4_3 is just L2*sin(q5)*(L4*cos(q3) + L3*sin(q3))*(L4*cos(q2 + q3) + L3*sin(q2 + q3) + L2*sin(q2))

% determinant of J4_3_eric = determinant of JT_0_eric 

%A = L2*sin(q5)*((L4^2*cos(q2 + 2*q3))/2 - (L3^2*cos(q2 + 2*q3))/2 + (L3^2*cos(q2))/2 + (L4^2*cos(q2))/2 + (L2*L4*sin(q2 - q3))/2 + L3*L4*sin(q2 + 2*q3) - (L2*L3*cos(q2 + q3))/2 + (L2*L4*sin(q2 + q3))/2 + (L2*L3*cos(q2 - q3))/2)
%B = L2*sin(q5)*((L4^2*cos(q2 + 2*q3))/2 - (L3^2*cos(q2 + 2*q3))/2 + (L3^2*cos(q2))/2 + (L4^2*cos(q2))/2 + (L2*L4*sin(q2 - q3))/2 + L3*L4*sin(q2 + 2*q3) - (L2*L3*cos(q2 + q3))/2 + (L2*L4*sin(q2 + q3))/2 + (L2*L3*cos(q2 - q3))/2)














