%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% a %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fwdkiniter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% b %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

ex = [1 0 0]';ey = [0 1 0]';ez = [0 0 1]';zz = zeros(3,1);

L1=399.1;L2=448;L3=42;L4=451;L5=82;

p01=L1*ez;p12=zz;p23=L2*ez;p34=L3*ez+L4*ex;p45=zz;p56=zz;p6T=L5*ex;

h1=ez;h2=ey;h3=ey;h4=ex;h5=ey;h6=ex;

syms L1 L2 L3 L4 L5 positive
syms q1 q2 q3 q4 q5 q6 real

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

J4_3_eric = simplify(J4_3_eric);
%pretty(J4_3_eric);

JT_0_eric=phi(R03,R30*R0T*p6T)*J4_3_eric; %
JT_0_eric = simplify(expand(JT_0_eric));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% c/d %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% i/ii %%%%%%%%%%%%%%%%%%%%%%%

irb1200.P=[p01 p12 p23 p34 p45 p56 p6T];

irb1200.H=[h1 h2 h3 h4 h5 h6];
irb1200.joint_type=[0 0 0 0 0 0];
radius=0.1;
[irb1200_rbt,~]=defineRobot(irb1200,radius);

N = 100;

errorA = zeros(1,N);
errorB = errorA;

for i=1:N
    %generate random joint angles
    q=(rand(6,1)-.5)*2*pi;
    
    % iterative Jacobian
    tstart = tic;
    irb1200.q=q;
    irb1200=fwdkiniter(irb1200);
    t_iter(i) = toc(tstart);
    J_iter = irb1200.J;
    
    % MATLAB function
    tstart = tic;
    J_MATLAB = geometricJacobian(irb1200_rbt,q,'body7');
    t_MATLAB(i) = toc(tstart);
    
    % Symbolic expression for Jacobian
    tstart = tic;
    R0T=irb1200.T(1:3,1:3);
    J4_3_num=makeJ4_3(q);
    R03=rot(irb1200.H(:,1),q(1))*rot(irb1200.H(:,2),q(2))*rot(irb1200.H(:,3),q(3));
    R3T=rot(irb1200.H(:,4),q(4))*rot(irb1200.H(:,5),q(5))*rot(irb1200.H(:,6),q(6));
    J_syms = phi(R03,R3T*irb1200.P(:,end))*J4_3_num;
    t_syms(i) = toc(tstart);
    
    % Compute difference
    errorA(i) = norm(J_iter-J_syms);
    errorB(i) = norm(J_MATLAB-J_syms);
end

figure(1);
semilogy((1:N),t_iter,'^',(1:N),t_MATLAB,'o',(1:N),t_syms,'x','linewidth',2);
xlabel('Run #');ylabel('Time (sec)');
grid on
title('Jacobian computation time');
legend('Iterative','MATLAB','Symbolic');

figure(2);
semilogy((1:N),errorA,'^',(1:N),errorB,'o','linewidth',2);
xlabel('Run #');ylabel('norm(J-J_{syms})');
grid on
title('Jacobian Error - Euclidean Norm');
legend('Iterative','MATLAB');

function J=makeJ4_3(q)
L2=448; %L2=350;
L3=42;
L4=451; %L4=351;

J=[-sin(q(2) + q(3)),0, 0, 1, 0, cos(q(5));
    0, 1,   1, 0, cos(q(4)),  sin(q(4))*sin(q(5));
    cos(q(2) + q(3)), 0,   0, 0, sin(q(4)), -cos(q(4))*sin(q(5));
    0, L3 + L2*cos(q(3)),  L3, 0,       0, 0;
    L4*cos(q(2) + q(3)) + L3*sin(q(2) + q(3)) + L2*sin(q(2)), 0, 0, 0,  0, 0;
    0, L2*sin(q(3)) - L4, -L4, 0, 0, 0];
end
