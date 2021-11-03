%
% comparing iterative Jacobian calculation in class vs. MATLAB
%

N=500;

for i=1:N
    
% random joint angles
q=(rand(6,1)-.5)*2*pi;

% iterative Jacobian calculation
irb1200.q=q;
tstart=tic;
irb1200=fwdkiniter(irb1200);
t(i)=toc(tstart);
J=irb1200.J;

% MATLAB function
tstart=tic;
J1=geometricJacobian(irb1200_rbt,q,'body7');
t1(i)=toc(tstart);

% check difference
diffJ1(i)=(norm(J-J1));

% using symbolic expression for Jacobian computation
R0T=irb1200.T(1:3,1:3);
tstart=tic;
J4_3_num=J4_3func(q);
R03=rot(irb1200.H(:,1),q(1))*rot(irb1200.H(:,2),q(2))*rot(irb1200.H(:,3),q(3));
R3T=rot(irb1200.H(:,4),q(4))*rot(irb1200.H(:,5),q(5))*rot(irb1200.H(:,6),q(6));
JT_0=phi(R03,R3T*irb1200.P(:,end))*J4_3_num;
t2(i)=toc(tstart);

diffJ2(i)=norm(J-JT_0);

end

figure(10);semilogy((1:N),t,'^',(1:N),t1,'o',(1:N),t2,'x','linewidth',2);
xlabel('run #');ylabel('sec');
grid on
title('Jacobian computation time');
legend('Iterative method','MATLAB','symbolic');



function J=J4_3func(q)
L2=448;L3=42;L4=451;
J=[-sin(q(2) + q(3)),0, 0, 1, 0, cos(q(5));
    0, 1,   1, 0, cos(q(4)),  sin(q(4))*sin(q(5));
   cos(q(2) + q(3)), 0,   0, 0, sin(q(4)), -cos(q(4))*sin(q(5));
   0, L3 + L2*cos(q(3)),  L3, 0,       0, 0;
L4*cos(q(2) + q(3)) + L3*sin(q(2) + q(3)) + L2*sin(q(2)), 0, 0, 0,  0, 0;
   0, L2*sin(q(3)) - L4, -L4, 0, 0, 0];
end

% function J=test()
% syms q2 q3 q4 q5 L2 L3 L4 real
% J = [-sin(q2 + q3),0, 0, 1, 0, cos(q5); 0, 1, 1, 0, cos(q4),  sin(q4)*sin(q5);cos(q2 + q3), 0,   0, 0, sin(q4), -cos(q4)*sin(q5);0, L3 + L2*cos(q3),  L3, 0, 0, 0; L4*cos(q2 + q3) + L3*sin(q2 + q3) + L2*sin(q2), 0, 0, 0,  0, 0; 0, L2*sin(q3) - L4, -L4, 0, 0, 0];
% end














