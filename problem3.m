
% forward/inverse fkinematics checker for an elbow arm
% using MATLAB inversekinematics solver vs. subproblem decomposition
% 
% make sure fwdkin_example_rbt is in the path to generate the elbow_rbt
% rigid body tree
% 
% output plots compare EE solution accuracy and computation times

clear all; close all;

zz=zeros(3,1); ex = [1;0;0]; ey = [0;1;0]; ez = [0;0;1];
L1=399.1;
L2=448; %L2=350;
L3=42;
L4=451; %L4=351;
L5=82;

p01=0*ex+L1*ez;p12=zz;p23=L2*ez;p34=L3*ez+L4*ex;p45=zz;p56=zz;p6T=L5*ex;
h1=ez;h2=ey;h3=ey;h4=ex;h5=ey;h6=ex;
irb1200.P=[p01 p12 p23 p34 p45 p56 p6T];
irb1200.H=[h1 h2 h3 h4 h5 h6];
irb1200.joint_type=[0 0 0 0 0 0];
radius=.01;
[irb1200_rbt,colLink]=defineRobot(irb1200,radius);

% MATLAB robotics toolbox inverse kinematics object
ik = inverseKinematics('RigidBodyTree',irb1200_rbt);

% set storage
N=50;

n=length(irb1200.H);
q=zeros(n,N);
qsol=zeros(n,N);
qsol1=zeros(n,8,N);
qsol2=zeros(n,N);
T=zeros(4,4,N);
Tsol=zeros(4,4,N);
Tsol1=zeros(4,4,8,N);
Jsol2=zeros(6,6,N);
errT=zeros(N,1);
errT1=zeros(N,8);
errT2=zeros(N,1);
telapsed=zeros(1,N);
telapsed1=zeros(1,N);
telapsed2=zeros(1,N);

for i=1:N
    % random arm configuration
    q(:,i)=(rand(6,1)-.5)*pi;
    irb1200.q=q(:,i);
    % forward kinematics
    irb1200=fwddiffkiniter(irb1200);
    % nominal end effector pose
    T(:,:,i)=irb1200.T;
    % MATLAB's inverse kinematics
    % keep track of the computation time
    tstart=tic;
    [qsol(:,i),solnInfo]=...
        ik('body7',irb1200.T,ones(1,6),q(:,i)+q(:,i)*.1*randn);
    telapsed(i) = toc(tstart);    
    % forward kinematics again to compare with EE pose
    irb1200.q=qsol(:,i);
    irb1200=fwddiffkiniter(irb1200);
    Tsol(:,:,i)=irb1200.T;
    % find EE pose error
    errT(i)=norm(T(:,:,i)-Tsol(:,:,i),'fro');   
    % now use our own exact inverse kinematics
    irb1200.q=q(:,i);
    irb1200.T=T(:,:,i);
    % keep track of the computation time
    tstart=tic;
    irb1200=invkinelbow(irb1200);
    telapsed1(i) = toc(tstart);
    % there are 8 solutions
    qsol1(:,:,i)=irb1200.q;
    % find all the EE pose
    for j=1:8
        irb1200.q=qsol1(:,j,i);
        irb1200=fwddiffkiniter(irb1200);
        Tsol1(:,:,j,i)=irb1200.T;
        Jsol1(:,:,j,i)=irb1200.J;
        errT1(i,j)=norm(T(:,:,i)-Tsol1(:,:,j,i),'fro');
    end
    % our own iterative inverse kinematics
    irb1200.MaxIter=100;
    irb1200.StepSize=.4;
    weight1 = 100;
    weight2 = 500;
    irb1200.Weights=[weight1;weight1;weight1;weight2;weight2;weight2];
    irb1200.q=q(:,i)+q(:,i).*.1.*randn(6,1);
    irb1200.T=T(:,:,i);
    tstart=tic;
    irb1200=invkin_iterJ_3D(irb1200);
    telapsed2(i) = toc(tstart);
    qsol2(:,i)=irb1200.q;
    irb1200=fwddiffkiniter(irb1200);
    Tsol2(:,:,i)=irb1200.T;
    Jsol2(:,:,i)=irb1200.J;
    errT2(i)=norm(T(:,:,i)-Tsol2(:,:,i),'fro');
end 

figure(10);semilogy((1:N),errT,'bx',(1:N),errT2,'linewidth',2); hold on;
for j=1:8
    plot((1:N),errT1(:,j),'ro','linewidth',2);
end
grid on
hold off;
xlabel('random test number');ylabel('end effector error');
title('End effector pose error || T - T_{solve} ||_F');
legend('MATLAB','iterative Jacobian','exact by subproblem');

fprintf('max MATLAB EE error: %g \n',max(errT));
fprintf('max subproblem EE error: %g \n',max(max(errT1)));
fprintf('max subproblem EE error: %g \n',max(max(errT2)));

figure(20);semilogy((1:N),telapsed,'bx',(1:N),telapsed1,'ro',...,
    (1:N),telapsed2,'g^','linewidth',2);
grid on
xlabel('random test number');ylabel('elapsed time (sec)');
title('inverse kinematics computation time');
legend('MATLAB','exact by subproblem','iterative Jacobian');

fprintf('max and average MATLAB invkin computation time: %g, %g \n',...
    max(telapsed),mean(telapsed));
fprintf('max and average subproblem invkin computation time: %g, %g \n',...
    max(telapsed1),mean(telapsed1));
fprintf('max and average iterative invkin computation time: %g, %g \n',...
    max(telapsed2),mean(telapsed2));


