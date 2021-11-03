% fwdkiniter.m
% 
% forward kinematics using POE approach
%
% input:    robot = robot object with
%                   .q = joint displacement vector
%                   .H = [h1, ..., hn] (3xn)
%                   .P = [p01, ..., p_{n-1,n},p_{n,T}] (3x(n+1))
%                   .joint_type = 1xn (0 for revolute, 1 for prismatic)
%
% output:   robot = robot object with 
%                   .T = homogeneous transformation T_{0T}
%                   .J = Jacobian J_T in the 0 frame
%
% usage:
%       robot = fwddiffkiniter(robot);   
%
% 

function robot=fwdkiniter(robot)
q=robot.q;
n=length(robot.q)-2; %THIS IS JSUT TO MAKE PROJ3 work
P = robot.P;
T=eye(4,4);
J=zeros(6,n);

for i=1:n
    h=robot.H(1:3,i);
    R=expm(hat(h)*q(i)); %2 loops?
    p=P(1:3,i);
    
    J = phi(eye(3,3),T(1:3,1:3)*p)*J; %phi {i,i-1} with oldJ
    J(:,i) = [T(1:3, 1:3)*h; zeros(3,1)]; % inserting revolute joint
    T=T*[R p;zeros(1,3) 1];
end    

robot.J = phi(eye(3,3), T(1:3,1:3)*P(:,n+1))*J; %computes J_T from J_6
robot.T=T*[eye(3,3) P(:,n+1); zeros(1,3) 1];
end

function phimat = phi(R,p)
    phimat = [R zeros(3,3); -R*hat(p) R];
end

% cross production function
function khat = hat(k)
  
  khat=[0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
  
end
