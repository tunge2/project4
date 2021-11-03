%
% invkin_iterJ.m
%
% inverse kinematics using Jacobian iteration (3D)


function robot = invkin_iterJ_3D(robot)
    %input parameters
    N=robot.MaxIter;
    alpha=robot.StepSize;
    w=robot.Weights;
    
    %Final Target R and p
    R0Td=robot.T(1:3,1:3);
    p0Td=robot.T(1:3,4);
    
    q0=robot.q;
    n=length(q0);
    q=zeros(n,N+1);
    q(:,1)=q0;
    p0T=zeros(3,N+1);
    
    for i=1:N
        robot.q=q(:,i);
        robot=fwdkiniter(robot);
        p0T(:,i)=robot.T(1:3,4);
        %s=R2qv(robot.T(1:3,1:3)*R0Td');
        quat=R2q(robot.T(1:3,1:3)*R0Td');
        %quat(1)*quat(2:4)
        %[k, th]=R2kth(robot.T(1:3,1:3)*R0Td');
        %k*theta
        dX=[quat(1)*quat(2:4);p0T(:,i)-p0Td];
        % Jacobian update
        qq = q(:,i) - alpha*robot.J'*inv(robot.J*robot.J' + .01*diag(1./w))*dX;
        q(:,i+1) = (qq>pi).*(-2*pi+qq) + (qq<-pi).*(2*pi+qq) + (qq<pi).*(qq>-pi).*qq;
    end
    
    %fwd kinematics with final q
    robot.q=q(:,N+1);
    robot=fwdkiniter(robot);
end 


% function robot = invkin_iterJ_3D(robot,N,alpha,reg_term)
%     % N is number of iterations
%     % alpha is step size
%     
%     % define unit z vector
%     % target (R,p)
%     posedez = robot.T;
%     Rd = robot.T(1:3,1:3);
%     pTd = robot.T(1:3,4);
% 
%     % set up storage space
%     q0=robot.q;
%     n=length(q0); % # of joints
%     q=zeros(n,N+1);
%     q(:,1)=q0; % output joint displacements
%     pT = zeros(3,N+1); % output p
%     qT = zeros(3,N+1); % output quaternion
%     % first iteration
%     robot.q = q(:,1);
%     robot = fwdkiniter(robot);
%     % iterative update
% %     boo = true;
%     for i=1:N
%         % forward kinematics
%         
%         R0T = robot.T(1:3,1:3);
%         s = R2qv(R0T * Rd');
%         qT(:,i) = s(:,1);
%         pT(:,i)=robot.T(1:3,4);
%         % task space error: angular error and position error
%         w = [50;50;50;1;1;1];
%         %dX = [qT(:,i); 1000*(pT(:,i) - pTd)].*robot.Weights;
%         [k, th] = R2kth(robot.T(1:3,1:3)*Rd');
%         dX = [k*th;pT(:,i)-pTd];
%         reg_term = reg_term+1;
%         %qq = q(:,i) - alpha * robot.J' * inv(robot.J*robot.J' + reg_term * eye(6))*dX;
%         qq = q(:,i) - alpha * robot.J' * inv(robot.J*robot.J'+0.01*diag(1./w))*dX;
%         q(:,i+1) = (qq>pi).*(-2*pi+qq) + (qq<-pi).*(2*pi+qq) + (qq<pi).*(qq>-pi).*qq;
%         
%         
%         robot.q = q(:,i+1);
%         robot = fwdkiniter(robot);
%         if norm(posedez - robot.T,'fro') < 0
%             break
%         end
%     end
%     % final iteration
%     robot.q = q(:,i+1); 
% end
% 
function qv=R2qv(R)
    q01 = (1/2)*sqrt(trace(R)+1);
    q02 = -(1/2)*sqrt(trace(R)+1);
    if abs(q01)<1e-5
        [k, ~] = R2kth(R);
        qv = k;
    else
        qv(:,1)=vee(R-R')/4/q01;
        qv(:,2)=vee(R-R')/4/q02;
    end
end


function q=R2q(R)
  
  q=zeros(4,1);
  q(1)=.5*sqrt(trace(R)+1);
  if abs(q(1))<1e-5
    [k,~]=R2kth(R);
    q(2:4)=k;
  else
    q(2:4)=vee(R-R')/4/q(1);
  end
  
end









