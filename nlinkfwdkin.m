
function robot = nlinkfwdkin(robot)

    % Extract # of links
    n=length(robot.q);
    
    % ********** calculate T_{0T} ***********
    % orientation of end effector
    qT=sum(robot.q);
    % for position, start with the first link p01 (stationary)
    p0T=robot.P(1:2,1);
    % calculate R_{0i} p_{i,i+1}, i=1,..,n
    for i=1:n
        P(:,i)=rot2(sum(robot.q(1:i)))*robot.P(1:2,i+1);
    end
    % p0T = sum_{i=0}^n R_{0i} p_{i,i+1},
    p0T=p0T+sum(P,2);
    
    % output the homogeneous matrix
    %robot.T=eye(4,4);
    robot.T(1:2,4)=p0T;
    robot.T(1:2,1:2)=rot2(qT);
    robot.T(3,3)=1;
    robot.T(4,4)=1;
    
    %
    % ********** calculate J_T ***********
    % 
    % set up cross product 
    ez_cross=[0 -1;1 0];
    % ith column of J is [1;ez x \sum_{k=i}^n R_{0k} p_{k,k+1}]
    for i=1:n
        robot.J(2:3,i)=ez_cross*sum(P(:,i:n),2);
        robot.J(1,i)=1;
    end

end
