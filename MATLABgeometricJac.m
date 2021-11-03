ex = [1 0 0]'; ey = [0 1 0]'; ez = [0 0 1]'; zz = zeros(3,1);
L1=399.1;L2=448;L3=42;L4=451;L5=82;
p01=L1*ez;p12=zz;p23=L2*ez;p34=L3*ez+L4*ex;p45=zz;p56=zz;p6T=L5*ex;
h1=ez; h2=ey; h3=ey; h4=ex; h5=ey; h6=ex;

irb1200.P=[p01 p12 p23 p34 p45 p56 p6T];
irb1200.H=[h1 h2 h3 h4 h5 h6];
irb1200.joint_type=[0 0 0 0 0 0];
radius=0.1;
[irb1200_rbt,~]=defineRobot(irb1200,radius);

q=(rand(6,1)-.5)*2*pi;
J_MATLAB = geometricJacobian(irb1200_rbt,q,'body7');

