function [p, xtraj, utraj, prog] = runRobustWObs

% simple planning demo which takes the quadrotor from hover at y=0m to a new hover at
% y=10m with minimal thrust. Creates a forest set up to have the Quadrotor
% swerve between trees blocking its path to get to y=10m.

p = Quadrotor();

% The important trees to create swerving path
p = addTree(p, [.8,.45,1.25], [.20;2.5], pi/4);
p = addTree(p, [.5,.35,1.65], [-.25;5], -pi/6);
p = addTree(p, [.55,.65,1.5], [.25;7.5], pi/4);
p = addTree(p, [.55,.85,1.6], [-1.35;8.5], pi/3.7);
p = addTree(p, [.85,.95,1.65], [-1.85;5.2], -pi/3.7);
p = addTree(p, [.75,.9,1.75], [2;4.4], -pi/5);

% Random trees to make forest bigger and more dense
%p = addTrees(p, 25);

% initial conditions: all-zeros
x0 = Point(getStateFrame(p));
x0.base_y = -1.5;
x0.base_z = .5;
x0.base_ydot = 5;
u0 = double(nominalThrust(p));

% final conditions: translated in x
xf = x0;
xf.base_y = 11;

v = constructVisualizer(p); %, struct('use_collision_geometry',true));
v.draw(0,double(x0));

tf0 = 4; % initial guess at duration
traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));
traj_init.u = ConstantTrajectory(u0);

N = 41;
duration = [0 4];

%prog = DirtranTrajectoryOptimization(p, N, duration);

D = 1/(.25^2)*eye(3);
Q = blkdiag(10*eye(6), 1*eye(6));
R = .1*eye(4);
Qf = Q;
prog = RobustDirtranTrajectoryOptimization(p,N,D,Q,R,Qf,duration);
prog = prog.addRobustCost(Q,R,Qf);
prog = prog.addRobustInputConstraint();

prog = addPlanVisualizer(p,prog);

prog = prog.addStateConstraint(ConstantConstraint(double(x0)),1);
prog = prog.addInputConstraint(ConstantConstraint(u0),1);

prog = prog.addStateConstraint(ConstantConstraint(double(xf)),N);
prog = prog.addInputConstraint(ConstantConstraint(u0),N-1);

prog = prog.addRunningCost(@cost);
prog = prog.addFinalCost(@final_cost);

collision_constraint = generateConstraint(MinDistanceConstraint(p,0.1),0);
prog = prog.addStateConstraint(collision_constraint{1},1:N,1:getNumPositions(p));

prog = prog.setSolverOptions('snopt','majoroptimalitytolerance',1e-3);
prog = prog.setSolverOptions('snopt','majorfeasibilitytolerance',1e-3);

tic
[xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
toc

v.playback(xtraj, struct('slider',true));

end

function [g,dg] = cost(dt,x,u)

R = eye(4);
g = u'*R*u;
%g = sum((R*u).*u,1);
dg = [zeros(1,1+size(x,1)),2*u'*R];
%dg = zeros(1, 1 + size(x,1) + size(u,1));

end

function [h,dh] = final_cost(t,x)

h = t;
dh = [1,zeros(1,size(x,1))];

end

