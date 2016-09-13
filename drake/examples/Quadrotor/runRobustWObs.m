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

N = 31;
duration = [0 4];

prog1 = DirtranTrajectoryOptimization(p, N, duration);

D = .01^2*eye(3);
E0 = .001^2*eye(12);
Q = blkdiag(10*eye(6), 1*eye(6));
R = .1*eye(4);
Qf = Q;
prog2 = RobustDirtranTrajectoryOptimization(p,N,D,E0,Q,R,Qf,duration);
prog2 = prog2.addRobustCost(Q,R,Qf);

% prog1 = addPlanVisualizer(p,prog1);
prog2 = addPlanVisualizer(p,prog2);

prog1 = prog1.addStateConstraint(ConstantConstraint(double(x0)),1);
prog1 = prog1.addInputConstraint(ConstantConstraint(u0),1);
prog2 = prog2.addStateConstraint(ConstantConstraint(double(x0)),1);
prog2 = prog2.addInputConstraint(ConstantConstraint(u0),1);

prog1 = prog1.addStateConstraint(ConstantConstraint(double(xf)),N);
prog1 = prog1.addInputConstraint(ConstantConstraint(u0),N-1);
prog2 = prog2.addStateConstraint(ConstantConstraint(double(xf)),N);
prog2 = prog2.addInputConstraint(ConstantConstraint(u0),N-1);
%prog2 = prog2.addRobustInputConstraint();

prog1 = prog1.addRunningCost(@cost);
prog1 = prog1.addFinalCost(@final_cost);
prog2 = prog2.addRunningCost(@cost);
prog2 = prog2.addFinalCost(@final_cost);

collision_constraint = generateConstraint(MinDistanceConstraint(p,0.1),0);
prog1 = prog1.addStateConstraint(collision_constraint{1},1:N,1:getNumPositions(p));
prog2 = prog2.addStateConstraint(collision_constraint{1},1:N,1:getNumPositions(p));
prog2 = prog2.addRobustStateConstraint(collision_constraint{1},1:getNumPositions(p));

prog1 = prog1.setSolverOptions('snopt','majoroptimalitytolerance',1e-2);
prog1 = prog1.setSolverOptions('snopt','majorfeasibilitytolerance',1e-3);
prog2 = prog2.setSolverOptions('snopt','majoroptimalitytolerance',1e-2);
prog2 = prog2.setSolverOptions('snopt','majorfeasibilitytolerance',1e-3);

tic
[xtraj1,utraj1,z1,F1,info1] = prog1.solveTraj(tf0,traj_init);
toc

tic
[xtraj2,utraj2,z2,F2,info2] = prog2.solveTraj(tf0,traj_init);
toc

tsamp1 = xtraj1.tspan(1):.1:xtraj1.tspan(2);
xsamp1 = xtraj1.eval(tsamp1);
usamp1 = utraj1.eval(tsamp1);

tsamp2 = xtraj2.tspan(1):.1:xtraj2.tspan(2);
xsamp2 = xtraj2.eval(tsamp2);
usamp2 = utraj2.eval(tsamp2);

figure();
subplot(3,1,1);
plot(tsamp1, xsamp1(1,:));
hold on;
plot(tsamp2, xsamp2(1,:));
subplot(3,1,2);
plot(tsamp1, xsamp1(2,:));
hold on;
plot(tsamp2, xsamp2(2,:));
subplot(3,1,3);
plot(tsamp1, xsamp1(3,:));
hold on;
plot(tsamp2, xsamp2(3,:));

figure();
subplot(4,1,1);
plot(tsamp1, usamp1(1,:));
hold on;
plot(tsamp2, usamp2(1,:));
subplot(4,1,2);
plot(tsamp1, usamp1(2,:));
hold on;
plot(tsamp2, usamp2(2,:));
subplot(4,1,3);
plot(tsamp1, usamp1(3,:));
hold on;
plot(tsamp2, usamp2(3,:));
subplot(4,1,4);
plot(tsamp1, usamp1(4,:));
hold on;
plot(tsamp2, usamp2(4,:));

%Simulate


v.playback(xtraj1); %, struct('slider',true));
v.playback(xtraj2); %, struct('slider',true));

end

function [g,dg] = cost(h,x,u)

R = .1*eye(4);
g = h*(u'*R*u);
dg = [u'*R*u, zeros(1,12), 2*h*u'*R];

end

function [h,dh] = final_cost(t,x)

h = t;
dh = [1,zeros(1,size(x,1))];

end

