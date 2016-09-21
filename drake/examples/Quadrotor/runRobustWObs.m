function [p, xtraj, utraj, prog] = runRobustWObs

% simple planning demo which takes the quadrotor from hover at y=0m to a new hover at
% y=10m with minimal thrust. Creates a forest set up to have the Quadrotor
% swerve between trees blocking its path to get to y=10m.

p = Quadrotor();

% The important trees to create swerving path
p = addTree(p, [.5,.5,1.5], [0;0], pi/4);
trees = [0;0];
% p = addTree(p, [.8,.45,1.25], [.20;2.5], pi/4);
% p = addTree(p, [.5,.35,1.65], [-.25;5], -pi/6);
% p = addTree(p, [.55,.65,1.5], [.25;7.5], pi/4);
% p = addTree(p, [.55,.85,1.6], [-1.35;8.5], pi/3.7);
% p = addTree(p, [.85,.95,1.65], [-1.85;5.2], -pi/3.7);
% p = addTree(p, [.75,.9,1.75], [2;4.4], -pi/5);
% trees = [[.20;2.5], [-.25;5], [.25;7.5], [-1.35;8.5], [-1.85;5.2], [2;4.4]];

% Random trees to make forest bigger and more dense
%p = addTrees(p, 25);

% initial conditions
x0 = Point(getStateFrame(p));
x0.base_y = -2;
x0.base_z = .5;
x0.base_x = .5;
x0.base_ydot = 5;
u0 = double(nominalThrust(p));

% final conditions: translated in x
xf = x0;
xf.base_y = 2;
xg = double(xf);

v = constructVisualizer(p); %, struct('use_collision_geometry',true));
v.draw(0,double(x0));

N = 31;
tf0 = 1.25; % duration
duration = [tf0 tf0];
traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));
traj_init.u = ConstantTrajectory(u0);

options.integration_method = DirtranTrajectoryOptimization.MIDPOINT;
prog1 = DirtranTrajectoryOptimization(p, N, duration,options);

prog1 = prog1.addStateConstraint(ConstantConstraint(double(x0)),1);
prog1 = prog1.addInputConstraint(ConstantConstraint(u0),1);
prog1 = prog1.addStateConstraint(ConstantConstraint(double(xf)),N);
prog1 = prog1.addInputConstraint(ConstantConstraint(u0),N);

prog1 = prog1.addRunningCost(@cost);

collision_constraint = FunctionHandleConstraint(1.2*ones(1,1),inf*ones(1,1),2,@treeDistance);
collision_constraint.grad_method = 'user';
collision_constraint.grad_level = 1;
prog1 = prog1.addStateConstraint(collision_constraint,2:N-1,1:2);

prog1 = prog1.setSolverOptions('snopt','majoroptimalitytolerance',1e-3);
prog1 = prog1.setSolverOptions('snopt','majorfeasibilitytolerance',1e-3);
prog1 = prog1.setSolverOptions('snopt','minorfeasibilitytolerance',1e-3);

prog1 = addPlanVisualizer(p,prog1);
tic
[xtraj1,utraj1,z1,F1,info1] = prog1.solveTraj(tf0,traj_init);
toc

v.playback(xtraj1); %, struct('slider',true));

function [g,dg] = cost(h,x,u)

Q = h*.1*eye(12);
R = h*1*eye(4);
g = .5*(u'*R*u + (x-xg)'*Q*(x-xg));
dg = [g, (x-xg)'*Q, u'*R];

end

function [c, dc] = treeDistance(q)
    c = zeros(size(trees,2),1);
    dc = zeros(size(trees,2),2);
    for k = 1:length(c)
        c(k) = (q(1:2)-trees(1:2,k))'*(q(1:2)-trees(1:2,k));
        dc(k,:) = 2*(q(1:2)-trees(1:2,k))';
    end
end

D = .2^2*eye(3);
E0 = .01^2*eye(12);
Q = blkdiag(10*eye(6), 1*eye(6));
R = .1*eye(4);
Qf = Q;
options.integration_method = DirtranTrajectoryOptimization.MIDPOINT;
prog2 = RobustDirtranTrajectoryOptimization(p,N,D,E0,Q,R,Qf,duration,options);

prog2 = prog2.addRobustCost(Q,10*R,Qf);
%prog2 = prog2.addRobustInputConstraint();
prog2 = prog2.addRobustStateConstraint(collision_constraint,2:(N-1),1:2);

prog2 = addPlanVisualizer(p,prog2);

prog2 = prog2.addStateConstraint(ConstantConstraint(double(x0)),1);
prog2 = prog2.addInputConstraint(ConstantConstraint(u0),1);
prog2 = prog2.addStateConstraint(ConstantConstraint(double(xf)),N);
prog2 = prog2.addInputConstraint(ConstantConstraint(u0),N-1);

prog2 = prog2.addStateConstraint(collision_constraint,2:N-1,1:2);

prog2 = prog2.addRunningCost(@cost);

prog2 = prog2.setSolverOptions('snopt','majoroptimalitytolerance',1e-2);
prog2 = prog2.setSolverOptions('snopt','majorfeasibilitytolerance',1e-3);
prog2 = prog2.setSolverOptions('snopt','minorfeasibilitytolerance',1e-3);

tic
[xtraj2,utraj2,z2,F2,info2] = prog2.solveTraj(tf0,traj_init);
toc

v.playback(xtraj2); %, struct('slider',true));

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



end



