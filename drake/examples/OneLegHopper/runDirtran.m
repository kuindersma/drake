function runDirtran(N)

warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits');
options.twoD = true;
options.view = 'right';
options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.ignore_self_collisions = true;
s = 'OneLegHopper.urdf';
dt = 0.001;
r = TimeSteppingRigidBodyManipulator(s,dt,options);
r = r.setStateFrame(OneLegHopperState(r));
r = r.setOutputFrame(OneLegHopperState(r));

nq = r.getNumPositions;
nx = r.getNumStates;
nu = r.getNumInputs;
nc = r.getNumContactPairs;

v = r.constructVisualizer();

tf0 = 1.0;

q0 = [0;0;.6;-1.2;.6+pi/2];
x0 = [q0;0*q0];
x0 = double(r.resolveConstraints(x0));

xG =  x0; 
xG(1) = x0(1) - 0.2;


traj_opt = SmoothContactImplicitTrajectoryOptimization(r,N,tf0*[(1-0.2) (1+0.2)],options);
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0),1);
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(xG),N);
traj_opt = traj_opt.addRunningCost(@cost);
traj_opt = traj_opt.addFinalCost(@finalCost);
traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xG)]));
      
      

% for i=1:10
% h = rand;
% x=randn(nx,1);
% u=randn(nu,1);
% l=rand(nc*2,1);
% 
% 
% tmp1 = @(h,x,u,l,a) traj_opt.lambda_constraint_fun(h,x,u,l,a); 
% 
% [f1,df1] = tmp1(h,x,u,l,0.1*l);
% [f2,df2] = geval(tmp1,h,x,u,l,0.1*l,struct('grad_method','numerical'));
% 
% valuecheck(df1,df2,1e-3);
% 
% end
% 
% keyboard

      
% add a display function to draw the trajectory on every iteration
function displayStateTrajectory(t,x,u)
  ts = [0,cumsum(t)'];
  xtraj = PPTrajectory(foh(ts,x));
  xtraj = xtraj.setOutputFrame(r.getStateFrame);
  v.playback(xtraj);
end
traj_opt = addTrajectoryDisplayFunction(traj_opt,@displayStateTrajectory);
            

tic;
[xtraj,utraj,z,F,info] = traj_opt.solveTraj(tf0,traj_init);
toc

v.playback(xtraj,struct('slider',true));

keyboard




function [g,dg] = cost(h,x,u)
  Q = 0*diag([100*ones(nx/2,1); 1*ones(nx/2,1)]);
  R = 0.01*eye(nu);
 
  g = (x-xG)'*Q*(x-xG) + u'*R*u;
  dg = [0, 2*(x'*Q -xG'*Q), 2*u'*R];
end

function [g,dg] = finalCost(t,x)
  Q = 0*diag([100*ones(nx/2,1); 1*ones(nx/2,1)]);

  g = (x-xG)'*Q*(x-xG);
  dg = [0,2*(x'*Q - xG'*Q)];
end



end