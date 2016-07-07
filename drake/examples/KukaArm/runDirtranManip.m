function runDirtranManip(N)

options.floating = false;
options.terrain = RigidBodyFlatTerrain();
w = warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits');
r = RigidBodyManipulator('urdf/iiwa14.urdf',options);
warning(w);

nx = r.getNumStates;
nq = r.getNumPositions;
nu = r.getNumInputs;

v=r.constructVisualizer;

x0 = double(r.resolveConstraints(zeros(nx,1)));
xG = x0;
xG(1) = xG(1)+0.1;
xG(4) = xG(4)+0.3;

tf0 = 1.0;

traj_init.x = PPTrajectory(foh([0,tf0],[x0,xG]));
traj_init.u = ConstantTrajectory(zeros(nu,1));
  
traj_opt = DirtranTrajectoryOptimization(r,N,tf0*[(1-0) (1+0)],options);
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0),1);
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(xG),N);
traj_opt = addTrajectoryDisplayFunction(traj_opt,@displayStateTrajectory);
% traj_opt = traj_opt.addRunningCost(@cost);
% traj_opt = traj_opt.addFinalCost(@finalCost);

[jlmin,jlmax] = r.getJointLimits;
xmin = [jlmin;-inf(nq,1)];
xmax = [jlmax;inf(nq,1)];
traj_opt = traj_opt.addConstraint(BoundingBoxConstraint(repmat(xmin,N,1),repmat(xmax,N,1)),traj_opt.x_inds);


tic;
[xtraj,utraj,z,F,info,infeasible] = traj_opt.solveTraj(tf0,traj_init);
toc

x = z(traj_opt.x_inds);
u = z(traj_opt.u_inds);

v.playback(xtraj,struct('slider','true'))

keyboard

  % add a display function to draw the trajectory on every iteration
  function displayStateTrajectory(t,x,u)
    ts = [0,cumsum(t)'];
    xtraj = PPTrajectory(foh(ts,x));
    xtraj = xtraj.setOutputFrame(r.getStateFrame);
    v.playback(xtraj);
  end

  function [g,dg] = cost(h,x,u)
    Q = diag(xcost);
    R = 0.0*eye(nu);

    g = (x-xG)'*Q*(x-xG) + u'*R*u;
    dg = [0, 2*(x'*Q -xG'*Q), 2*u'*R];
  end

  function [g,dg] = finalCost(T,x)
    Q = 10*eye(nx);

    g = (x-xG)'*Q*(x-xG);
    dg = [0, 2*(x'*Q -xG'*Q)];
  end


end