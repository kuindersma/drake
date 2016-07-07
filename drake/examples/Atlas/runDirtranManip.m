function runDirtranManip(N)

% silence some warnings
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints')
warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits')

options.floating = true;
options.use_bullet = false;
options.dt = 0.005;
r = Atlas('urdf/atlas_minimal_contact_with_box.urdf',options);
r = r.removeCollisionGroupsExcept({'heel','toe'});
r = compile(r);

nx = r.getNumStates;
nq = r.getNumPositions;
nu = r.getNumInputs;

v=r.constructVisualizer;

% set initial state to fixed point
load('data/xstar_no_arm.mat','xstar');
x0=xstar;
xG=x0;

q0 = x0(1:nq);
kinsol = r.doKinematics(q0);
rfoot_pos = forwardKin(r,kinsol,r.foot_body_id.right,[0;0;0],1);
lfoot_pos = forwardKin(r,kinsol,r.foot_body_id.left,[0;0;0],1);
lhand_pos = forwardKin(r,kinsol,findLinkId(r,'l_hand'),[0;0;0]);

lhand_pos(3) = lhand_pos(3)+0.05;

xcost = Point(r.getStateFrame,1);
xcost.base_z = 10;
xcost.base_roll = 100;
xcost.base_pitch = 100;
xcost.base_yaw = 100;
xcost.back_bky = 10;
xcost.back_bkx = 10;
xcost.back_bkz = 10;
xcost = double(xcost);
xcost(nq+(1:nq)) = 100;
tf0 =0.02*N;


traj_init.x = ConstantTrajectory(x0);
traj_init.u = ConstantTrajectory(zeros(nu,1));
  
traj_opt = DirtranExplicitForcesTrajectoryOptimization(r,N,tf0*[(1-0) (1+0)],rfoot_pos,lfoot_pos,lhand_pos,options);
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0),1);
% traj_opt = traj_opt.addStateConstraint(ConstantConstraint(xG),N);
% traj_opt = addTrajectoryDisplayFunction(traj_opt,@displayStateTrajectory);
traj_opt = traj_opt.addRunningCost(@cost);
traj_opt = traj_opt.addFinalCost(@finalCost);

tic;
[xtraj,utraj,z,F,info,infeasible] = traj_opt.solveTraj(tf0,traj_init);
toc

x = z(traj_opt.x_inds);
u = z(traj_opt.u_inds);
l = z(traj_opt.l_inds);

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
    Q = diag(xcost);

    g = (x-xG)'*Q*(x-xG);
    dg = [0, 2*(x'*Q -xG'*Q)];
  end


end