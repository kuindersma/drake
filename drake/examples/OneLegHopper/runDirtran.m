function runDirtran(N)

warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits');
options.twoD = true;
options.view = 'right';
options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.ignore_self_collisions = true;
options.use_bullet = false;
s = 'OneLegHopper.urdf';
dt = 0.001;
r = TimeSteppingRigidBodyManipulator(s,dt,options);
r = r.setStateFrame(OneLegHopperState(r));
r = r.setOutputFrame(OneLegHopperState(r));

nq = r.getNumPositions;
nx = r.getNumStates;
nu = r.getNumInputs;
nc = r.getNumContactPairs;
nl = nc*2;

v = r.constructVisualizer();
v.display_dt=0.01;

qd_init = zeros(nq,1);
distance = 0.15;

q0 = [0;0;.6;-1.2;.6+pi/2];
phi_f = r.contactConstraints(q0);
q0(2) = -phi_f(1);
x0 = [q0;qd_init];

q1 = [-distance/2;0;.6;-1.2;.2+pi/2];
phi_f = r.contactConstraints(q1);
q1(2) = -phi_f(1) + 0.15;
x1 = [q1;qd_init];

qf = [-distance;0;.6;-1.2;.6+pi/2];
phi_f = r.contactConstraints(qf);
qf(2) = -phi_f(1);
xG = [qf;qd_init];

tf0 = 0.6;
t_init = linspace(0,tf0,N);
xtraj = PPTrajectory(foh(t_init,[linspacevec(x0,x1,floor(N/2)),linspacevec(x1,xG,N-floor(N/2))]));
utraj = PPTrajectory(foh(t_init,zeros(nu,N)));
l = zeros(nl*N,1);
alpha = zeros(nl*N,1);

slack = 10;

[xtraj,utraj,z,F,info,infeasible,traj_opt] = doTrajopt(tf0,xtraj,utraj,l,alpha,slack);
v.playback(xtraj,struct('slider',true));

h = z(traj_opt.h_inds);
x = z(traj_opt.x_inds);
u = z(traj_opt.u_inds);
l = z(traj_opt.l_inds);
a = z(traj_opt.alpha_inds);
keyboard

% for i=1:6
%   slack = 0.1*slack;
%   [xtraj,utraj,z,F,info,infeasible,traj_opt] = doTrajopt(xtraj.tspan(2),xtraj,utraj,z(traj_opt.l_inds),z(traj_opt.alpha_inds),slack);
%   
%   h = z(traj_opt.h_inds);
%   x = z(traj_opt.x_inds);
%   u = z(traj_opt.u_inds);
%   l = z(traj_opt.l_inds);
%   a = z(traj_opt.alpha_inds);
%   keyboard
% 
% end

v.playback(xtraj,struct('slider',true));
keyboard




% tmp1 = @(h,x,u,l) traj_opt.forward_dynamics_fun(h,x,u,l); 
% for i=1:N
%   
%   [f1,df1] = tmp1(h(1),x(:,i),u(:,i),l(:,i));
%   [f2,df2] = geval(tmp1,h(1),x(:,i),u(:,i),l(:,i),struct('grad_method','numerical'));
%   
%   try
%     valuecheck(df1,df2,1e-3);
%   catch
%     keyboard
%   end
% end





% 
% keyboard
% 
% % open loop
% v.playback(simulate(cascade(utraj,r),xtraj.tspan,xtraj.eval(0)),struct('slider',true));

keyboard

  function [g,dg] = cost(h,x,u)
    Q = 0*diag([100*ones(nx/2,1); 1*ones(nx/2,1)]);
    R = 0.1*eye(nu);

    g = (x-xG)'*Q*(x-xG) + u'*R*u;
    dg = [0, 2*(x'*Q -xG'*Q), 2*u'*R];
  end

  function [g,dg] = finalCost(t,x)
    Q = diag([10*ones(nx/2,1); 10*ones(nx/2,1)]);
%     Q = 0*eye(nx);
  %   Q(2,2) = 100;

    g = (x-xG)'*Q*(x-xG);
    dg = [0,2*(x'*Q - xG'*Q)];
  end

  function [xtraj,utraj,z,F,info,infeasible,traj_opt] = doTrajopt(tf0,xtraj,utraj,l,alpha,slack)

    t_init = linspace(0,tf0,N);
    traj_init.x = xtraj;
    traj_init.u = utraj;
    traj_init.l = l;
    traj_init.alpha = alpha;
  
    options.linc_slack = slack;
    traj_opt = SmoothContactImplicitTrajectoryOptimization(r,N,tf0*[(1-0.2) (1+0.2)],options);
    traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0),1);
    traj_opt = traj_opt.addStateConstraint(ConstantConstraint(xG),N);
%     traj_opt = traj_opt.addRunningCost(@cost);
%     traj_opt = traj_opt.addFinalCost(@finalCost);
    traj_opt = addTrajectoryDisplayFunction(traj_opt,@displayStateTrajectory);

    tic;
    [xtraj,utraj,z,F,info,infeasible] = traj_opt.solveTraj(tf0,traj_init);
    toc
  
  end


  % add a display function to draw the trajectory on every iteration
  function displayStateTrajectory(t,x,u)
    ts = [0,cumsum(t)'];
    xtraj = PPTrajectory(foh(ts,x));
    xtraj = xtraj.setOutputFrame(r.getStateFrame);
    v.playback(xtraj);
  end

end