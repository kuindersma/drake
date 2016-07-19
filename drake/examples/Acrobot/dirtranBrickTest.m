function dirtranBrickTest(N)

options.use_bullet = false;
options.update_convex = true;
options.view = 'right';
options.terrain = RigidBodyFlatTerrain();
options.floating = true;

r = TimeSteppingRigidBodyManipulator(PlanarRigidBodyManipulator('brick_point_contact.urdf',options),.001,options);

nx = r.getNumStates;
nu = r.getNumInputs;

v = r.constructVisualizer();

brick_pos = [0;1.1;0];
brick_vel = [0;0;0];
x0 = [brick_pos; brick_vel];

brick_goal_pos = [0;0.9;0];
brick_goal_vel = [0;0;0];
xf = [brick_goal_pos; brick_goal_vel];

tf0 = 1.5;
% 
% traj = simulate(r,[0,tf0],x0);
% v.playback(traj,struct('slider',true));
% keyboard

options.linc_slack = 1;
traj_opt = NewContactImplicitTrajectoryOptimization(r,N,tf0*[(1-0.1) (1+0.1)],options);
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0),1);
% traj_opt = traj_opt.addStateConstraint(ConstantConstraint(xf),N);
% traj_opt = traj_opt.addRunningCost(@cost);
% traj_opt = traj_opt.addFinalCost(@final_cost);
traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));

% add a display function to draw the trajectory on every iteration
function displayStateTrajectory(hs,x,u)
  xtraj = PPTrajectory(foh(cumsum([0;hs]),x));
  xtraj = xtraj.setOutputFrame(r.getStateFrame);
  v.playback(xtraj);
end
  
% traj_opt = addTrajectoryDisplayFunction(traj_opt,@displayStateTrajectory);
     
% 
% for i=1:20
%   h = rand;
%   x=double(r.resolveConstraints(randn(nx,1)));
%   x2=double(r.resolveConstraints(randn(nx,1)));
%   u=randn(nu,1);
%   nc = r.getNumContactPairs();
%   nl = nc*2;
%   l=rand(nl,1);
%   a = rand(nl,1);
%   b = rand(nl,1);
% 
%   i
% 
%   tmp1 = @(h,x,u,l) traj_opt.forward_dynamics_fun(h,x,u,l); 
%   [f1,df1] = tmp1(h(1),x(:,1),u(:,1),l(:,1));
%   [f2,df2] = geval(tmp1,h(1),x(:,1),u(:,1),l(:,1),struct('grad_method','numerical'));
%   try
%     valuecheck(df1,df2,1e-3);
%   catch
%     keyboard
%   end
% 
%   tmp1 = @(h,x,x2,u,l) traj_opt.forward_constraint_fun(h,x,x2,u,l); 
%   [f1,df1] = tmp1(h,x,x2,u,l);
%   [f2,df2] = geval(tmp1,h,x,x2,u,l,struct('grad_method','numerical'));
%   try
%     valuecheck(df1,df2,1e-3);
%   catch
%     keyboard
%   end
%   
%   
%   tmp1 = @(x) traj_opt.phi_bound(x); 
%   [f1,df1] = tmp1(x);
%   [f2,df2] = geval(tmp1,x,struct('grad_method','numerical'));
%   try
%     valuecheck(df1,df2,1e-3);
%   catch
%     keyboard
%   end
%   
%     
%   tmp1 = @(x,l) traj_opt.phi_comp(x,l); 
%   [f1,df1] = tmp1(x,l);
%   [f2,df2] = geval(tmp1,x,l,struct('grad_method','numerical'));
%   try
%     valuecheck(df1,df2,1e-3);
%   catch
%     keyboard
%   end
%   
% 
%   tmp1 = @(h,x,u,l,a,b) traj_opt.lambda_constraint_fun(h,x,u,l,a,b); 
%   [f1,df1] = tmp1(h,x,u,l,a,b);
%   [f2,df2] = geval(tmp1,h,x,u,l,a,b,struct('grad_method','numerical'));
%   try
%     valuecheck(df1,df2,1e-3);
%   catch
%     keyboard
%   end
% 
% end
% 
% keyboard


tic
[xtraj,utraj,z,F,info,infeasible] = traj_opt.solveTraj(tf0,traj_init);
toc

h = z(traj_opt.h_inds);
x = z(traj_opt.x_inds);
u = z(traj_opt.u_inds);
l = z(traj_opt.l_inds);
a = z(traj_opt.alpha_inds);
b = z(traj_opt.beta_inds);

v.playback(xtraj,struct('slider',true));
keyboard

  function [g,dg] = cost(h,x,u)
    Q = zeros(nx);
    R = zeros(nu);

    xerr=(x-xf);
    g = xerr'*Q*xerr + u'*R*u;
    dg = [0,2*xerr'*Q, 2*u'*R];
  end

  function [g,dg] = final_cost(tf,x)
    Q = zeros(nx);
  
    xerr=(x-xf);
    g = xerr'*Q*xerr;
    dg = [0,2*xerr'*Q];
  end

end