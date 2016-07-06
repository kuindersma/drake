function checkDirtranGradients

N=10;
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits');
options.twoD = true;
options.view = 'right';
options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.use_bullet = false;
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
      
    

for i=1:20
  h = rand;
  x=double(r.resolveConstraints(randn(nx,1)));
  x2=double(r.resolveConstraints(randn(nx,1)));
  u=randn(nu,1);
  nl = nc*2;
  l=rand(nl,1);
  a = rand(nl,1);
  b = rand(nc,1);

  i

  tmp1 = @(h,x,u,l) traj_opt.forward_dynamics_fun(h,x,u,l); 
  [f1,df1] = tmp1(h(1),x(:,1),u(:,1),l(:,1));
  [f2,df2] = geval(tmp1,h(1),x(:,1),u(:,1),l(:,1),struct('grad_method','numerical'));
  try
    valuecheck(df1,df2,1e-3);
  catch
    keyboard
  end



  tmp1 = @(h,x,x2,u,l) traj_opt.forward_constraint_fun(h,x,x2,u,l); 
  [f1,df1] = tmp1(h,x,x2,u,l);
  [f2,df2] = geval(tmp1,h,x,x2,u,l,struct('grad_method','numerical'));
  try
    valuecheck(df1,df2,1e-3);
  catch
    keyboard
  end
  
  
%   tmp1 = @(x) traj_opt.phi_bound(x); 
%   [f1,df1] = tmp1(x);
%   [f2,df2] = geval(tmp1,x,struct('grad_method','numerical'));
%   try
%     valuecheck(df1,df2,1e-3);
%   catch
%     keyboard
%   end

  tmp1 = @(h,x,u,l,a,b) traj_opt.lambda_constraint_fun(h,x,u,l,a,b); 
  [f1,df1] = tmp1(h,x,u,l,a,b);
  [f2,df2] = geval(tmp1,h,x,u,l,a,b,struct('grad_method','numerical'));
  try
    valuecheck(df1,df2,1e-3);
  catch
    keyboard
  end

  
  tmp1 = @(h,x,u,l,b) traj_opt.vmin_lagrange_eq(h,x,u,l,b); 
  [f1,df1] = tmp1(h,x,u,l,b);
  [f2,df2] = geval(tmp1,h,x,u,l,b,struct('grad_method','numerical'));
  try
    valuecheck(df1,df2,1e-3);
  catch
    keyboard
  end
  
  tmp1 = @(h,x,u,l) traj_opt.vmin_inequality(h,x,u,l); 
  [f1,df1] = tmp1(h,x,u,l);
  [f2,df2] = geval(tmp1,h,x,u,l,struct('grad_method','numerical'));
  try
    valuecheck(df1,df2,1e-3);
  catch
    keyboard
  end
  
end



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