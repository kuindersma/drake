function checkDirtranGradients

N=10;

options.use_bullet = false;
options.update_convex = true;
options.view = 'right';
options.terrain = RigidBodyFlatTerrain();
options.floating = true;

r = TimeSteppingRigidBodyManipulator(PlanarRigidBodyManipulator('brick_point_contact.urdf',options),.001,options);

nx = r.getNumStates;
nu = r.getNumInputs;
nc = r.getNumContactPairs;

v = r.constructVisualizer();

tf0 = 1.0;

traj_opt = SmoothContactImplicitTrajectoryOptimization(r,N,tf0*[(1-0.2) (1+0.2)],options);

      
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




end