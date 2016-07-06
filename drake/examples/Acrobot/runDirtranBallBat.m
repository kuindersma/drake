function runDirtranBallBat(N)

options.use_bullet = true;
options.update_convex = true;
options.view = 'right';
% options.ignore_self_collisions = false;
options.terrain = RigidBodyFlatTerrain();

r = TimeSteppingRigidBodyManipulator(PlanarRigidBodyManipulator('AcrobotCollision.urdf'),.001,options);
options.floating = true;
r = r.addRobotFromURDF('../../systems/plants/test/ball.urdf',zeros(3,1),zeros(3,1), options);

nx = r.getNumStates;
nq = r.getNumPositions;
nu = r.getNumInputs;

v = r.constructVisualizer();

acrobot_pos = [0;0];
acrobot_vel = -[.5;5];
ball_pos = [1.5;.5;0];
ball_vel = [0;0;0];

x0 = [acrobot_pos; ball_pos; acrobot_vel; ball_vel];

acrobot_goal_pos = [-pi;0];
acrobot_goal_vel = [0;0];
ball_goal_pos = [-2;2;0];
ball_goal_vel = [0;0;0];

xf = [acrobot_goal_pos; ball_goal_pos; acrobot_goal_vel; ball_goal_vel];

tf0 = 2;

options.linc_slack = 1e-3;
options.scale_factor = 50;
traj_opt = SmoothContactImplicitTrajectoryOptimization(r,N,tf0*[(1-0.1) (1+0.1)],options);
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0),1);
% traj_opt = traj_opt.addStateConstraint(ConstantConstraint(xf),N);
traj_opt = traj_opt.addRunningCost(@cost);
traj_opt = traj_opt.addFinalCost(@final_cost);
traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));

% add a display function to draw the trajectory on every iteration
function displayStateTrajectory(hs,x,u)
  xtraj = PPTrajectory(foh(cumsum([0;hs]),x));
  xtraj = xtraj.setOutputFrame(r.getStateFrame);
  v.playback(xtraj);
end
  
traj_opt = addTrajectoryDisplayFunction(traj_opt,@displayStateTrajectory);
      
% for j=1:10
%   h = rand;
%   x = randn(nx,1);
%   u = randn(nu,1);
% 
%   [~,df1] = geval(@r.update,h,x,u,struct('grad_method','numerical'));
%   [~,df2] = r.update(h,x,u);
% 
%   valuecheck(df1,df2,1e-5);
% end
% 
% keyboard
     



tic
[xtraj,utraj,z,F,info] = traj_opt.solveTraj(tf0,traj_init);
toc


v.playback(xtraj,struct('slider',true));
keyboard

  function [g,dg] = cost(h,x,u)
    Q = zeros(nx,1);
%     Q(3) = 10;
%     Q(4) = 10;
    Q = diag(Q);
    R = 0.0001*eye(nu);

    xerr=(x-xf);
    g = xerr'*Q*xerr + u'*R*u;
    dg = [0,2*xerr'*Q, 2*u'*R];
  end

  function [g,dg] = final_cost(tf,x)
    Q = zeros(nx,1);
    Q(3) = 100;
    Q(4) = 100;
    Q = diag(Q);

    xerr=(x-xf);
    g = xerr'*Q*xerr;
    dg = [0,2*xerr'*Q];
  end

end