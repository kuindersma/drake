function runiLQRWithBall()

warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits');
% options.use_bullet = true;
options.view = 'right';
% options.ignore_self_collisions = false;
% options.terrain = RigidBodyFlatTerrain();

dt = 0.005;
r = TimeSteppingRigidBodyManipulator(PlanarRigidBodyManipulator('AcrobotCollision.urdf'),dt,options);
options.floating = true;
r = r.addRobotFromURDF('../../systems/plants/test/ball.urdf',zeros(3,1),zeros(3,1), options);

nx = r.getNumStates;
nu = r.getNumInputs;

v = r.constructVisualizer();

x_ind = 1:nx;
u_ind = nx+(1:nu);

% control limits
Op.lims = [10*r.umin, 10*r.umax];
Op.parallel = false;
Op.plotFn = @plotfn;

  function plotfn(x)
    ts = linspace(0,T,N+1);
    xtraj = PPTrajectory(foh(ts,x));
    xtraj = xtraj.setOutputFrame(r.getStateFrame);
    v.playback(xtraj);
  end

% optimization problem
DYNCST  = @(x,u,i) dyn_cost(x,u);
T = 1.0; % traj time
N = T/dt; % horizon

acrobot_pos = [0;0];
acrobot_vel = -[.5;5];
ball_pos = [1.5;.5;0];
ball_vel = [0;0;0];

x0 = [acrobot_pos; ball_pos; acrobot_vel; ball_vel];
% 
% q0 = [0;0];
% x0 = [q0;0*q0];
% x0 = double(r.resolveConstraints(x0));
u0 = .01*randn(nu,N);    % initial controls

acrobot_goal_pos = [pi;0];
acrobot_goal_vel = [0;0];
ball_goal_pos = [-2;2;0];
ball_goal_vel = [0;0;0];

xG = [acrobot_goal_pos; ball_goal_pos; acrobot_goal_vel; ball_goal_vel];


% test gradients
% 
% for j=1:10
%   xr = randn(nx,1);
%   ur = randn(nu,1);
% 
%   [~,df1,ddf1] = geval(@cost,xr,ur,struct('grad_method','taylorvar'));
%   [~,df2,ddf2] = cost(xr,ur);
% 
%   valuecheck(df1,df2,1e-4);
%   valuecheck(reshape(ddf1,nx+nu,nx+nu),ddf2,1e-4);
% 
%   [~,df1,ddf1] = geval(@final_cost,xr,struct('grad_method','taylorvar'));
%   [~,df2,ddf2] = final_cost(xr);
% 
%   valuecheck(df1,df2,1e-4);
%   valuecheck(reshape(ddf1,nx,nx),ddf2,1e-4);
% 
%   xr = randn(nx,N+1);
%   ur = randn(nu,N+1);
% 
%   [f1,df1] = geval(@tmp1,xr,struct('grad_method','numerical'));
%   [f2,df2] = tmp1(xr);
% keyboard
% 
%   [f1,df1] = geval(@tmp2,ur,struct('grad_method','numerical'));
%   [f2,df2] = tmp2(ur);
% keyboard
% 
% 
% end


% run the optimization
[xtraj, utraj, L, Vx, Vxx, total_cost, trace, stop]  = iLQG(DYNCST, x0, u0, Op);


ts = linspace(0,T,N+1);
xtraj = PPTrajectory(foh(ts,xtraj));
xtraj = xtraj.setOutputFrame(r.getStateFrame);
v.playback(xtraj,struct('slider',true));

function [g,dg,d2g] = cost(x,u)
  Q = 0*diag([10*ones(nx/2,1); 0.001*ones(nx/2,1)]);
  R = 0.001*eye(nu);
  
  g = (x-xG)'*Q*(x-xG) + u'*R*u;
  if nargout > 1
    dg = [2*(x'*Q -xG'*Q), 2*u'*R];
    d2g = [2*Q, zeros(nx,nu); 
           zeros(nu,nx), 2*R];
  end
end

function [g,dg,d2g] = final_cost(x)
  Q = zeros(nx,1);
  Q(3) = 200;
  Q(4) = 200;
  Q = diag(Q);
%   Q = 5*diag([100*ones(nx/2,1); 10*ones(nx/2,1)]);

  g = (x-xG)'*Q*(x-xG);
  if nargout > 1
    dg = 2*(x'*Q - xG'*Q);
    d2g = 2*Q;
  end
end

% function [f,fx,fxx] = tmp1(x)
%   f = [];
%   for i=1:N+1
%     [f_,~] = dyn_cost(x(:,i),ur(:,i));
%     f = [f,f_];
%   end
%   [~,~,fx,~,fxx] = dyn_cost(x,ur);  
% end
% 
% function [f,fu,fuu] = tmp2(u)
%   f = [];
%   for i=1:N+1
%     [f_,~] = dyn_cost(xr(:,i),u(:,i));
%     f = [f,f_];
%   end
%   [~,~,~,fu,~,fuu] = dyn_cost(xr,u);  
% end

function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = dyn_cost(x,u,~)
  if nargout == 2
    if size(x,2) > 1
      error('Dynamics are not vectorized');
    end
    if isnan(u) % final cost
      f = [];
      c = final_cost(x);
    else
      f = r.updateConvex(0,x,u);
      c = cost(x,u);
    end
  else
    % x should be nx x N+1 where N is the trajectory length
    % u should be nu x N+1 where N is the trajectory length. the last element
    %    is nan
    
    fx = zeros(nx,nx,N+1);
    fu = zeros(nx,nu,N+1);
    cx = zeros(nx,N+1);
    cu = zeros(nu,N+1);
    cxx = zeros(nx,nx,N+1);
    cxu = zeros(nx,nu,N+1);
    cuu = zeros(nu,nu,N+1);
    
    for i=1:N+1
      xi = x(:,i);
      ui = u(:,i);
      if isnan(ui)
        [~,dg,d2g] = final_cost(xi);
        % cu is 0
      else
        [~,dg,d2g] = cost(xi,ui);
        cu(:,i) = dg(u_ind);
        [~,df] = geval(@r.updateConvex,0,xi,ui,struct('grad_method','numerical'));
        fx(:,:,i) = full(df(:,1+x_ind)); % +1 to skip time argument
        fu(:,:,i) = full(df(:,1+u_ind));
        cxu(:,:,i) = d2g(x_ind,u_ind);
        cuu(:,:,i) = d2g(u_ind,u_ind);
      end
      cx(:,i) = dg(x_ind);
      cxx(:,:,i) = d2g(x_ind,x_ind);
    end
    [f,c,fxx,fxu,fuu] = deal([]); % full DP not implemented 
  end
end

end