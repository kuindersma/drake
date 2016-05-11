function runiLQRWithBall(robust)

warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits');
options.view = 'right';
options.update_convex = true;

dt = 0.005;
r = TimeSteppingRigidBodyManipulator(PlanarRigidBodyManipulator('AcrobotCollision.urdf'),dt,options);
options.floating = true;
r = r.addRobotFromURDF('../../systems/plants/test/ball.urdf',zeros(3,1),zeros(3,1), options);

% create LCP copy
options.update_convex = false;
rlcp = TimeSteppingRigidBodyManipulator(PlanarRigidBodyManipulator('AcrobotCollision.urdf'),dt,options);
options.floating = true;
rlcp = rlcp.addRobotFromURDF('../../systems/plants/test/ball.urdf',zeros(3,1),zeros(3,1), options);

nx = r.getNumStates;
nu = r.getNumInputs;
nw = 4; % num_c * num_d

v = r.constructVisualizer();

x_ind = 1:nx;
u_ind = nx+(1:nu);
w_ind = nx+nu+(1:nw);

Op.parallel = false;
Op.plotFn = @plotfn;

  function plotfn(x)
    ts = linspace(0,T,N+1);
    xtraj = PPTrajectory(foh(ts,x));
    xtraj = xtraj.setOutputFrame(r.getStateFrame);
    v.playback(xtraj);
  end

% optimization problem
T = 1.0; % traj time
N = T/dt; % horizon

acrobot_pos = [0;0];
acrobot_vel = -[.5;5];
ball_pos = [1.5;.5;0];
ball_vel = [0;0;0];

x0 = [acrobot_pos; ball_pos; acrobot_vel; ball_vel];

acrobot_goal_pos = [pi;0];
acrobot_goal_vel = [0;0];
ball_goal_pos = [-2;2;0];
ball_goal_vel = [0;0;0];

xG = [acrobot_goal_pos; ball_goal_pos; acrobot_goal_vel; ball_goal_vel];

if robust
%   Op.lims = [[10*r.umin;-10*ones(nw,1)],[10*r.umax; 10*ones(nw,1)]];
  DYNCST  = @(x,u,i) dyn_cost(x,u(1:nu,:),u(nu+(1:nw),:));
  u0 = .01*randn(nu+nw,N);    % initial controls
else
  Op.lims = [10*r.umin, 10*r.umax];
  DYNCST  = @(x,u,i) dyn_cost(x,u);
  u0 = .01*randn(nu,N);    % initial controls
end

% run the optimization
[xtraj, utraj, L, Vx, Vxx, total_cost, trace, stop]  = iLQG(DYNCST, x0, u0, Op);

ts = linspace(0,T,N+1);
xtraj = PPTrajectory(foh(ts,xtraj));
xtraj = xtraj.setOutputFrame(r.getStateFrame);
v.playback(xtraj,struct('slider',true));

keyboard

ts = linspace(0,T,N);
utraj = PPTrajectory(foh(ts,utraj));
utraj = utraj.setOutputFrame(r.getInputFrame);
sys1 = cascade(utraj,r);

traj1 = sys1.simulate([0,T],x0);
v.playback(traj1,struct('slider',true));

keyboard

utraj = utraj.setOutputFrame(rlcp.getInputFrame);
sys2 = cascade(utraj,rlcp);
traj2 = sys2.simulate([0,T],x0);
traj2 = traj2.setOutputFrame(r.getStateFrame);
v.playback(traj2,struct('slider',true));



function [g,dg,d2g] = cost(x,u,w)
  Q = 0*diag([10*ones(nx/2,1); 0.001*ones(nx/2,1)]);
  R = 0.001*eye(nu);

  g = (x-xG)'*Q*(x-xG) + u'*R*u;
  
  if nargin > 2
    G = 50*ones(nw);
    g = g - w'*G*w;
  end

  if nargout > 1
    if nargin < 3
      dg = [2*(x'*Q -xG'*Q), 2*u'*R];
      d2g = [2*Q, zeros(nx,nu); 
             zeros(nu,nx), 2*R];      
    else
      dg = [2*(x'*Q -xG'*Q), 2*u'*R, 2*w'*G];
      d2g = [2*Q, zeros(nx,nu+nw); 
             zeros(nu,nx), 2*R, zeros(nu,nw);
             zeros(nw,nx+nu), 2*G];      
    end
  end
end

function [g,dg,d2g] = final_cost(x)
  Q = zeros(nx,1);
  Q(3) = 200;
  Q(4) = 200;
  Q = diag(Q);
  
  g = (x-xG)'*Q*(x-xG);
  if nargout > 1
    dg = 2*(x'*Q - xG'*Q);
    d2g = 2*Q;
  end
end

function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = dyn_cost(x,u,w)

  if nargout == 2
    if size(x,2) > 1
      error('Dynamics are not vectorized');
    end
    if isnan(u) % final cost
      f = [];
      c = final_cost(x);
    else
      if nargin < 3
        f = r.updateConvex(0,x,u);
        c = cost(x,u);
      else
        f = r.updateConvex(0,x,u,w);
        c = cost(x,u,w);
      end
    end
  else
    % x should be nx x N+1 where N is the trajectory length
    % u should be nu x N+1 where N is the trajectory length. the last element
    %    is nan
    
    if nargin < 3
      nuw = nu;
      uw_ind = u_ind;
    else        
      nuw = nu+nw;
      uw_ind = nx+(1:(nu+nw));
    end
    
    fx = zeros(nx,nx,N+1);
    fu = zeros(nx,nuw,N+1);
    cx = zeros(nx,N+1);
    cu = zeros(nuw,N+1);
    cxx = zeros(nx,nx,N+1);
    cxu = zeros(nx,nuw,N+1);
    cuu = zeros(nuw,nuw,N+1);
    
    for i=1:N+1
      xi = x(:,i);
      ui = u(:,i);
      
      if isnan(ui)
        [~,dg,d2g] = final_cost(xi);
        % cu is 0
      else
        if nargin < 3
          [~,dg,d2g] = cost(xi,ui);
          [~,df] = geval(@r.updateConvex,0,xi,ui,struct('grad_method','numerical'));
        else
          wi = w(:,i);
          [~,dg,d2g] = cost(xi,ui,wi);
          [~,df] = geval(@r.updateConvex,0,xi,ui,wi,struct('grad_method','numerical'));
        end
        cu(:,i) = dg(uw_ind);
        fx(:,:,i) = full(df(:,1+x_ind)); % +1 to skip time argument
        fu(:,:,i) = full(df(:,1+uw_ind));
        cxu(:,:,i) = d2g(x_ind,uw_ind);
        cuu(:,:,i) = d2g(uw_ind,uw_ind);
      end
      cx(:,i) = dg(x_ind);
      cxx(:,:,i) = d2g(x_ind,x_ind);
    end
    [f,c,fxx,fxu,fuu] = deal([]); % full DP not implemented 
  end
end

end