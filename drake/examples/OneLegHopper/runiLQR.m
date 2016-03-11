function runiLQR()

warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits');
options.twoD = true;
options.view = 'right';
options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.ignore_self_collisions = true;
options.use_bullet = false;
options.enable_fastqp = false;
s = 'OneLegHopper.urdf';
dt = 0.01;
r = TimeSteppingRigidBodyManipulator(s,dt,options);
r = r.setStateFrame(OneLegHopperState(r));
r = r.setOutputFrame(OneLegHopperState(r));

nx = r.getNumStates;
nu = r.getNumInputs;

v = r.constructVisualizer();

x_ind = 1:nx;
u_ind = nx+(1:nu);

% control limits
Op.lims = [r.umin, r.umax];
Op.parallel = false;

% optimization problem
DYNCST  = @(x,u,i) hopper_dyn_cost(x,u);
T = 0.5; % traj time
N = T/dt; % horizon
q0 = [0;0;.6;-1.2;.6+pi/2];
x0 = [q0;0*q0];
x0 = double(r.resolveConstraints(x0));
u0 = .1*randn(nu,N);    % initial controls


xG = x0;
xG(2) = xG(2) + 0.4;

% run the optimization
[xtraj, utraj, L, Vx, Vxx, total_cost, trace, stop]  = iLQG(DYNCST, x0, u0, Op);

ddp = DifferentialDynamicProgramming(r,false);
xtraj = ddp.reconstructStateTrajectory(xtraj,0,dt);
utraj = ddp.reconstructInputTrajectory(utraj,0,dt);
v.playback(xtraj,struct('slider',true));

save('hopper_iLQR_traj.mat','xtraj','utraj');

function [g,dg,d2g] = cost(x,u)
  Q = 0*diag([100*ones(nx/2,1); 1.0*ones(nx/2,1)]);
  R = 0.01*eye(nu);
  
  g = (x-xG)'*Q*(x-xG) + u'*R*u;
  if nargout > 1
    dg = [2*(x'*Q -xG'*Q), 2*u'*R];
    d2g = [2*Q, zeros(nx,nu); 
           zeros(nu,nx), 2*R];
  end
end

function [g,dg,d2g] = final_cost(x)
  Q = diag([100*ones(nx/2,1); 1.0*ones(nx/2,1)]);

  g = (x-xG)'*Q*(x-xG);
  if nargout > 1
    dg = 2*(x'*Q - xG'*Q);
    d2g = 2*Q;
  end
end

function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = hopper_dyn_cost(x,u,~)
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