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

dt = 0.0025;
m = PlanarRigidBodyManipulator('KneedCompassGait.urdf', options);
r = TimeSteppingRigidBodyManipulator(m,dt);

nq = r.getNumPositions;
nx = r.getNumStates;
nu = r.getNumInputs;

v = r.constructVisualizer();

x_ind = 1:nx;
u_ind = nx+(1:nu);

% control limits
Op.lims = [r.umin, r.umax];
Op.parallel = false;
Op.regType = 1; % regularization type 1: q_uu+lambda*eye(); 2: V_xx+lambda*eye()
Op.lambda = 1;
Op.plotFn = @plotfn;

  function plotfn(x)
    ts = linspace(0,T,N+1);
    xtraj = PPTrajectory(foh(ts,x));
    xtraj = xtraj.setOutputFrame(r.getStateFrame);
    v.playback(xtraj);
  end


% optimization problem
DYNCST  = @(x,u,i) dyn_cost(x,u,i);
T = 1.0; % traj time
N = T/dt; % horizon


% q0 = [0; 1; 0; 0;0;0;];
q0 = [0.009;0.997;-0.066;0.161;0.274;0.237];
v0 = 0*q0;

qG = [0.3+0.009; 0.997; 2*0.066;0.161; -0.274;0.237];
vG = 0*qG;

x0 = [q0;v0];
xG = [qG;vG];

xtraj_kin = PPTrajectory(foh([0,T],[x0,xG]));

% solve for fixed point
kinsol = doKinematics(r, q0);
[phi,~,~,~,~,~,~,~,n] = contactConstraints(r,kinsol,false);
[~,C,B] = manipulatorDynamics(r,q0,v0);

nc = length(phi<5e-3);
nparams = nc+nu;

gurobi_options.outputflag = 0; % verbose flag
gurobi_options.method = 1; % -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier

Iu = zeros(nu,nparams); Iu(:,1:nu) = eye(nu);
If = zeros(nc,nparams); If(:,nu+(1:nc)) = eye(nc);

model.Q=sparse(zeros(nparams));
model.obj = zeros(nparams,1);
model.A = sparse(B*Iu + n'*If);
model.rhs = C;
model.sense = repmat('=',length(C),1);
model.lb = [r.umin;zeros(nc,1)];
model.ub = [r.umax;inf(nc,1)];
result = gurobi(model,gurobi_options);
u0 = Iu*result.x;
u0 = repmat(u0,1,N) + 0.0*randn(nu,N);    % initial controls

% u0 = 0.1*randn(nu,N);    % initial controls

qG = q0;
qG(1) = q0(1)+0.5;
% qG(3) = -2*q0(3);
% qG(5) = -q0(5);
vG = v0;
xG = [qG;vG];

% run the optimization
[xtraj, utraj, L, Vx, Vxx, total_cost, trace, stop]  = iLQG(DYNCST, x0, u0, Op);

ts = linspace(0,T,N+1);
xtraj = PPTrajectory(foh(ts,xtraj));
xtraj = xtraj.setOutputFrame(r.getStateFrame);
v.playback(xtraj,struct('slider',true));

function [g,dg,d2g] = cost(x,u,i)
  Q = diag([10*ones(nx/2,1); 0*ones(nx/2,1)]);
  R = 0.01*eye(nu);

  xt = xtraj_kin.eval(i./N * T);

  g = (x-xt)'*Q*(x-xt) + u'*R*u;
  if nargout > 1
    dg = [2*(x'*Q -xt'*Q), 2*u'*R];
    d2g = [2*Q, zeros(nx,nu); 
           zeros(nu,nx), 2*R];
  end
end

function [g,dg,d2g] = final_cost(x)
  Q = diag([100*ones(nx/2,1); 0.0001*ones(nx/2,1)]);

  g = (x-xG)'*Q*(x-xG);
  if nargout > 1
    dg = 2*(x'*Q - xG'*Q);
    d2g = 2*Q;
  end
end

function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = dyn_cost(x,u,n)
  if nargout == 2
    if size(x,2) > 1
      error('Dynamics are not vectorized');
    end
    if isnan(u) % final cost
      f = [];
      c = final_cost(x);
    else
      f = r.updateConvex(0,x,u);
      c = cost(x,u,n);
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
        [~,dg,d2g] = cost(xi,ui,i);
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