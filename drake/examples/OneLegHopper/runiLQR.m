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
dt = 0.005;
r = TimeSteppingRigidBodyManipulator(s,dt,options);
r = r.setStateFrame(OneLegHopperState(r));
r = r.setOutputFrame(OneLegHopperState(r));

nq = r.getNumPositions;
nx = r.getNumStates;
nu = r.getNumInputs;
nc = r.getNumContactPairs;
nz = nc*2;

v = r.constructVisualizer();

x_ind = 1:nx;
uz_ind = nx+(1:nu+nz);

% control limits
Op.lims = [];%[r.umin, r.umax];
Op.parallel = false;

% optimization problem
DYNCST  = @(x,u,i) hopper_dyn_cost(x,u);
T = 1.0; % traj time
N = T/dt; % horizon
q0 = [0;0;.6;-1.2;.6+pi/2];
x0 = [q0;0*q0];
x0 = double(r.resolveConstraints(x0));
u0 = .1*randn(nu+nz,N);    % initial controls

xG =  zeros(nx,1); 
xG(1) = 0.0;
xG(2) = 1.0;

% run the optimization
[xtraj, u, L, Vx, Vxx, total_cost, trace, stop]  = iLQG(DYNCST, x0, u0, Op, energy_cost);

ddp = DifferentialDynamicProgramming(r,false);
xtraj = ddp.reconstructStateTrajectory(xtraj,0,dt);
v.playback(xtraj,struct('slider',true));


function [Q,c] = contact_cost_terms(x,u)

  R = 1e-8*eye(nq); % regularization for smoothing contacts
  q=x(1:nq);
  qd=x(nq+(1:nq));
  kinsol = doKinematics(obj,q,nargout>1);
  [H,C,B] = manipulatorDynamics(r,q,qd);
  [~,~,J] = contactConstraintsBV(obj,kinsol,obj.multiple_contacts);
  J = vertcat(J{:}); % 2km x n
  Hinv = inv(H);
  
  A = J*Hinv*J';
  c = J*qd + J*Hinv*(C+B*u)*dt;
  Q = A+R;
end

function [g,dg,d2g] = cost(x,u,z)
  Q = diag([100*ones(nx/2,1); 1.0*ones(nx/2,1)]);
  Ru = 0.01*eye(nu);
  Rz = 0.01*eye(nz);

  g = (x-xG)'*Q*(x-xG) + u'*Ru*u + z'*Rz*z;
  if nargout > 1
    dg = [2*(x'*Q -xG'*Q), 2*u'*Ru, 2*z'*Rz];
    d2g = [2*Q, zeros(nx,nu+nz); 
           zeros(nu,nx), 2*Ru, zeros(nu,nz); 
           zeros(nz,nx+nu),2*Rz];
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

function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = hopper_dyn_cost(x,uz,I)
  if nargout == 2
    if size(x,2) > 1
      error('Dynamics are not vectorized');
    end
    if isnan(uz) % final cost
      f = [];
      c = final_cost(x);
    else
      u = uz(1:nu);
      z = uz(nu+(1:nz));
      f = r.updateWithBasisForces(0,x,u,z);
      c = cost(x,u,z);
    end
  else
    % x should be nx x N+1 where N is the trajectory length
    % uz should be nu+nz x N+1 where N is the trajectory length. the last element
    %    is nan
    
    fx = zeros(nx,nx,N+1);
    fu = zeros(nx,nu+nz,N+1);
    cx = zeros(nx,N+1);
    cu = zeros(nu+nz,N+1);
    cxx = zeros(nx,nx,N+1);
    cxu = zeros(nx,nu+nz,N+1);
    cuu = zeros(nu+nz,nu+nz,N+1);
    
    for i=1:N+1
      xi = x(:,i);
      ui = uz(1:nu,i);
      zi = uz(nu+(1:nz),i);
      if isnan(ui)
        [~,dg,d2g] = final_cost(xi);
        % cu is 0
      else
        [~,dg,d2g] = cost(xi,ui,zi);
        cu(:,i) = dg(uz_ind);
        [~,df] = r.updateWithBasisForces(0,xi,ui,zi);
        fx(:,:,i) = full(df(:,1+x_ind)); % +1 to skip time argument
        fu(:,:,i) = full(df(:,1+uz_ind));
        cxu(:,:,i) = d2g(x_ind,uz_ind);
        cuu(:,:,i) = d2g(uz_ind,uz_ind);
      end
      cx(:,i) = dg(x_ind);
      cxx(:,:,i) = d2g(x_ind,x_ind);
    end
    [f,c,fxx,fxu,fuu] = deal([]); % full DP not implemented 
  end
end

end