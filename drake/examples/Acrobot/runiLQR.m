function runiLQR()

warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits');
options.twoD = true;
options.view = 'right';
s = 'Acrobot.urdf';
dt = 0.01;
r = TimeSteppingRigidBodyManipulator(s,dt,options);

nx = r.getNumStates;
nu = r.getNumInputs;

v = r.constructVisualizer();

x_ind = 1:nx;
u_ind = nx+(1:nu);

% control limits
Op.lims = [r.umin, r.umax];
Op.parallel = false;

% optimization problem
DYNCST  = @(x,u,i) dyn_cost(x,u);
T = 0.1; % traj time
N = 1;%T/dt; % horizon
q0 = [0;0];
x0 = [q0;0*q0];
x0 = double(r.resolveConstraints(x0));
u0 = .1*randn(nu,N);    % initial controls

xG = [pi;0;0;0];


% test gradients

for j=1:10
  xr = randn(nx,1);
  ur = randn(nu,1);

  [~,df1] = geval(@cost,xr,ur,struct('grad_method','numerical'));
  [~,df2] = cost(xr,ur);

  valuecheck(df1,df2,1e-4);

  [~,df1] = geval(@final_cost,xr,struct('grad_method','numerical'));
  [~,df2] = final_cost(xr);

  valuecheck(df1,df2,1e-4);

  xr = randn(nx,N+1);
  ur = randn(nu,N+1);

  [f1,df1] = geval(@tmp1,xr,struct('grad_method','numerical'));
  [f2,df2] = tmp1(xr);
keyboard

  [f1,df1] = geval(@tmp2,xr,struct('grad_method','numerical'));
  [f2,df2] = tmp2(xr);
keyboard


end


% run the optimization
[xtraj, utraj, L, Vx, Vxx, total_cost, trace, stop]  = iLQG(DYNCST, x0, u0, Op);


ts = linspace(0,T,N+1);
xtraj = PPTrajectory(foh(ts,xtraj));
xtraj = xtraj.setOutputFrame(r.getStateFrame);
v.playback(xtraj,struct('slider',true));

function [g,dg,d2g] = cost(x,u)
  Q = diag([100*ones(nx/2,1); 0.001*ones(nx/2,1)]);
  R = 0.01*eye(nu);
  
  g = (x-xG)'*Q*(x-xG) + u'*R*u;
  if nargout > 1
    dg = [2*(x'*Q -xG'*Q), 2*u'*R];
    d2g = [2*Q, zeros(nx,nu); 
           zeros(nu,nx), 2*R];
  end
end

function [g,dg,d2g] = final_cost(x)
  Q = diag([100*ones(nx/2,1); 10*ones(nx/2,1)]);

  g = (x-xG)'*Q*(x-xG);
  if nargout > 1
    dg = 2*(x'*Q - xG'*Q);
    d2g = 2*Q;
  end
end

function [f,fx,fxx] = tmp1(x)
  f = [];
  for i=1:N+1
    [f_,~] = dyn_cost(x(:,i),ur(:,i));
    f = [f,f_];
  end
  [~,~,fx,~,fxx] = dyn_cost(x,ur);  
end

function [f,fu,fuu] = tmp2(u)
  f = [];
  for i=1:N+1
    [f_,~] = dyn_cost(xr(:,i),u(:,i));
    f = [f,f_];
  end
  [~,~,~,fu,~,fuu] = dyn_cost(xr,u);  
end

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