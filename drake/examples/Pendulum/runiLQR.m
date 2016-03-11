function runiLQR()

p = PendulumPlantDT;
pv = PendulumVisualizer();

nx=2;
nu=1;

% control limits
Op.lims = [p.umin, p.umax];
Op.parallel = false;

% optimization problem
DYNCST  = @(x,u,i) pend_dyn_cost(x,u);
T = 3.0; % traj time
N = T/p.dt; % horizon
x0 = zeros(nx,1);       % initial state
u0 = .1*randn(nu,N);    % initial controls

% run the optimization
[xtraj, u, L, Vx, Vxx, cost, trace, stop]  = iLQG(DYNCST, x0, u0, Op);

ts = linspace(0,T,N+1);
xtraj = PPTrajectory(foh(ts,xtraj));
xtraj = xtraj.setOutputFrame(p.getStateFrame);

pv.playback(xtraj,struct('slider',true));

function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = pend_dyn_cost(x,u,I)
  if nargout == 2
    if size(x,2) > 1
      error('Dynamics are not vectorized');
    end
    if isnan(u) % final cost
      f = [];
      c = p.finalCost(x);
    else
      f = p.update(0,x,u);
      c = p.cost(x,u);
    end
  else
    % x should be nxN+1 where N is the trajectory length
    % u should be mxN+1 where N is the trajectory length. the last element
    %    is nan
    x_ind = 1:nx;
    u_ind = nx+(1:nu);
    
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
        [~,dg,d2g] = p.finalCost(xi);
        % cu is 0
      else
        [~,dg,d2g] = p.cost(xi,ui);
        cu(:,i) = dg(u_ind);
        [~,df] = p.update(0,xi,ui);
        fx(:,:,i) = full(df(:,1+x_ind)); % +1 to skip time argument
        try
        fu(:,:,i) = full(df(:,1+u_ind));
        catch
          keyboard
        end
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