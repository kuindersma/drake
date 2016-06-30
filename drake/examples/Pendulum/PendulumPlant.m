classdef PendulumPlant < SecondOrderSystem
% Defines the dynamics for the Pendulum.
  
  properties
    m = 1;   % kg
    l = .5;  % m
    b = 0.1; % kg m^2 /s
    lc = .5; % m
    I = .25; % m*l^2; % kg*m^2
    g = 9.81; % m/s^2
    
    xG;
    uG;
    u_pert_max = 0; % input disturbance bound
    disturb_inputs = false; 
  end
  
  methods
    function obj = PendulumPlant(b,options)

      if nargin < 2
        options = struct();
      end
      
      % Construct a new PendulumPlant
      obj = obj@SecondOrderSystem(1,1,true);

      if nargin>0 && ~isempty(b) % accept damping as optional input
        obj.b = b;
      end
      
      if isfield(options,'torque_limit')
        torque_limit = options.torque_limit;
      else
        torque_limit = 3;
      end

      if isfield(options,'mass')
        obj.m = options.mass;
      end

      if isfield(options,'length')
        obj.lc = options.length;
      end
      
      if isfield(options,'damping')
        obj.b = options.damping;
      end
      
      if isfield(options,'disturb_inputs')
        obj.disturb_inputs = options.disturb_inputs;
      else
        obj.disturb_inputs = false;
      end

      obj = setInputFrame(obj,PendulumInput);
      obj = setInputLimits(obj,-torque_limit,torque_limit);
      
      obj = setStateFrame(obj,PendulumState);
      obj = setOutputFrame(obj,PendulumState);
      
      obj.xG = Point(getStateFrame(obj),[pi;0]);
      obj.uG = Point(getInputFrame(obj),0);
    end
    
    function qdd = sodynamics(obj,~,q,qd,u)
      % Implement the second-order dynamics
      qdd = (u - obj.m*obj.g*obj.lc*sin(q) - obj.b*qd)/obj.I;
    end
    
    function obj = addInputDisturbanceBound(obj,dist)
      obj.u_pert_max = dist;
    end

    function obj = setMass(obj,mass)
      obj.m = mass;
    end
    
    function obj = setDamping(obj,b)
      obj.b = b;
    end
    
    function obj = setCOMLength(obj,length)
      obj.lc = length;
    end
    
    function obj = setDisturbInputs(obj,flag)
      obj.disturb_inputs = flag;
    end

    function n = getNumDisturbances(~)
      n = 1;
    end
    
%     function [f,df] = dynamics_w(obj,~,x,u,w)
% 
%       % w(1) is mass
%       % w(2) is length of dist to COM
%       % w(3) is damping
%       
%       q=x(1:obj.num_q); 
%       qd=x((obj.num_q+1):end);
%       
%       m_ = obj.m + w(1);
%       l_ = obj.lc + w(2);
%       b_ = obj.b + w(3);
%       
%       qdd = (u - m_*obj.g*l_*sin(q) - b_*qd)/obj.I;
%       f = [qd;qdd];
%       
%       if nargout>1
%         ii = 1.0/obj.I;
%         dfdt = zeros(2,1);
%         dfdx = [0, 1; -ii*m_*obj.g*l_*cos(q),-b_*ii];
%         dfdu = [0;ii];
%         dfdw = [zeros(1,3);-ii*obj.g*l_*sin(q), -ii*m_*obj.g*sin(q), -ii*qd];
%         df = [dfdt, dfdx, dfdu, dfdw];
%       end
%     end
    
    function [f,df] = dynamics_w(obj,x,u,w)

      % w is mass
      
      q=x(1:obj.num_q); 
      qd=x((obj.num_q+1):end);
      
      m_ = obj.m + w;
      l_ = obj.lc;
      b_ = obj.b;
      
      qdd = (u - m_*obj.g*l_*sin(q) - b_*qd)/(m_*l_*l_);
      f = [qd;qdd];
      
      if nargout>1
        ii = 1.0/(m_*l_*l_);
        dfdt = zeros(2,1);
        dfdx = [0, 1; -ii*m_*obj.g*l_*cos(q),-b_*ii];
        dfdu = [0;ii];
        dfdw = [0;-ii*obj.g*l_*sin(q)];
        df = [dfdt, dfdx, dfdu, dfdw];
      end
    end
    
    function [f,df,d2f,d3f]=dynamics(obj,t,x,u)
      if obj.disturb_inputs
        w = rand(obj.getNumInputs,1)*(2*obj.u_pert_max) - obj.u_pert_max;
      else
        w=zeros(obj.getNumInputs,1);
      end
      f = dynamics@SecondOrderSystem(obj,t,x,u+w);
      if (nargout>1)
        [df,d2f,d3f]= dynamicsGradients(obj,t,x,u,nargout-1);
      end
    end
    
    function [T,U] = energy(obj,x)
      theta = x(1);
      thetadot = x(2);
      T = .5*obj.m*obj.l^2*thetadot^2;
      U = -obj.m*obj.g*obj.l*cos(theta);
    end
    
    function x = getInitialState(obj)
      % Start me anywhere!
      x = randn(2,1);
    end
  end  
  
  methods 
    function [c,V]=balanceLQR(obj)
      Q = diag([10 1]); R = 1;
      if (nargout<2)
        c = tilqr(obj,obj.xG,obj.uG,Q,R);
      else
        if any(~isinf([obj.umin;obj.umax]))
          error('currently, you must disable input limits to estimate the ROA');
        end
        [c,V] = tilqr(obj,obj.xG,obj.uG,Q,R);
        pp = feedback(obj.taylorApprox(0,obj.xG,obj.uG,3),c);
        options.method='levelSet';
        V=regionOfAttraction(pp,V,options);
      end
    end  
     
function [utraj,xtraj,z,traj_opt]=swingUpTrajectory(obj,N,options)
      x0 = [0;0]; 
      xf = double(obj.xG);
      tf = 4;

      if nargin == 1
          N = 21;
          options = struct();
      elseif nargin == 2
          options = struct();
      end
      traj_opt = DirtranTrajectoryOptimization(obj,N,tf,options);
      traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0),1);
      traj_opt = traj_opt.addStateConstraint(ConstantConstraint(xf),N);
      traj_opt = traj_opt.addRunningCost(@cost);
      traj_opt = traj_opt.addFinalCost(@finalCost);
      traj_init.x = PPTrajectory(foh([0,tf],[double(x0),double(xf)]));
      
      
      function [g,dg] = cost(dt,x,u);
        xg = double(obj.xG);
        R = .1;
        Q = [.1 0; 0 .01];
        pmpi = [pi; -pi];
        [e1, i] = min(abs([x(1)+pi x(1)-pi]));
        e = [e1; x(2)-xg(2)];
        g = e'*Q*e + u'*R*u;
        if (nargout>1)
          dg = [0, 2*[x(1)+pmpi(i), x(2)]*Q, 2*u'*R];
        end
      end
      
      function [h,dh] = finalCost(tf,x)
        xg = double(obj.xG);
        Qf = 0*eye(2);
        pmpi = [pi; -pi];
        [e1, i] = min(abs([x(1)+pi x(1)-pi]));
        e = [e1; x(2)-xg(2)];
        h = e'*Qf*e;
        if (nargout>1)
          dh = [0, 2*[x(1)+pmpi(i), x(2)]*Qf];
        end
      end
      
      % add a display function to draw the trajectory on every iteration
      function displayStateTrajectory(t,x,u)
        subplot(2,1,1);
        plot(x(1,:),x(2,:),'b.-','MarkerSize',10);
        subplot(2,1,2);
        plot([0; cumsum(t(1:end-1))],u(1:end-1),'r.-','MarkerSize',10);
        drawnow;
      end
      traj_opt = addTrajectoryDisplayFunction(traj_opt,@displayStateTrajectory);
      
      info=0;
      while (info~=1)
        tic
        [xtraj,utraj,z,F,info] = traj_opt.solveTraj(tf,traj_init);
        toc
      end
    end
    
       
    function [g,dg] = integrationCost(obj,h,x,u)
      Q = 0.1;
      [xdot,dxdot,ddxdot] = obj.dynamics(0,x,u);
      nx = obj.getNumStates;
      nu = obj.getNumInputs;

      e = 0.5*h^2*dxdot(:,1+(1:nx))*xdot;
      g = e'*Q*e;
      if (nargout>1)
        de_dh = h*dxdot(:,1+(1:nx))*xdot;

        ddxdot_dx = [ddxdot(:,4+(2:3)),ddxdot(:,8+(2:3))];
        M = full(ddxdot_dx(:,1:2)*xdot(1) + ddxdot_dx(:,3:4)*xdot(2));
        de_dx = 0.5*h^2*(dxdot(:,1+(1:nx))*dxdot(:,1+(1:nx)) + M);

        de_du = 0.5*h^2*dxdot(:,1+(1:nx))*dxdot(:,1+nx+(1:nu));
        dg = [2*e'*Q*de_dh,2*e'*Q*de_dx,2*e'*Q*de_du];
      end
    end
    
    function c=trajectorySwingUpAndBalance(obj)
      [ti,Vf] = balanceLQR(obj);
      Vf = 5*Vf;  % artificially prune, since ROA is solved without input limits

      c = LQRTree(obj.xG,obj.uG,ti,Vf);
      [utraj,xtraj]=swingUpTrajectory(obj);  
      
      Q = diag([10 1]);  R=1;
      [tv,Vswingup] = tvlqr(obj,xtraj,utraj,Q,R,Vf);
      psys = taylorApprox(feedback(obj,tv),xtraj,[],3);
      options.degL1=2;
      Vswingup=sampledFiniteTimeVerification(psys,xtraj.getBreaks(),Vf,Vswingup,options);

      c = c.addTrajectory(xtraj,utraj,tv,Vswingup);
      
      c = setInputFrame(c,c.getInputFrame.constructFrameWithAnglesWrapped([1;0]));
    end
    
    function c=balanceLQRTree(obj)
      Q = diag([10 1]); R = 1;

      options.num_branches=5;
      options.stabilize=true;
      options.verify=false;
      options.xs = [0;0];
      options.Tslb = 2;
      options.Tsub = 6;
      options.degL1=4;
      c = LQRTree.buildLQRTree(obj,obj.xG,obj.uG,@()rand(2,1).*[2*pi;10]-[pi;5],Q,R,options);
    end
 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    function [utraj,xtraj,z,prog]=robustSwingUpTrajectory(obj,N,D)
      x0 = [0;0]; 
      xf = double(obj.xG);
      tf = 4;
      
      nw=1;
      
      options.integration_method = DirtranTrajectoryOptimization.FORWARD_EULER;
      prog = RobustDirtranTrajectoryOptimization(obj,N,D,tf,options);
      prog = prog.addStateConstraint(ConstantConstraint(x0),1);
      prog = prog.addStateConstraint(ConstantConstraint(xf),N);
      prog = prog.addRunningCost(@cost);
      prog = prog.addFinalCost(@finalCost);
      prog = prog.addRobustCost(@robust_cost);
      prog = prog.addRobustConstraints(@robust_cost);
      
      disp('constructor done');

%       nx=2; nu=1; nw=1;
%       for j=1:20
%         hr = rand();
%         xr = 2*randn(nx,1);
%         x1r = 2*randn(nx,1);
%         dx1r = 0.1*randn(nx,1);
%         dxr = 0.1*randn(nx,1);
%         ur = randn(nu,1);
%         dur = randn(nu,1);
%         u1r = randn(nu,1);
%         wr = randn(nw,1);
%         gr = randn();
%         
%         [~,df1] = geval(@cost,hr,xr,ur,struct('grad_method','numerical'));
%         [~,df2] = cost(hr,xr,ur);
%       
%         valuecheck(df1,df2,1e-3);
% 
%         [~,df1] = geval(@finalCost,hr,xr,struct('grad_method','numerical'));
%         [~,df2] = finalCost(hr,xr);
%       
%         valuecheck(df1,df2,1e-3);
% 
%         
%         [~,df1] = geval(@robust_cost,xr,ur,wr,struct('grad_method','numerical'));
%         [~,df2] = robust_cost(xr,ur,wr);
%       
%         valuecheck(df1,df2,1e-3);        
%         
%         [~,df1] = geval(@prog.forward_robust_dynamics_fun,hr,xr,ur,dur,wr,struct('grad_method','numerical'));
%         [~,df2] = prog.forward_robust_dynamics_fun(hr,xr,ur,dur,wr);
%       
%         valuecheck(df1,df2,1e-3);
% 
%         [~,df1] = geval(@prog.midpoint_robust_dynamics_fun,hr,xr,dxr,x1r,ur,dur,u1r,wr,struct('grad_method','numerical'));
%         [~,df2] = prog.midpoint_robust_dynamics_fun(hr,xr,dxr,x1r,ur,dur,u1r,wr);
%       
%         valuecheck(df1,df2,1e-3);
% 
%         tmp1 = @(gamma,h,x0,dx0,x1,u,du) prog.forward_robust_bound_fun(@robust_cost,1,gamma,h,x0,dx0,x1,u,du);
%         [~,df1] = geval(tmp1,gr,hr,xr,dxr,x1r,ur,dur,struct('grad_method','numerical'));
%         [~,df2] = tmp1(gr,hr,xr,dxr,x1r,ur,dur);
%       
%         valuecheck(df1,df2,1e-3);
% 
%         tmp1 = @(gamma,h,x0,dx0,x1,u,du,u1) prog.midpoint_robust_bound_fun(@robust_cost,1,gamma,h,x0,dx0,x1,u,du,u1);
%         [~,df1] = geval(tmp1,gr,hr,xr,dxr,x1r,ur,dur,u1r,struct('grad_method','numerical'));
%         [~,df2] = tmp1(gr,hr,xr,dxr,x1r,ur,dur,u1r);
%       
%         valuecheck(df1,df2,1e-3);
%         
%         
%         tmp1 = @(gamma,h,x0,dx0,x1,u,du,w0) prog.forward_w_bound_fun(@robust_cost,gamma,h,x0,dx0,x1,u,du,w0);
%         [~,df1] = geval(tmp1,gr,hr,xr,dxr,x1r,ur,dur,wr,struct('grad_method','numerical'));
%         [~,df2] = tmp1(gr,hr,xr,dxr,x1r,ur,dur,wr);
%       
%         valuecheck(df1,df2,1e-3);
% 
%         tmp1 = @(gamma,h,x0,dx0,x1,u0,du0,u1,w0) prog.midpoint_w_bound_fun(@robust_cost,gamma,h,x0,dx0,x1,u0,du0,u1,w0);
%         [~,df1] = geval(tmp1,gr,hr,xr,dxr,x1r,ur,dur,u1r,wr,struct('grad_method','numerical'));
%         [~,df2] = tmp1(gr,hr,xr,dxr,x1r,ur,dur,u1r,wr);
%       
%         valuecheck(df1,df2,1e-3);
%         
% 
%         tmp1 = @(h,x0,dx0,x1,dx1,u,du,w) prog.forward_delta_x_constraint(h,x0,dx0,x1,dx1,u,du,w);
%         [~,df1] = geval(tmp1,hr,xr,dxr,x1r,dx1r,ur,dur,wr,struct('grad_method','numerical'));
%         [~,df2] = tmp1(hr,xr,dxr,x1r,dx1r,ur,dur,wr);
%       
%         valuecheck(df1,df2,1e-3);
% 
%         tmp1 = @(h,x0,dx0,x1,dx1j,u0,du0,u1,w) prog.midpoint_delta_x_constraint(h,x0,dx0,x1,dx1j,u0,du0,u1,w);
%         [~,df1] = geval(tmp1,hr,xr,dxr,x1r,dx1r,ur,dur,u1r,wr,struct('grad_method','numerical'));
%         [~,df2] = tmp1(hr,xr,dxr,x1r,dx1r,ur,dur,u1r,wr);
%       
%         valuecheck(df1,df2,1e-3);        
%         
%         [~,df1] = geval(@prog.linear_delta_u_constraint,xr,ur,struct('grad_method','numerical'));
%         [~,df2] = prog.linear_delta_u_constraint(xr,ur);
%       
%         valuecheck(df1,df2,1e-3);
%         
%         
%         
%       end
% 
%       keyboard
             
      % add a display function to draw the trajectory on every iteration
      function displayStateTrajectory(t,x,u)
        subplot(2,1,1);
        plot(x(1,:),x(2,:),'b.-','MarkerSize',10);
        subplot(2,1,2);
        plot([0; cumsum(t(1:end-1))],u(1:end-1),'r.-','MarkerSize',10);
        drawnow;
      end
      prog = addTrajectoryDisplayFunction(prog,@displayStateTrajectory);
      
      traj_init.x = PPTrajectory(foh([0,tf],[double(x0),double(xf)]));
      
      for attempts=1:1
        disp('Running solve');
        tic
        [xtraj,utraj,z,F,info] = prog.solveTraj(tf,traj_init);
        toc
        if info==1, break; end
      end
      
%       h = z(prog.h_inds);
%       x = z(prog.x_inds);
%       dx = z(prog.dx_inds);
%       u = z(prog.u_inds);
%       du = z(prog.du_inds);
%       w = z(prog.w_inds);
      
      %keyboard
      
      function [g,dg] = cost(dt,x,u);
        xg = double(obj.xG);
        R = .1;
        Q = [.1 0; 0 .01];
        pmpi = [pi; -pi];
        [e1, i] = min(abs([x(1)+pi x(1)-pi]));
        e = [e1; x(2)-xg(2)];
        g = e'*Q*e + u'*R*u;
        if (nargout>1)
          dg = [0, 2*[x(1)+pmpi(i), x(2)]*Q, 2*u'*R];
        end
      end
      
      function [h,dh] = finalCost(tf,x)
        xg = double(obj.xG);
        Qf = 0*eye(2);
        pmpi = [pi; -pi];
        [e1, i] = min(abs([x(1)+pi x(1)-pi]));
        e = [e1; x(2)-xg(2)];
        h = e'*Qf*e;
        if (nargout>1)
          dh = [0, 2*[x(1)+pmpi(i), x(2)]*Qf];
        end
      end
      
      function [g,dg,ddg] = robust_cost(dx,du,w)
        W = 0*eye(length(w));
        Qw = [1 0; 0 1];
        Rw = .1;
        g = dx'*Qw*dx + du'*Rw*du + w'*W*w;
        dg = [2*dx'*Qw, 2*du'*Rw, 2*w'*W];
        ddg = blkdiag(2*Qw, 2*Rw, 2*W);
      end
       
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  end

end
