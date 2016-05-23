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
    
    function [f,df] = dynamics_w(obj,~,x,u,w)

      % w is damping
      
      q=x(1:obj.num_q); 
      qd=x((obj.num_q+1):end);
      
      m_ = obj.m + w;
      l_ = obj.lc;
      b_ = obj.b;
      
      qdd = (u - m_*obj.g*l_*sin(q) - b_*qd)/obj.I;
      f = [qd;qdd];
      
      if nargout>1
        ii = 1.0/obj.I;
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
     
    function [utraj,xtraj]=swingUpTrajectory(obj,options)
      x0 = [0;0]; 
      xf = double(obj.xG);
      tf0 = 4;

      N = 21;
      traj_opt = DircolTrajectoryOptimization(obj,N,[2 6]);
      traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0),1);
      traj_opt = traj_opt.addStateConstraint(ConstantConstraint(xf),N);
      traj_opt = traj_opt.addRunningCost(@cost);
      traj_opt = traj_opt.addFinalCost(@finalCost);
      traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));
      
      
      function [g,dg] = cost(dt,x,u);
        R = 10;
        g = (R*u).*u;
        
        if (nargout>1)
          dg = [zeros(1,3),2*u'*R];
        end
      end
      
      function [h,dh] = finalCost(tf,x)
        h = tf;
        if (nargout>1)
          dh = [1, zeros(1,2)];
        end
      end
      
      info=0;
      while (info~=1)
        tic
        [xtraj,utraj,z,F,info] = traj_opt.solveTraj(tf0,traj_init);
        toc
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
    
    
    function [utraj,xtraj]=robustSwingUpTrajectory(obj,N,M,lb,ub)
      x0 = [0;0]; 
      xf = double(obj.xG);
      tf0 = 4;
      
      nw=1;
%       
%       if M>2
%         d = rand(nw,M).*repmat(ub-lb,1,M) + repmat(lb,1,M);
%       elseif M==1
%         d = lb;
%       end

%       d = [lb,ub];
        d = linspace(ub,lb,M);

%       d = zeros(nw,M);
%       d(:,1) = [lb(1);lb(2);lb(3)];
%       d(:,2) = [ub(1);lb(2);lb(3)];
%       d(:,3) = [lb(1);ub(2);lb(3)];
%       d(:,4) = [ub(1);ub(2);lb(3)];
%       d(:,5) = [lb(1);lb(2);ub(3)];
%       d(:,6) = [ub(1);lb(2);ub(3)];
%       d(:,7) = [lb(1);ub(2);ub(3)];
%       d(:,8) = [ub(1);ub(2);ub(3)];
      

      
      options.integration_method = DirtranTrajectoryOptimization.MIDPOINT;
      prog = RobustDirtranTrajectoryOptimization(obj,N,M,[2 10],options);
      disp('constructor done');
      prog = prog.setDisturbances(d);
      prog = prog.addStateConstraint(ConstantConstraint(x0),1);
      prog = prog.addStateConstraint(ConstantConstraint(xf),N);
      prog = prog.addRunningCost(@cost);
%       prog = prog.addFinalCost(@finalCost);
      prog = prog.addRobustConstraints(@robust_cost);
      
%       
% 
%    
%       nx=2; nu=1; nw=1;
%       for j=1:10
%         hr = randn();
%         xr = randn(nx,1);
%         xr2 = randn(nx,1);
%         ur = randn(nu,1);
%         ur2 = randn(nu,1);
%         wr = randn(nw,1);
%         gr = randn();
%       
%         [~,df1] = geval(@cost,hr,xr,ur,struct('grad_method','taylorvar'));
%         [~,df2] = cost(hr,xr,ur);
%       
%         valuecheck(df1,df2,1e-4);
%       
%         [~,df1] = geval(@robust_cost,xr,ur,wr,struct('grad_method','taylorvar'));
%         [~,df2] = robust_cost(xr,ur,wr);
%       
%         valuecheck(df1,df2,1e-4);
% 
%         
%         [~,df1] = geval(@obj.dynamics_w,hr,xr,ur,wr,struct('grad_method','taylorvar'));
%         [~,df2] = obj.dynamics_w(hr,xr,ur,wr);
%       
%         valuecheck(df1,df2,1e-4);
% 
%         
%         
%         [~,df1] = geval(@prog.forward_robust_dynamics_fun,hr,xr,ur,ur2,wr,struct('grad_method','taylorvar'));
%         [~,df2] = prog.forward_robust_dynamics_fun(hr,xr,ur,ur2,wr);
%       
%         valuecheck(df1,df2,1e-4);
%       
%         [~,df1] = geval(@prog.backward_robust_dynamics_fun,hr,xr,xr2,ur,ur2,wr,struct('grad_method','taylorvar'));
%         [~,df2] = prog.backward_robust_dynamics_fun(hr,xr,xr2,ur,ur2,wr);
%       
%         valuecheck(df1,df2,1e-4);
%       
%         [~,df1] = geval(@prog.midpoint_robust_dynamics_fun,hr,xr,xr2,ur,ur2,ur*.1,ur2*.1,wr,wr*.3,struct('grad_method','taylorvar'));
%         [~,df2] = prog.midpoint_robust_dynamics_fun(hr,xr,xr2,ur,ur2,ur*.1,ur2*.1,wr,wr*.3);
%       
%         valuecheck(df1,df2,1e-4);
%       
%         
%         tmp1 = @(gamma,h,x0,x1,u,du) prog.robust_forward_bound_fun(@robust_cost,1,gamma,h,x0,x1,u,du);
% 
%         [~,df1] = geval(tmp1,gr,hr,xr,xr2,ur,ur2,struct('grad_method','taylorvar'));
%         [~,df2] = tmp1(gr,hr,xr,xr2,ur,ur2);
%       
%         valuecheck(df1,df2,1e-4);
% 
%         tmp2 = @(gamma,h,x0,x1,u,du) prog.robust_backward_bound_fun(@robust_cost,1,gamma,h,x0,x1,u,du);
% 
%         [~,df1] = geval(tmp2,gr,hr,xr,xr2,ur,ur2,struct('grad_method','taylorvar'));
%         [~,df2] = tmp2(gr,hr,xr,xr2,ur,ur2);
%       
%         valuecheck(df1,df2,1e-4);
%         
%         tmp3 = @(gamma,h,x0,x1,u0,du0,u1) prog.robust_midpoint_bound_fun(@robust_cost,1,gamma,h,x0,x1,u0,du0,u1);
% 
%         [~,df1] = geval(tmp3,gr,hr,xr,xr2,ur,ur2,ur*.1,struct('grad_method','taylorvar'));
%         [~,df2] = tmp3(gr,hr,xr,xr2,ur,ur2,ur*.1);
%       
%         valuecheck(df1,df2,1e-4);
% 
% 
%       end
% 
%       keyboard
% 
%       
%       
%             
      % add a display function to draw the trajectory on every iteration
      function displayStateTrajectory(t,x,u)
        plot(x(1,:),x(2,:),'b.-','MarkerSize',10);
        drawnow;
      end
      prog = addTrajectoryDisplayFunction(prog,@displayStateTrajectory);
   
      
      traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));
      
      for attempts=1:1
        disp('Running solve');
        tic
        [xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
        toc
        if info==1, break; end
      end
      keyboard
    
      
      function [g,dg] = cost(dt,x,u);
        R = 1;
        g = (R*u).*u;
        
        if (nargout>1)
          dg = [zeros(1,3),2*u'*R];
        end
      end
      
      function [h,dh] = finalCost(tf,x)
        h = tf;
        if (nargout>1)
          dh = [1, zeros(1,2)];
        end
      end
      
      function [g,dg] = robust_cost(dx,du,w)
        W = 0*eye(length(w));
        Qw = 100*eye(2);
        Rw = 1;
        g = dx'*Qw*dx + du'*Rw*du + w'*W*w;
        dg = [2*dx'*Qw 2*du'*Rw, 2*w'*W];
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
