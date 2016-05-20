classdef PendulumPlant < SecondOrderSystem
% Defines the dynamics for the Pendulum.
  
  properties
    m = 1;   % kg
    l = .5;  % m
    b = 0.1; % kg m^2 /s
    lc = .5; % m
    I = .25; %m*l^2; % kg*m^2
    g = 9.81; % m/s^2
    
    xG;
    uG;
    w_max =0;

  end
  
  methods
    function obj = PendulumPlant(b)
      % Construct a new PendulumPlant
      obj = obj@SecondOrderSystem(1,1,true);

      if nargin>0 && ~isempty(b) % accept damping as optional input
        obj.b = b;
      end
      
      obj = setInputFrame(obj,PendulumInput);
      torque_limit = 3;
      obj = setInputLimits(obj,-torque_limit,torque_limit);
      
      obj = setStateFrame(obj,PendulumState);
      obj = setOutputFrame(obj,PendulumState);
      
      obj.xG = Point(getStateFrame(obj),[pi;0]);
      obj.uG = Point(getInputFrame(obj),0);
    end
    
    function qdd = sodynamics(obj,t,q,qd,u)
      % Implement the second-order dynamics
      qdd = (u - obj.m*obj.g*obj.lc*sin(q) - obj.b*qd)/obj.I;
    end
    
    function obj = addDisturbanceBound(obj,dist)
      obj.w_max = dist;
    end

    function n = getNumDisturbances(obj)
      n =1;
    end
    
    function [f,df] = dynamics_w(obj,t,x,u,w)
      if (nargout>1)
        [f,df] = dynamics(obj,t,x,u+w);
        nx = obj.getNumStates;
        nu = obj.getNumInputs;
        df = [df,df(:,1+nx+(1:nu))];
      else
        f = dynamics(obj,t,x,u+w);
      end
    end
    
    function [f,df,d2f,d3f]=dynamics(obj,t,x,u)
      w = rand*(2*obj.w_max) - obj.w_max;
      f=dynamics@SecondOrderSystem(obj,t,x,u+w);
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    function [utraj,xtraj]=robustSwingUpTrajectory(obj,N,M,disturbances)
      x0 = [0;0]; 
      xf = double(obj.xG);
      tf0 = 4;
      
      d = linspace(-disturbances,disturbances,M);
      options.integration_method = DirtranTrajectoryOptimization.MIDPOINT;
      prog = RobustDirtranTrajectoryOptimization(obj,N,M,[2 6],options);
      disp('constructor done');
      prog = prog.setDisturbances(d);
      prog = prog.addStateConstraint(ConstantConstraint(x0),1);
      prog = prog.addStateConstraint(ConstantConstraint(xf),N);
      prog = prog.addRunningCost(@cost);
      prog = prog.addFinalCost(@finalCost);
      prog = prog.addRobustConstraints(@robust_cost);
            
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
      
      
      function [g,dg] = robust_cost(dx,du,w)
        W = 0*eye(length(w));
        Qw = 500*eye(2);
        Rw = 0.1;
        g = dx'*Qw*dx + du'*Rw*du + w'*W*w;
        dg = [2*dx'*Qw 2*du'*Rw, 2*w'*W];
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
%      options.stabilize=true;
%      options.verify=false;
      options.xs = [0;0];
      options.Tslb = 2;
      options.Tsub = 6;
      options.degL1=4;
      c = LQRTree.buildLQRTree(obj,obj.xG,obj.uG,@()rand(2,1).*[2*pi;10]-[pi;5],Q,R,options);
    end

  end

end
