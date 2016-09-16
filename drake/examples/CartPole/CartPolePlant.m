classdef CartPolePlant < Manipulator

  properties
    mc = 1;   % mass of the cart in kg
    mp = .2;    % mass of the pole (point mass at the end) in kg
    l = 0.5;   % length of the pole in m
    g = 9.81;  % gravity m/s^2
    mu = 1; %Coulomb friction for cart
    friction_on = 0; %turn on Coulomb friction
    
    xG;
    uG;
  end

  methods
    function obj = CartPolePlant
      obj = obj@Manipulator(2,1);
      obj = setInputLimits(obj,-30,30);
      obj = setOutputFrame(obj,obj.getStateFrame);  % allow full-state feedback

      obj.xG = Point(obj.getStateFrame,[0;pi;0;0]);
      obj.uG = Point(obj.getInputFrame,0);
    end

    function [H,C,B] = manipulatorDynamics(obj,q,qd)
      mc=obj.mc;  mp=obj.mp;  l=obj.l;  g=obj.g;
      mu = obj.mu;
      s = sin(q(2)); c = cos(q(2));

      H = [mc+mp, mp*l*c; mp*l*c, mp*l^2];
      C = [0 -mp*qd(2)*l*s; 0 0];
      G = [0; mp*g*l*s];
      B = [1; 0];
      
      if obj.friction_on
          G(1) = sign(qd(1))*(mu*(mc+mp)*g);
      end

      C = C*qd + G;
    end
    
    function [f,df,d2f,d3f] = dynamics(obj,t,x,u)
      f = dynamics@Manipulator(obj,t,x,u);
      
      if (nargout>1)
        [df,d2f,d3f]= dynamicsGradients(obj,t,x,u,nargout-1);
      end
    end
    
    function [f,df,d2f] = dynamics_w(obj,t,x,u,w)
        q = x(1:2);  qd = x(3:4);
        
        mc=obj.mc;  mp=obj.mp;  l=obj.l;  g=obj.g;
        s = sin(q(2)); c = cos(q(2));
        
        H = [mc+mp, mp*l*c; mp*l*c, mp*l^2];
        C = [0 -mp*qd(2)*l*s; 0 0];
        G = [0; mp*g*l*s];
        Bu = [1; 0];
        Bw = Bu;
        
        qdd = -H\(C*qd + G - Bu*u - Bw*w);
        
        f = [qd; qdd];
        
        if (nargout>1)
            [df,d2f] = dynamicsGradients_w(obj,t,x,u,w,nargout-1);
        end
    end

    function nW = getNumDisturbances(obj)
        nW = 1;
    end
    
    function x0 = getInitialState(obj)
      x0 = randn(4,1);
    end

  end

  methods
    function [c,V]=balanceLQR(obj)
      Q = diag([1 50 1 50]);
      R = .2;

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

    function [c,V]=balanceHinf(obj)
      Bw = [zeros(2,2); eye(2)];  %uncertainty affects velocities only
      Q = diag([1 50 1 50]);
      R = .1;
      gamma = 20;
      if (nargout<2)
        c = tiHinf(obj,obj.xG,obj.uG,Q,R,Bw,gamma);
      else
        if any(~isinf([obj.umin;obj.umax]))
          error('currently, you must disable input limits to estimate the ROA');
        end
        [c,V]=tiHinf(obj,obj.xG,obj.uG,Q,R,Bw,gamma);
        pp = feedback(obj.taylorApprox(0,obj.xG,obj.uG,3),c);
        options.method='levelSet';
        V=regionOfAttraction(pp,V,options);
      end
    end

    function [utraj,xtraj]=swingUpTrajectory(obj)
      x0 = zeros(4,1); tf0 = 4; xf = double(obj.xG);
      N = 41; % This controls the number of time samples (knot points) used
              % for trajectory optimization. The more time samples you have
              % the more accurate the trajectory will be. But the
              % computation time will increase as you increase this. I've
              % increased the number of time samples to make it more
              % accurate.

      obj = setInputLimits(obj,-inf,inf);
      prog = DircolTrajectoryOptimization(obj,N,[2 6]);
      prog = prog.addStateConstraint(ConstantConstraint(x0),1);
      prog = prog.addStateConstraint(ConstantConstraint(xf),N);
      prog = prog.addRunningCost(@cost);
      prog = prog.addFinalCost(@finalCost);

      function [g,dg] = cost(dt,x,u)
        R = 1;
        g = sum((R*u).*u,1);
        dg = [zeros(1,1+size(x,1)),2*u'*R];
      end

      function [h,dh] = finalCost(t,x)
        h = t;
        dh = [1,zeros(1,size(x,1))];
      end

      traj_init.x = PPTrajectory(foh([0,tf0],[x0,xf]));
      for attempts=1:10
        tic
        [xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
        toc
        if info==1, break; end
      end
    end
    
    function [utraj,xtraj,z,prog]=swingUpDirtran(obj,N)
        x0 = zeros(4,1); tf0 = 5; xf = double(obj.xG);
        
        obj = setInputLimits(obj,-10,10);
        prog = RobustDirtranTrajectoryOptimization(obj,N,1,zeros(4),eye(4),1,eye(4),[tf0 tf0]);
        prog = prog.addStateConstraint(ConstantConstraint(x0),1);
        prog = prog.addStateConstraint(ConstantConstraint(xf),N);
        prog = prog.addRunningCost(@cost);
        %prog = prog.addFinalCost(@finalCost);
        
        prog = prog.setSolverOptions('snopt','majoroptimalitytolerance', 1e-4);
        prog = prog.setSolverOptions('snopt','majorfeaasibilitytolerance', 1e-4);
        prog = prog.setSolverOptions('snopt','minorfeaasibilitytolerance', 1e-4);
        
        function [g,dg] = cost(dt,x,u)
            Q = dt*1*eye(4);
            R = dt*.1;
            g = .5*(x-xf)'*Q*(x-xf) + .5*u'*R*u;
            dg = [g, (x-xf)'*Q, u'*R];
        end
        
        function [h,dh] = finalCost(t,x)
            Qf = 100*eye(4);
            h = .5*(x-xf)'*Qf*(x-xf);
            dh = [0, (x-xf)'*Qf];
        end
        
        traj_init.x = PPTrajectory(foh([0,tf0],[x0,xf]));
        tic
        [xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
        toc
    end

    function [utraj,xtraj,z,prog]=robustSwingUp(obj,N,D,E0,Q,R,Qf)
        x0 = zeros(4,1); tf0 = 5; xf = double(obj.xG);
        
        obj = setInputLimits(obj,-10,10);
        prog = RobustDirtranTrajectoryOptimization(obj,N,D,E0,Q,R,Qf,[tf0 tf0]);
        prog = prog.addStateConstraint(ConstantConstraint(x0),1);
        prog = prog.addStateConstraint(ConstantConstraint(xf),N);
        prog = prog.addRunningCost(@cost);
        %prog = prog.addFinalCost(@finalCost);
        
        prog = prog.setSolverOptions('snopt','majoroptimalitytolerance', 1e-5);
        prog = prog.setSolverOptions('snopt','majorfeaasibilitytolerance', 1e-5);
        prog = prog.setSolverOptions('snopt','minorfeaasibilitytolerance', 1e-5);
        
        prog = prog.addRobustCost(Q,1,Qf);
        prog = prog.addRobustInputConstraint();
        
        function [g,dg] = cost(dt,x,u)
            Q = dt*1*eye(4);
            R = dt*.1;
            g = .5*(x-xf)'*Q*(x-xf) + .5*u'*R*u;
            dg = [g, (x-xf)'*Q, u'*R];
        end
        
        function [h,dh] = finalCost(t,x)
            Qf = 100*eye(4);
            h = .5*(x-xf)'*Qf*(x-xf);
            dh = [0, (x-xf)'*Qf];
        end
        
        %traj_init.u = uguess;
        %traj_init.x = xguess;
        traj_init.x = PPTrajectory(foh([0,tf0],[x0,xf]));
        tic
        [xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
        toc
    end
    
    function [tv,Vtraj]=trajectorySwingUpAndBalance(obj)
      [ti,Vf] = balanceLQR(obj);

      [utraj,xtraj]=swingUpTrajectory(obj);

      Q=diag([10,10,1,1]); R=.1;
      [tv,Vtraj] = tvlqr(obj,xtraj,utraj,Q,R,Vf);
      
      if nargout>1
        psys = taylorApprox(feedback(obj,tv),xtraj,[],3);
        
        options.rho0_tau = 5; % The funnel is defined as the rho(t) sublevel set of the Lyapunov function V(t).
        % This parameter specifies the initial guess
        % for rho(t). In particular, the initial rho(t)
        % is specified as an exponential function with
        % options.rho0_tau being the exponent (after
        % time is scaled). See line 119 in
        % sampledFiniteTimeVerification.m for the exact
        % function. There's no magic formula for
        % choosing this. You have to play around with
        % it a bit (but in general, if a particular
        % value doesn't work, increasing it is a
        % reasonable choice).
        
        options.max_iterations = 3; % This controls the maximum number of iterations the bilinear alternations run for.
        % It's set relatively low right now.
        
        ts = xtraj.getBreaks();
        ts = linspace(ts(1),ts(end),12);  % These are time samples at which the Lyapunov conditions for the funnel are checked.
        % I've chosen a slightly coarse
        % sampling here to make things
        % faster. You can increase it if
        % you need more accuracy.
        Vtraj=sampledFiniteTimeVerification(psys,ts,Vf,Vtraj,options);
      end
    end
  end
end
