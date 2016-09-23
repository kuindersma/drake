classdef CartDoublePendulumPlant < Manipulator

  properties
    mc = 1;     % mass of the cart in kg
    mp1 = .2;   % mass of the first pendulum (point mass at the end) in kg
    mp2 = .2;   % mass of the second pendulum (point mass at the end) in kg
    l1 = 0.5;    % length of the first pole in m
    l2 = 0.5;   % length of the second pole in m
    g = 9.81;  % gravity m/s^2
    
    xG;
    uG;
  end

  methods
    function obj = CartDoublePendulumPlant
      obj = obj@Manipulator(3,1);
      obj = setInputLimits(obj,-50,50);
      obj = setOutputFrame(obj,obj.getStateFrame);  % allow full-state feedback

      obj.xG = Point(obj.getStateFrame,[0;pi;pi;0;0;0]);
      obj.uG = Point(obj.getInputFrame,0);
    end

    function [H,C,B] = manipulatorDynamics(obj,q,qd)
      mc=obj.mc;  mp1=obj.mp1;  mp2=obj.mp2;
      l1=obj.l1;  l2=obj.l2;    g=obj.g;

      t1 = q(2)-pi;
      t2 = q(3)-pi;
      t1d = qd(2);
      t2d = qd(3);
      
      s1 = sin(t1);
      c1 = cos(t1);
      s2 = sin(t2);
      c2 = cos(t2);
      s12 = sin(t1-t2);
      c12 = cos(t1-t2);

      H = [mc+mp1+mp2,      l1*(mp1+mp2)*c1,    mp2*l2*c2;
           l1*(mp1+mp2)*c1, l1*l1*(mp1+mp2),    l1*l2*c12;
           l2*mp2*c2,       l1*l2*mp2*c12,      l2*l2*mp2];
      
      C = -[l1*(mp1+mp2)*t1d*t1d*s1 + mp2*l2*t2d*t2d*s2;
           -l1*l2*mp2*t2*t2*s12 + g*(mp1+mp2)*l1*s1;
            l1*l2*mp2*t1d*t1d*s12 + g*l2*mp2*s2];

      B = [1; 0; 0];
      
    end
    
    function [f,df,d2f] = dynamics(obj,t,x,u)
      f = dynamics@Manipulator(obj,t,x,u);
      
      if (nargout>1)
        [df,d2f]= dynamicsGradients(obj,t,x,u,nargout-1);
      end
    end
    
    function x0 = getInitialState(obj)
      x0 = randn(4,1);
    end

  end

  methods

    function [utraj,xtraj]=swingUpTrajectory(obj)
      x0 = zeros(6,1); tf0 = 8; xf = double(obj.xG);
      N = 41; % This controls the number of time samples (knot points) used
              % for trajectory optimization. The more time samples you have
              % the more accurate the trajectory will be. But the
              % computation time will increase as you increase this. I've
              % increased the number of time samples to make it more
              % accurate.

      obj = setInputLimits(obj,-inf,inf);
      prog = DircolTrajectoryOptimization(obj,N,[6 10]);
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
        
        prog = prog.setSolverOptions('snopt','majoroptimalitytolerance', 1e-3);
        prog = prog.setSolverOptions('snopt','majorfeaasibilitytolerance', 1e-4);
        prog = prog.setSolverOptions('snopt','minorfeaasibilitytolerance', 1e-4);
        
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
  end
end
