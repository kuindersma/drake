classdef PendulumPlantDT < DrakeSystem
% Defines the dynamics for the Pendulum.
  
  properties
    p % PendulumPlant
    dt;
    S = eye(2);
    xG = [pi;0];
  end
  
  methods
    function obj = PendulumPlantDT(h)
      % Construct a new PendulumPlant
      obj = obj@DrakeSystem(0,2,1,2,false,true);
      obj = setSampleTime(obj,[h;0]); 
      
      obj.p = PendulumPlant;
      
      obj = setInputFrame(obj,PendulumInput);

      obj.dt = h;
      torque_limit = 3;
      obj.p = setInputLimits(obj.p,-torque_limit,torque_limit);
      obj = setInputLimits(obj,-torque_limit,torque_limit);
      
      obj = setStateFrame(obj,PendulumState);
      obj = setOutputFrame(obj,PendulumState);
% 
%       [~,V] = obj.runLQR();
%       obj.S = V.S;
    end
    
    function obj = setMass(obj,mass)
      obj.p.m = mass;
    end
    
    function [xdn,dxdn,d2xdn]=update(obj,h,x,u)
      [f,df,d2f] = obj.p.dynamics(0,x,u);
      xdn=x+h*f;
      if nargout > 1
        nx = getNumStates(obj);
        nu = getNumInputs(obj);
        dxdn = [f,eye(nx)+h*df(:,1+(1:nx)),h*df(:,1+nx+(1:nu))];
        d2xdn = obj.dt*d2f;
      end
    end
   
    function y=output(obj,t,x,u)
      y=x;
    end
    
    function x = getInitialState(obj)
      % Start me anywhere!
      x = randn(2,1);
    end
    
    function [g,dg,d2g] = cost(obj,x,u)
      R = 0.1;
      Q = diag([10 0.1]);
      
      g = (x-obj.xG)'*Q*(x-obj.xG) + (R*u).*u;
        
      if (nargout>1)
        dg = [2*(x'*Q -obj.xG'*Q), 2*u'*R];
        d2g = [2*Q, zeros(2,1); zeros(1,2), 2*R];
      end
    end
      
    function [h,dh,d2h] = finalCost(obj,x)
      h = (x-obj.xG)'*obj.S*(x-obj.xG);
      if (nargout>1)
        dh = 2*(x'*obj.S -obj.xG'*obj.S);
        d2h = 2*obj.S;
      end
    end
    
    function [c,V] = runLQR(obj)
      Q = diag([50 2]); R = 0.1;
      [c,V] = tilqr(obj,obj.xG,0,Q,R);
    end
     
    
    function [utraj,xtraj]=swingUpTrajectory(obj,N,options)
      x0 = [0;0]; 
      xf = double(obj.xG);
      tf0 = 4;

      options.integration_method = DirtranTrajectoryOptimization.DT_SYSTEM;
      traj_opt = DirtranTrajectoryOptimization(obj,N,[3 8],options);
      traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0),1);
      traj_opt = traj_opt.addStateConstraint(ConstantConstraint(xf),N);
      traj_opt = traj_opt.addRunningCost(@cost);
      traj_opt = traj_opt.addFinalCost(@finalCost);
      traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));
      
      % add a display function to draw the trajectory on every iteration
      function displayStateTrajectory(t,x,u)
        plot(x(1,:),x(2,:),'b.-','MarkerSize',10);
        drawnow;
      end
      traj_opt = addTrajectoryDisplayFunction(traj_opt,@displayStateTrajectory);
      
      function [g,dg] = cost(h,x,u)
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
   
      
%       info=0;
%       while (info~=1)
      tic
      [xtraj,utraj,z,F,info] = traj_opt.solveTraj(tf0,traj_init);
      toc
%       end
    end
    
    
    
    
  end
end
