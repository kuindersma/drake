classdef PendulumPlantDT < DrakeSystem
% Defines the dynamics for the Pendulum.
  
  properties
    p % PendulumPlant
    dt = 0.005;
    S = eye(2);
    xG = [pi;0];
  end
  
  methods
    function obj = PendulumPlantDT
      % Construct a new PendulumPlant
      obj = obj@DrakeSystem(0,2,1,2,false,true);
      obj = setSampleTime(obj,[obj.dt;0]); 
      
      obj.p = PendulumPlant;
      
      obj = setInputFrame(obj,PendulumInput);

      torque_limit = 3;
      obj.p = setInputLimits(obj.p,-torque_limit,torque_limit);
      obj = setInputLimits(obj,-torque_limit,torque_limit);
      
      obj = setStateFrame(obj,PendulumState);
      obj = setOutputFrame(obj,PendulumState);

      [~,V] = obj.runLQR();
      obj.S = V.S;
    end
    
    function [xdn,df,d2f]=update(obj,t,x,u)
      [f,df,d2f] = obj.p.dynamics(t,x,u);
      xdn=x+obj.dt*f;
      if nargout > 1
        nx = getNumStates(obj);
        nu = getNumInputs(obj);
        df = [zeros(nx,1),eye(nx),zeros(nx,nu)]+obj.dt*df;
        d2f = obj.dt*d2f;
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
        
  end
end
