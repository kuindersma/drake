classdef DirtranTrajectoryOptimization < DirectTrajectoryOptimization
  % Direct transcription trajectory optimization
  %  implements multiple possible integration schemes for the dynamics
  %  constraints xdot = f(x,u) and for for integrating the running cost
  %  term.
  %
  %  For forward euler integratino:
  %    dynamics constraints are: x(k+1) = x(k) + h(k)*f(x(k),u(k))
  %    integrated cost is sum of g(h(k),x(k),u(k))
  %  For backward euler integration:
  %    dynamics constraints are: x(k+1) = x(k) + h(k)*f(x(k+1),u(k))
  %    integrated cost is sum of g(h(k),x(k+1),u(k))
  %  For midpoint integration:
  %    dynamics constraints are: x(k+1) = x(k) + h(k)*f(.5*x(k)+.5*x(k+1),.5*u(k)+.5*u(k+1))
  %    integrated cost is sum of g(h(k),.5*x(k)+.5*x(k+1),.5*u(k)+.5*u(k+1))
  properties (Constant)
    FORWARD_EULER = 1;
    BACKWARD_EULER = 2;
    MIDPOINT = 3;  % DEFAULT
    DT_SYSTEM = 4;
  end
  
  methods
    function obj = DirtranTrajectoryOptimization(plant,N,duration,options)
      if nargin < 4
        options = struct();
      end
      if ~isfield(options,'integration_method')
        options.integration_method = DirtranTrajectoryOptimization.FORWARD_EULER;
      end
      
      obj = obj@DirectTrajectoryOptimization(plant,N,duration,options);
    end
    
    function obj = addDynamicConstraints(obj)
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
      N = obj.N;
      
      constraints = cell(N-1,1);
      dyn_inds = cell(N-1,1);
      
      switch obj.options.integration_method
        case DirtranTrajectoryOptimization.FORWARD_EULER
          n_vars = 2*nX + nU + 1;
          cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.forward_constraint_fun);
        case DirtranTrajectoryOptimization.BACKWARD_EULER
          n_vars = 2*nX + nU + 1;
          cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.backward_constraint_fun);
        case DirtranTrajectoryOptimization.MIDPOINT
          n_vars = 2*nX + 2*nU + 1;
          cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.midpoint_constraint_fun);
        case DirtranTrajectoryOptimization.DT_SYSTEM
          n_vars = 2*nX + nU + 1;
          cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.discrete_time_constraint_fun);
        otherwise
          error('Drake:DirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
      end
      
      for i=1:obj.N-1,
        switch obj.options.integration_method
          case DirtranTrajectoryOptimization.FORWARD_EULER
            dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i)};
          case DirtranTrajectoryOptimization.BACKWARD_EULER
            dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i)};
          case DirtranTrajectoryOptimization.MIDPOINT
            dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.u_inds(:,i+1)};
          case DirtranTrajectoryOptimization.DT_SYSTEM
            dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i)};
          otherwise
            error('Drake:DirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
        end
        cnstr = cnstr.setName(sprintf('dynamics_constr_%d_',i));
        constraints{i} = cnstr;
        
        obj = obj.addConstraint(constraints{i}, dyn_inds{i});
        
        if 0
        nc = obj.plant.getNumContactPairs;
        inds = {obj.x_inds(:,i)};
        constraint = FunctionHandleConstraint(zeros(nc,1),inf(nc,1),nX,@obj.phi_bound);
        constraint = constraint.setName(sprintf('phi_bound_%d',i));
        obj = obj.addConstraint(constraint, inds);
            end
      end
    end
    
    function obj = addRunningCost(obj,running_cost_function)
      % Adds an integrated cost to all time steps, which is
      % numerical implementation specific (thus abstract)
      % this cost is assumed to be time-invariant
      % @param running_cost_function a function handle
      %  of the form running_cost_function(dt,x,u)
      
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
      
      for i=1:obj.N-1,
        switch obj.options.integration_method
          case DirtranTrajectoryOptimization.FORWARD_EULER
            running_cost = FunctionHandleObjective(1+nX+nU, running_cost_function);
            inds_i = {obj.h_inds(i);obj.x_inds(:,i);obj.u_inds(:,i)};
          case DirtranTrajectoryOptimization.BACKWARD_EULER
            running_cost = FunctionHandleObjective(1+nX+nU, running_cost_function);
            inds_i = {obj.h_inds(i);obj.x_inds(:,i+1);obj.u_inds(:,i)};
          case DirtranTrajectoryOptimization.MIDPOINT
            running_cost = FunctionHandleObjective(1+2*nX+2*nU,...
              @(h,x0,x1,u0,u1) obj.midpoint_running_fun(running_cost_function,h,x0,x1,u0,u1));
            inds_i = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.u_inds(:,i+1)};
          case DirtranTrajectoryOptimization.DT_SYSTEM
            running_cost = FunctionHandleObjective(1+nX+nU, running_cost_function);
            inds_i = {obj.h_inds(i);obj.x_inds(:,i);obj.u_inds(:,i)};
          otherwise
            error('Drake:DirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
        end
        
        obj = obj.addCost(running_cost,inds_i);
      end
    end
    
    function [f,df] = phi_bound(obj,x)
      nq = obj.plant.getNumPositions;
      q = x(1:nq);
      kinsol = doKinematics(obj.plant, q);
      [f,~,~,~,~,~,~,~,dfdq] = obj.plant.contactConstraints(kinsol,obj.plant.multiple_contacts);
      df = [dfdq,0*dfdq];
    end

    function [f,df] = discrete_time_constraint_fun(obj,h,x0,x1,u)
      nX = obj.plant.getNumStates();
      [x1_,dx1_] = obj.plant.updateConvexOpt(h,x0,u);
      f = x1 - x1_;
      df = [-dx1_(:,1), -dx1_(:,2:1+nX), eye(nX), -dx1_(:,nX+2:end)];
    end

    function [f,df] = forward_constraint_fun(obj,h,x0,x1,u)
      nX = obj.plant.getNumStates();
      [xdot,dxdot] = obj.plant.dynamics(0,x0,u);
      f = x1 - x0 - h*xdot;
      df = [-xdot (-eye(nX) - h*dxdot(:,2:1+nX)) eye(nX) -h*dxdot(:,nX+2:end)];
    end
    
    function [f,df] = backward_constraint_fun(obj,h,x0,x1,u)
      nX = obj.plant.getNumStates();
      [xdot,dxdot] = obj.plant.dynamics(0,x1,u);
      f = x1 - x0 - h*xdot;
      df = [-xdot -eye(nX) (eye(nX) - h*dxdot(:,2:1+nX)) -h*dxdot(:,nX+2:end)];
    end
    
    function [f,df] = midpoint_constraint_fun(obj,h,x0,x1,u0,u1)
      nX = obj.plant.getNumStates();
      [xdot,dxdot] = obj.plant.dynamics(0,.5*(x0+x1),.5*(u0+u1));
      f = x1 - x0 - h*xdot;
      df = [-xdot (-eye(nX) - .5*h*dxdot(:,2:1+nX)) (eye(nX)- .5*h*dxdot(:,2:1+nX)) -.5*h*dxdot(:,nX+2:end) -.5*h*dxdot(:,nX+2:end)];
    end
    
    
    function [f,df] = midpoint_running_fun(obj,running_handle,h,x0,x1,u0,u1)
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
      [f,dg] = running_handle(h,.5*(x0+x1),.5*(u0+u1));
      
      df = [dg(:,1) .5*dg(:,2:1+nX) .5*dg(:,2:1+nX) .5*dg(:,2+nX:1+nX+nU) .5*dg(:,2+nX:1+nX+nU)];
    end
  end
end