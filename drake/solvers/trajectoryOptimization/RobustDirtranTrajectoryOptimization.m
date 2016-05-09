classdef RobustDirtranTrajectoryOptimization < DirtranTrajectoryOptimization
  properties
    M       % number of sample disturbance per timestep
    w_inds  % d x N indices for disturbances
    gamma_inds  % N-1 x 1 indices for uppper bound on cost gain
    z_inds  % M x N indices for indicator variables on disturbances---identifies which
    % of the sampled disturbances incurs the highest cost gain at each knot
    % point    
  end
  
  methods
    function obj = RobustDirtranTrajectoryOptimization(plant,N,M,duration,options)
      if nargin < 4
        options = struct();
      end
      if ~isfield(options,'integration_method')
        options.integration_method = DirtranTrajectoryOptimization.FORWARD_EULER;
      end
      obj = obj@DirtranTrajectoryOptimization(plant,N,duration,options);
      obj = obj.setupRobustVariables(N,M);
      obj = obj.addDynamicConstraints;
      obj = obj.addRunningGammaCost;
    end
    
    function obj = setupVariables(obj,~)
      return;
    end
    
    function obj = setupRobustVariables(obj, N, M)
      nH = N-1;
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
      nW = obj.plant.getNumDisturbances();
      nG = N;
      nZ = M*N;
      
      num_vars = nH + N*(nX+nU+nW) + nG + nZ;
      obj.h_inds = (1:nH)';
      obj.x_inds = reshape(nH + (1:nX*N),nX,N);
      obj.u_inds = reshape(nH + nX*N + (1:nU*N),nU,N);
      obj.w_inds = reshape(nH + (nX+nU)*N + (1:nW*N),nW,N);
      obj.gamma_inds = (nH + (nX+nU+nW)*N + (1:nG))';
      obj.z_inds = reshape(nH + (nX+nU)*N + nW*N + (1:nZ),M,N);
      
      obj.M = M;
      x_names = cell(num_vars,1);
      for i = 1:N
        if(i<N)
          x_names{i} = sprintf('h[%d]',i);
        end
        for j = 1:nX
          x_names{nH+(i-1)*nX+j}=sprintf('x%d[%d]',j,i);
        end
        for j = 1:nU
          x_names{nH+nX*N+(i-1)*nU+j} = sprintf('u%d[%d]',j,i);
        end
        for j = 1:nW
          x_names{nH+(nX+nU)*N+(i-1)*nW+j} = sprintf('w%d[%d]',j,i);
        end
        for j = 1:nG
          x_names{nH+(nX+nU+nW)*N+(i-1)*nG+j} = sprintf('g%d[%d]',j,i);
        end
        for j = 1:nZ
          x_names{nH+(nX+nU+nW)*N+nG+(i-1)*nZ+j} = sprintf('z%d[%d]',j,i);
        end
      end

      obj = obj.addDecisionVariable(num_vars,x_names);
    end
    
    function obj = addRobustCost(obj)
      for i=1:obj.N-1,
        running_cost = FunctionHandleObjective(1, @running_gamma_cost);
        inds_i = {obj.gamma_inds(i)};
        obj = obj.addCost(running_cost,inds_i);
      end
    end
    
    function obj = addRobustBoundConstraint(obj,running_cost_fun_w,running_cost_fun)
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
      nW = obj.plant.getNumDisturbances();
      N = obj.N;
      M = obj.M;
      
      constraints = cell(M*(N-1),1);
      dyn_inds = cell(M*(N-1),1);
      
      switch obj.options.integration_method
        case RobustDirtranTrajectoryOptimization.FORWARD_EULER
          n_vars = 2*nX + nU + nW + 1;
          cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.forward_constraint_fun);
        case RobustDirtranTrajectoryOptimization.BACKWARD_EULER
          n_vars = 2*nX + nU + nW + 1;
          cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.backward_constraint_fun);
        otherwise
          error('Drake:DirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
      end
      
      for i=1:obj.N-1,
        switch obj.options.integration_method
          case RobustDirtranTrajectoryOptimization.FORWARD_EULER
            dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.w_inds(:,i)};
          case RobustDirtranTrajectoryOptimization.BACKWARD_EULER
            dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.w_inds(:,i)};
          otherwise
            error('Drake:DirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
        end
        constraints{i} = cnstr;
        
        obj = obj.addConstraint(constraints{i}, dyn_inds{i});
      end
    end
    
    function obj = addDynamicConstraints(obj)
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
      nW = obj.plant.getNumDisturbances();
      N = obj.N;
      
      constraints = cell(N-1,1);
      dyn_inds = cell(N-1,1);
      
      switch obj.options.integration_method
        case RobustDirtranTrajectoryOptimization.FORWARD_EULER
          n_vars = 2*nX + nU + nW + 1;
          cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.forward_constraint_fun);
        case RobustDirtranTrajectoryOptimization.BACKWARD_EULER
          n_vars = 2*nX + nU + nW + 1;
          cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.backward_constraint_fun);
        otherwise
          error('Drake:DirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
      end
      
      for i=1:obj.N-1,
        switch obj.options.integration_method
          case RobustDirtranTrajectoryOptimization.FORWARD_EULER
            dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.w_inds(:,i)};
          case RobustDirtranTrajectoryOptimization.BACKWARD_EULER
            dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.w_inds(:,i)};
          otherwise
            error('Drake:DirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
        end
        constraints{i} = cnstr;
        
        obj = obj.addConstraint(constraints{i}, dyn_inds{i});
      end
    end
  end
  
  methods (Access=protected)
    function [f,df] = forward_constraint_fun(obj,h,x0,x1,u,w)
      nX = obj.plant.getNumStates();
      [xdot,dxdot] = obj.plant.dynamics(0,x0,u,w);
      f = x1 - x0 - h*xdot;
      df = [-xdot (-eye(nX) - h*dxdot(:,2:1+nX)) eye(nX) -h*dxdot(:,nX+2:end)];
    end
    
    function [f,df] = backward_constraint_fun(obj,h,x0,x1,u,w)
      nX = obj.plant.getNumStates();
      [xdot,dxdot] = obj.plant.dynamics(0,x1,u,w);
      f = x1 - x0 - h*xdot;
      df = [-xdot -eye(nX) (eye(nX) - h*dxdot(:,2:1+nX)) -h*dxdot(:,nX+2:end)];
    end
    
    function [f,df] = running_gamma_cost(~,gamma)
      f = gamma; % scalar
      df = 1;
    end

  end
end