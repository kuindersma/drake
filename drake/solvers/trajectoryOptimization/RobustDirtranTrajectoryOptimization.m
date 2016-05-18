classdef RobustDirtranTrajectoryOptimization < DirtranTrajectoryOptimization
  properties
    M       % number of sample disturbance per timestep
    xr_inds  % n x N indices for states in disturbed trajectory
    w_inds  % d x N indices for disturbances
    gamma_inds  % N-1 x 1 indices for uppper bound on cost gain
    z_inds  % M x N indices for indicator variables on disturbances---identifies which
    % of the sampled disturbances incurs the highest cost gain at each knot
    % point    
    disturbances
    nX
    nU
    nW
  end
  
  methods
    function obj = RobustDirtranTrajectoryOptimization(plant,N,M,duration,options)
      if nargin < 5
        options = struct();
      end
      if isscalar(duration), duration=[duration,duration]; end

      if ~isfield(options,'integration_method')
        options.integration_method = DirtranTrajectoryOptimization.FORWARD_EULER;
      end
      obj = obj@DirtranTrajectoryOptimization(plant,N,duration,options);
      obj.nX = plant.getNumStates();
      obj.nU = plant.getNumInputs();
      obj.nW = plant.getNumDisturbances();
      obj.N = N;
      obj.M = M;
      obj = obj.setupRobustVariables(N,M);
      obj = obj.addDynamicConstraints;
      obj = obj.addRobustDynamicConstraints;
      obj = obj.addGammaCost;
    
      if ~isfield(options,'time_option')
        options.time_option = 1;
      end
      
      
      % Construct total time linear constraint
      switch options.time_option
        case 1 % all timesteps are constant
          A_time = [ones(1,N-1);[eye(N-2) zeros(N-2,1)] - [zeros(N-2,1) eye(N-2)]];
          time_constraint = LinearConstraint([duration(1);zeros(N-2,1)],[duration(2);zeros(N-2,1)],A_time);
          obj = obj.addConstraint(time_constraint,obj.h_inds);
        case 2 % all timesteps independent
          A_time = ones(1,N-1);
          time_constraint = LinearConstraint(duration(1),duration(2),A_time);
          obj = obj.addConstraint(time_constraint,obj.h_inds);
      end

      % Ensure that all h values are non-negative
      obj = obj.addConstraint(BoundingBoxConstraint(zeros(N-1,1),inf(N-1,1)),obj.h_inds);

      % add control inputs as bounding box constraints
      if any(~isinf(plant.umin)) || any(~isinf(plant.umax))
        control_limit = BoundingBoxConstraint(repmat(plant.umin,N,1),repmat(plant.umax,N,1));
        obj = obj.addConstraint(control_limit,obj.u_inds(:));
      end

      
    end

    function obj = addRobustStateConstraint(obj,constraint,time_index,x_indices)
      if ~iscell(time_index)
        % then use { time_index(1), time_index(2), ... } ,
        % aka independent constraints for each time
        time_index = num2cell(reshape(time_index,1,[]));
      end
      if nargin<4, x_indices = 1:size(obj.xr_inds,1); end

      for j=1:length(time_index),
        cstr_inds = mat2cell(obj.xr_inds(x_indices,time_index{j}),numel(x_indices),ones(1,length(time_index{j})));

        % record constraint for posterity
        obj.constraints{end+1}.constraint = constraint;
        obj.constraints{end}.var_inds = cstr_inds;
        obj.constraints{end}.time_index = time_index;

        obj = obj.addConstraint(constraint,cstr_inds);
      end
    end

    
    function obj = setDisturbances(obj,d)
      obj.disturbances = d;
    end
    
    function obj = setupRobustVariables(obj, N, M)
      nH = N-1;
      nX = obj.nX;
      nU = obj.nU;
      nW = obj.nW;
      nG = N-1;
      nZ = M*(N-1);
      
      num_vars = nH + N*(2*nX+nU+nW) + nG + nZ;
      obj.h_inds = (1:nH)';
      obj.x_inds = reshape(nH + (1:nX*N),nX,N);
      obj.u_inds = reshape(nH + nX*N + (1:nU*N),nU,N);
      obj.xr_inds = reshape(nH + nX*N+ nU*N+(1:nX*N),nX,N);
      obj.w_inds = reshape(nH + (2*nX+nU)*N + (1:nW*N),nW,N);
      obj.gamma_inds = (nH + (2*nX+nU+nW)*N + (1:nG))';
      obj.z_inds =     (nH + (2*nX+nU+nW)*N + nG + (1:nZ))';
      
      obj.M = M;
      x_names = cell(num_vars,1);
      for i = 1:N
        if(i<N)
          x_names{i} = sprintf('h[%d]',i);
          x_names{nH+(2*nX+nU+nW)*N+i} = sprintf('g[%d]',i);
          for j = 1:M
            x_names{nH+(2*nX+nU+nW)*N+nG+(i-1)*M+j} = sprintf('z%d[%d]',j,i);
          end
        end
        for j = 1:nX
          x_names{nH+(i-1)*nX+j}=sprintf('x%d[%d]',j,i);
        end
        for j = 1:nU
          x_names{nH+nX*N+(i-1)*nU+j} = sprintf('u%d[%d]',j,i);
        end
        for j = 1:nX
          x_names{nH+nX*N+nU*N+(i-1)*nX+j}=sprintf('xr%d[%d]',j,i);
        end
        for j = 1:nW
          x_names{nH+(2*nX+nU)*N+(i-1)*nW+j} = sprintf('w%d[%d]',j,i);
        end
      end

      obj = obj.addDecisionVariable(num_vars,x_names);
    end
    
    function obj = addGammaCost(obj)
      running_cost = FunctionHandleObjective(obj.N-1, @obj.gamma_cost);
      inds_i = {obj.gamma_inds};
      obj = obj.addCost(running_cost,inds_i);
    end
    
    function obj = addRobustConstraints(obj,robust_cost)
      nX = obj.nX;
      nU = obj.nU;
      nW = obj.nW;
      N = obj.N;
      M = obj.M;
      nZ = M*(N-1);
      
      for i=1:N-1
        for j=1:M
          % gamma lower bound constraint: gamma - ell \ge 0
          n_vars = 4*nX + nW + 1;
          inds = {obj.gamma_inds(i);obj.x_inds(:,i);obj.xr_inds(:,i);obj.x_inds(:,i+1);obj.xr_inds(:,i+1);obj.w_inds(:,i+1)};
          robust_bound_function_j = @(gamma,x0,xr0,x1,xr1,w1) obj.robust_bound_fun(robust_cost,j,gamma,x0,xr0,x1,xr1,w1);
          constraint = FunctionHandleConstraint(0,inf,n_vars,robust_bound_function_j);
          obj = obj.addConstraint(constraint, inds);
        end
        % complementarity constraint: (gamma-ell)'z=0
        n_vars = 4*nX + nW + 1 + M;
        inds = {obj.gamma_inds(i);obj.x_inds(:,i);obj.xr_inds(:,i);obj.x_inds(:,i+1);obj.xr_inds(:,i+1);obj.w_inds(:,i+1);obj.z_inds((i-1)*M+(1:M))};
        complementarity_fun = @(gamma,x0,xr0,x1,xr1,w1,zi) obj.complementarity_fun(robust_cost,gamma,x0,xr0,x1,xr1,w1,zi);
        constraint = FunctionHandleConstraint(0,0,n_vars,complementarity_fun);
        obj = obj.addConstraint(constraint, inds);
        % Sum constraint on z: \sum_j zij = 1
        inds = {obj.z_inds((i-1)*M+(1:M))};
        constraint = FunctionHandleConstraint(1,1,M,@obj.z_sum_constr);
        obj = obj.addConstraint(constraint, inds);
        % equality constraint on wi: wi - \sum_j zij = 0
        inds = {obj.w_inds(:,i);obj.z_inds((i-1)*M+(1:M))};
        constraint = FunctionHandleConstraint(zeros(nW,1),zeros(nW,1),nW+M,@obj.w_equality);
        obj = obj.addConstraint(constraint, inds);        
      end
      obj = obj.addConstraint(BoundingBoxConstraint(zeros(nZ,1),ones(nZ,1)),obj.z_inds);
      obj = obj.addConstraint(ConstantConstraint(zeros(nW,1)),obj.w_inds(:,N));
      obj = obj.addConstraint(ConstantConstraint(zeros(nU,1)),obj.u_inds(:,N));
    end
    
    function obj = addRobustDynamicConstraints(obj)
      nX = obj.nX;
      nU = obj.nU;
      nW = obj.nW;
      N = obj.N;

      switch obj.options.integration_method
        case DirtranTrajectoryOptimization.FORWARD_EULER
          n_vars = 2*nX + nU + nW + 1;
          cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.forward_robust_dynamics_fun);
        case DirtranTrajectoryOptimization.BACKWARD_EULER
          n_vars = 2*nX + nU + nW + 1;
          cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.backward_robust_dynamics_fun);
        case DirtranTrajectoryOptimization.MIDPOINT
          n_vars = 2*nX + 2*nU + 2*nW + 1;
          cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.midpoint_robust_dynamics_fun);
        otherwise
          error('Drake:DirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
      end
      
      for i=1:N-1,
        switch obj.options.integration_method
          case DirtranTrajectoryOptimization.FORWARD_EULER
            dyn_inds = {obj.h_inds(i);obj.xr_inds(:,i);obj.xr_inds(:,i+1);obj.u_inds(:,i);obj.w_inds(:,i)};
          case DirtranTrajectoryOptimization.BACKWARD_EULER
            dyn_inds = {obj.h_inds(i);obj.xr_inds(:,i);obj.xr_inds(:,i+1);obj.u_inds(:,i);obj.w_inds(:,i)};
          case DirtranTrajectoryOptimization.MIDPOINT
            dyn_inds = {obj.h_inds(i);obj.xr_inds(:,i);obj.xr_inds(:,i+1);obj.u_inds(:,i);obj.u_inds(:,i+1);obj.w_inds(:,i);obj.w_inds(:,i+1)};
          otherwise
            error('Drake:DirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
        end

        obj = obj.addConstraint(cnstr, dyn_inds);
      end
    end
    
    
    function [f,df] = forward_robust_dynamics_fun(obj,h,x0,x1,u,w)
      [xdot,dxdot] = obj.plant.dynamics_w(0,x0,u,w);
      f = x1 - x0 - h*xdot;
      df = [-xdot (-eye(obj.nX) - h*dxdot(:,2:1+obj.nX)) eye(obj.nX) -h*dxdot(:,obj.nX+2:end)];
    end
    
    function [f,df] = backward_robust_dynamics_fun(obj,h,x0,x1,u,w)
      [xdot,dxdot] = obj.plant.dynamics_w(0,x1,u,w);
      f = x1 - x0 - h*xdot;
      df = [-xdot -eye(obj.nX) (eye(obj.nX) - h*dxdot(:,2:1+obj.nX)) -h*dxdot(:,obj.nX+2:end)];
    end
    
    function [f,df] = midpoint_robust_dynamics_fun(obj,h,x0,x1,u0,u1,w0,w1)
      [xdot,dxdot] = obj.plant.dynamics_w(0,.5*(x0+x1),.5*(u0+u1),.5*(w0+w1));
      f = x1 - x0 - h*xdot;
      df = [-xdot (-eye(obj.nX) - .5*h*dxdot(:,1+(1:obj.nX))) (eye(obj.nX) - .5*h*dxdot(:,1+(1:obj.nX))) -.5*h*dxdot(:,1+obj.nX+(1:obj.nU)) -.5*h*dxdot(:,1+obj.nX+(1:obj.nU))  -.5*h*dxdot(:,1+obj.nX+obj.nU+(1:obj.nW)) -.5*h*dxdot(:,1+obj.nX+obj.nU+(1:obj.nW))];
    end


    function [f,df] = robust_bound_fun(obj,robust_cost,j,gamma,x0,xr0,x1,xr1,w1)
      [g0,dg0] = robust_cost(x0,xr0,obj.disturbances(:,j));
      [g1,dg1] = robust_cost(x1,xr1,w1);
      f = gamma - g0 - g1;
      df = [1, ...
            -dg0((1:obj.nX)), ...
            -dg0(obj.nX+(1:obj.nX)), ...
            -dg1((1:obj.nX)), ...
            -dg1(obj.nX+(1:obj.nX)), ...
            -dg1(2*obj.nX+(1:obj.nW))];
    end

    function [f,df] = complementarity_fun(obj,robust_cost,gamma,x0,xr0,x1,xr1,w1,zi)
      v = zeros(obj.M,1);
      dv = [];
      for j=1:obj.M
        [g,dg] = robust_bound_fun(obj,robust_cost,j,gamma,x0,xr0,x1,xr1,w1);
        v(j) = g;
        dv = [dv;dg];
      end
      f = v'*zi;
      df = [zi'*dv, v'];
    end
   
    function [f,df] = gamma_cost(obj,gamma)
      f = -ones(1,obj.N-1)*gamma; 
      df = -ones(1,obj.N-1);
    end
    
    function [f,df] = z_sum_constr(~,zi)
      f = ones(1,length(zi))*zi;
      df = ones(1,length(zi));
    end

    function [f,df] = w_equality(obj,wi,zi)
      f = wi - obj.disturbances*zi;
      df = [eye(obj.nW), -obj.disturbances];
    end

    
  end
end