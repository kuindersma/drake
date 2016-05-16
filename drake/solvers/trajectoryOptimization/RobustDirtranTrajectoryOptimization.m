classdef RobustDirtranTrajectoryOptimization < DirtranTrajectoryOptimization
  properties
    M       % number of sample disturbance per timestep
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
      if ~isfield(options,'integration_method')
        options.integration_method = DirtranTrajectoryOptimization.MIDPOINT;
      end
      obj = obj@DirtranTrajectoryOptimization(plant,N,duration,options);
      obj.nX = plant.getNumStates();
      obj.nU = plant.getNumInputs();
      obj.nW = plant.getNumDisturbances();
      obj = obj.setupRobustVariables(N,M);
      obj = obj.addDynamicConstraints;
      obj = obj.addGammaCost;
    end

    function obj = setDisturbances(obj,d)
      obj.disturbances = d;
    end
    
    function obj = setupRobustVariables(obj, N, M)
      nH = N-1;
      nX = obj.nX;
      nU = obj.nU;
      nW = obj.nU;
      nG = N-1;
      nZ = M*(N-1);
      
      num_vars = nH + N*(nX+nU+nW) + nG + nZ;
      obj.h_inds = (1:nH)';
      obj.x_inds = reshape(nH + (1:nX*N),nX,N);
      obj.u_inds = reshape(nH + nX*N + (1:nU*N),nU,N);
      obj.w_inds = reshape(nH + (nX+nU)*N + (1:nW*N),nW,N);
      obj.gamma_inds = (nH + (nX+nU+nW)*N + (1:nG))';
      obj.z_inds =     (nH + (nX+nU+nW)*N + nG + (1:nZ))';
      
      obj.M = M;
      x_names = cell(num_vars,1);
      for i = 1:N
        if(i<N)
          x_names{i} = sprintf('h[%d]',i);
          x_names{nH+(nX+nU+nW)*N+i} = sprintf('g%d[%d]',i);
          for j = 1:M
            x_names{nH+(nX+nU+nW)*N+nG+(i-1)*M+j} = sprintf('z%d[%d]',j,i);
          end
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
      end

      obj = obj.addDecisionVariable(num_vars,x_names);
    end
    
    function obj = addGammaCost(obj)
      for i=1:obj.N-1,
        running_cost = FunctionHandleObjective(1, @obj.running_gamma_cost);
        inds_i = {obj.gamma_inds(i)};
        obj = obj.addCost(running_cost,inds_i);
      end
    end
    
    function obj = addRobustConstraints(obj,running_cost,running_cost_with_w)
      nX = obj.nX;
      nU = obj.nU;
      nW = obj.nW;
      N = obj.N;
      M = obj.M;
      nZ = M*(N-1);
      
      for i=1:N-1
        for j=1:M
          % gamma lower bound constraint: gamma - ell \ge 0
          n_vars = 2*nX + 2*nU + nW + 1;
          inds = {obj.gamma_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.u_inds(:,i+1);obj.w_inds(:,i+1)};
          robust_bound_function_j = @(gamma,x0,x1,u0,u1,w1) obj.robust_constraint_fun(running_cost,running_cost_with_w,j,gamma,x0,x1,u0,u1,w1);
          constraint = FunctionHandleConstraint(0,inf,n_vars,robust_bound_function_j);
          obj = obj.addConstraint(constraint, inds);
%           % complementarity constraint: (gamma-ell)'z=0
%           n_vars = 2*nX + 2*nU + nW + 2;
%           inds = {obj.gamma_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.u_inds(:,i+1);obj.w_inds(:,i+1);obj.z_inds((i-1)*M+j)};
%           complementarity_fun_j = @(gamma,x0,x1,u0,u1,w1,zij) obj.complementarity_fun(running_cost,running_cost_with_w,j,gamma,x0,x1,u0,u1,w1,zij);
%           constraint = FunctionHandleConstraint(0,0,n_vars,complementarity_fun_j);
%           obj = obj.addConstraint(constraint, inds);
        end
          % complementarity constraint: (gamma-ell)'z=0
          n_vars = 2*nX + 2*nU + nW + 1 + M;
          inds = {obj.gamma_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.u_inds(:,i+1);obj.w_inds(:,i+1);obj.z_inds((i-1)*M+(1:M))};
          complementarity_fun = @(gamma,x0,x1,u0,u1,w1,zi) obj.complementarity_fun(running_cost,running_cost_with_w,gamma,x0,x1,u0,u1,w1,zi);
          constraint = FunctionHandleConstraint(0,0,n_vars,complementarity_fun);
          obj = obj.addConstraint(constraint, inds);

%           % complementarity constraint: (gamma-ell)'z=0
%           inds = {obj.gamma_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.u_inds(:,i+1);obj.w_inds(:,i+1);obj.z_inds((i-1)*M+(1:M))};
%           complementarity_fun = @(gamma,x0,x1,u0,u1,w1) obj.complementarity_fun(running_cost,running_cost_with_w,gamma,x0,x1,u0,u1,w1);
%           constraint = NonlinearComplementarityConstraint(complementarity_fun,2*nX + 2*nU + nW + 1,M,2,0.1);
%           obj = obj.addConstraint(constraint, inds);

          
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
    end
    
    function obj = addDynamicConstraints(obj)
      nX = obj.nX;
      nU = obj.nU;
      nW = obj.nW;
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
        case DirtranTrajectoryOptimization.MIDPOINT
          n_vars = 2*nX + 2*nU + 2*nW + 1;
          cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.midpoint_constraint_fun);
        otherwise
          error('Drake:DirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
      end
      
      for i=1:obj.N-1,
        switch obj.options.integration_method
          case RobustDirtranTrajectoryOptimization.FORWARD_EULER
            dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.w_inds(:,i)};
          case RobustDirtranTrajectoryOptimization.BACKWARD_EULER
            dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.w_inds(:,i)};
          case DirtranTrajectoryOptimization.MIDPOINT
            dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.u_inds(:,i+1);obj.w_inds(:,i);obj.w_inds(:,i+1)};
          otherwise
            error('Drake:DirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
        end
        constraints{i} = cnstr;
        
        obj = obj.addConstraint(constraints{i}, dyn_inds{i});
      end
    end
    
    function [f,df] = forward_constraint_fun(obj,h,x0,x1,u,w)
      [xdot,dxdot] = obj.plant.dynamics_w(0,x0,u,w);
      f = x1 - x0 - h*xdot;
      df = [-xdot (-eye(obj.nX) - h*dxdot(:,2:1+obj.nX)) eye(obj.nX) -h*dxdot(:,obj.nX+2:end)];
    end
    
    function [f,df] = backward_constraint_fun(obj,h,x0,x1,u,w)
      [xdot,dxdot] = obj.plant.dynamics_w(0,x1,u,w);
      f = x1 - x0 - h*xdot;
      df = [-xdot -eye(obj.nX) (eye(obj.nX) - h*dxdot(:,2:1+obj.nX)) -h*dxdot(:,obj.nX+2:end)];
    end
    
    function [f,df] = midpoint_constraint_fun(obj,h,x0,x1,u0,u1,w0,w1)
      [xdot,dxdot] = obj.plant.dynamics_w(0,.5*(x0+x1),.5*(u0+u1),.5*(w0+w1));
      f = x1 - x0 - h*xdot;
      df = [-xdot (-eye(obj.nX) - .5*h*dxdot(:,1+(1:obj.nX))) (eye(obj.nX) - .5*h*dxdot(:,1+(1:obj.nX))) -.5*h*dxdot(:,1+obj.nX+(1:obj.nU)) -.5*h*dxdot(:,1+obj.nX+(1:obj.nU))  -.5*h*dxdot(:,1+obj.nX+obj.nU+(1:obj.nW)) -.5*h*dxdot(:,1+obj.nX+obj.nU+(1:obj.nW))];
    end


    function [f,df] = robust_constraint_fun(obj,running_cost,running_cost_with_w,j,gamma,x0,x1,u0,u1,w1)
      [gw0,dgw0] = running_cost_with_w(0,x0,u0,obj.disturbances(:,j));
      [gw1,dgw1] = running_cost_with_w(0,x1,u1,w1);
      [g0,dg0] = running_cost(0,x0,u0);
      [g1,dg1] = running_cost(0,x1,u1);
      f = gamma - gw0 - gw1 + g0 + g1;
      df = [1, ...
            -dgw0(1+(1:obj.nX))+dg0(1+(1:obj.nX)), ...
            -dgw1(1+(1:obj.nX))+dg1(1+(1:obj.nX)), ...
            -dgw0(1+obj.nX+(1:obj.nU))+dg0(1+obj.nX+(1:obj.nU)), ...
            -dgw1(1+obj.nX+(1:obj.nU))+dg1(1+obj.nX+(1:obj.nU)), ...
            -dgw1(1+obj.nX+obj.nU+(1:obj.nW))];
    end

%     function [f,df] = complementarity_fun(obj,running_cost,running_cost_with_w,j,gamma,x0,x1,u0,u1,w1,zij)
%       [g,dg] = robust_constraint_fun(obj,running_cost,running_cost_with_w,j,gamma,x0,x1,u0,u1,w1);
%       f = g'*zij;
%       df = [dg*zij, g];
%     end

    function [f,df] = complementarity_fun(obj,running_cost,running_cost_with_w,gamma,x0,x1,u0,u1,w1,zi)
      v = zeros(obj.M,1);
      dv = [];
      for j=1:obj.M
        [g,dg] = robust_constraint_fun(obj,running_cost,running_cost_with_w,j,gamma,x0,x1,u0,u1,w1);
        v(j) = g;
        dv = [dv;dg];
      end
      f = v'*zi;
      df = [zi'*dv, v'];
    end

%     function [f,df] = complementarity_fun(obj,running_cost,running_cost_with_w,gamma,x0,x1,u0,u1,w1)
%       f = zeros(obj.M,1);
%       df = [];
%       for j=1:obj.M
%         [g,dg] = robust_constraint_fun(obj,running_cost,running_cost_with_w,j,gamma,x0,x1,u0,u1,w1);
%         f(j) = g;
%         df = [df;dg];
%       end
%     end
    
    function [f,df] = running_gamma_cost(~,gamma)
      f = gamma; % scalar
      df = 1;
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