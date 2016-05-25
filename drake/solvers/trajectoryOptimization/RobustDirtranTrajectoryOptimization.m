classdef RobustDirtranTrajectoryOptimization < DirtranTrajectoryOptimization
  properties
    M       % number of sample disturbance per timestep
    dx_inds  % n x N indices for delta state in disturbed trajectory
    du_inds  % m x N indices for delta inputs in disturbed trajectory
    gamma_inds  % N-1 x 1 indices for uppper bound on cost gain
    % of the sampled disturbances incurs the highest cost gain at each knot
    % point    
    disturbances
    nX
    nU
    nW
    K
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
    
    function obj = setDisturbances(obj,d)
      obj.disturbances = d;
    end
    
    function obj = setupRobustVariables(obj, N, M)
      nH = N-1;
      nX = obj.nX;
      nU = obj.nU;
      nG = N-1;
      ndX = nX*M;
      ndU = nU*M;
      
      num_vars = nH + N*(nX+ndX+nU+ndU) + nG;
      obj.h_inds = (1:nH)';
      obj.x_inds = reshape(nH + (1:nX*N),nX,N);
      obj.u_inds = reshape(nH + nX*N + (1:nU*N),nU,N);
      obj.dx_inds = reshape(nH + (nX+nU)*N +(1:ndX*N),ndX,N);
      obj.du_inds = reshape(nH + (nX+ndX+nU)*N +(1:ndU*N),ndU,N);
      obj.gamma_inds = (nH + (nX+ndX+nU+ndU)*N + (1:nG))';
      
      x_names = cell(num_vars,1);
      for i = 1:N
        if(i<N)
          x_names{i} = sprintf('h[%d]',i);
          x_names{nH+(nX+ndX+nU+ndU)*N+i} = sprintf('g[%d]',i);
        end
        for j = 1:nX
          x_names{nH+(i-1)*nX+j}=sprintf('x%d[%d]',j,i);
        end
        for j = 1:nU
          x_names{nH+nX*N+(i-1)*nU+j} = sprintf('u%d[%d]',j,i);
        end
        for j = 1:ndX
          x_names{nH+(nX+nU)*N+(i-1)*ndX+j} = sprintf('dx%d[%d]',j,i);
        end
        for j = 1:ndU
          x_names{nH+(nX+ndX+nU)*N+(i-1)*ndU+j} = sprintf('du%d[%d]',j,i);
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
      N = obj.N;
      M = obj.M;
      
      for i=1:N-1
        for j=1:M
          switch obj.options.integration_method
            case DirtranTrajectoryOptimization.FORWARD_EULER
              % set delta_x1 = x1 - f(x0,u0+du0,w)*h - x0 
              n_vars = 1 + nX*M + 3*nX + 2*nU;
              inds = {obj.h_inds(i); ...
                      obj.x_inds(:,i); ...
                      obj.dx_inds(:,i); ...
                      obj.x_inds(:,i+1); ...
                      obj.dx_inds((j-1)*nX+(1:nX),i+1); ...
                      obj.u_inds(:,i); ...
                      obj.du_inds((j-1)*nU+(1:nU),i)};
             
                    
              forward_delta_x_constraint_j = @(h,x0,dx0,x1,dx1j,u,du) obj.forward_delta_x_constraint(j,h,x0,dx0,x1,dx1j,u,du);

              constraint = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,forward_delta_x_constraint_j);
              obj = obj.addConstraint(constraint, inds);
              
             
            case DirtranTrajectoryOptimization.BACKWARD_EULER
              error('not implemented yet');
           
            case DirtranTrajectoryOptimization.MIDPOINT
              error('not implemented yet');
            
            otherwise
              error('Drake:DirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
          end
          
          % gamma lower bound constraint: gamma - ell \ge 0
          n_vars = 1 + nX + nU;
          inds = {obj.gamma_inds(i); ...
                  obj.dx_inds((j-1)*nX+(1:nX),i+1); ...
                  obj.du_inds((j-1)*nU+(1:nU),i)};
              
          robust_bound_fun_j = @(gamma,dx,du) obj.robust_bound_fun(robust_cost,j,gamma,dx,du);

          constraint = FunctionHandleConstraint(0,inf,n_vars,robust_bound_fun_j);
          obj = obj.addConstraint(constraint, inds);
              
          
          if any(~isinf(obj.plant.umin)) || any(~isinf(obj.plant.umax))
            inds = {obj.u_inds(:,i); obj.du_inds((j-1)*nU+(1:nU),i)};
            constraint = FunctionHandleConstraint(obj.plant.umin,obj.plant.umax,2*nU,@obj.delta_u_bound);
            obj = obj.addConstraint(constraint, inds);  
          end        
        end
        
        
        
      end
      
      obj = obj.addConstraint(ConstantConstraint(zeros(nX,1)),obj.dx_inds(1:nX,1)); % first dx=0
      obj = obj.addConstraint(ConstantConstraint(zeros(nU,1)),obj.u_inds(:,N)); % last u=0
    end
    

    function obj = addStateRegularization(obj)

      for i=1:obj.N-1
        running_cost = FunctionHandleObjective(2*obj.nX, @obj.state_regularizer);
        inds = {obj.x_inds(:,i);obj.x_inds(:,i+1)};
        obj = obj.addCost(running_cost,inds);
      end        
      
    end
    
    
    function obj = addDeltaXEqualsZeroConstraint(obj)
      obj = obj.addConstraint(ConstantConstraint(zeros(obj.nX*obj.M*obj.N,1)),obj.dx_inds(:)); 
    end
  
    
    
    function obj = addLinearControlConstraint(obj,K)
      obj.K = K;
      n_vars = obj.nX + obj.nU;
      for i=1:obj.N-1
        for j=1:obj.M
          inds = {obj.dx_inds((j-1)*obj.nX+(1:obj.nX),i); ...
                  obj.du_inds((j-1)*obj.nU+(1:obj.nU),i)};

          constraint = FunctionHandleConstraint(zeros(obj.nU,1),zeros(obj.nU,1),n_vars,@obj.linear_delta_u_constraint);
          obj = obj.addConstraint(constraint, inds);
        end
      end
    end
    
    
    
    function [f,df] = forward_robust_dynamics_fun(obj,h,x0,u,du,w)
      [xdot,dxdot] = obj.plant.dynamics_w(0,x0,u+du,w);
      f = x0 + h*xdot;
      df = [xdot ... h
        (eye(obj.nX) + h*dxdot(:,2:1+obj.nX)) ... x0
        h*dxdot(:,1+obj.nX+(1:obj.nU)) ... u
        h*dxdot(:,1+obj.nX+(1:obj.nU)) ... du
        h*dxdot(:,1+obj.nX+obj.nU+(1:obj.nW))]; % w
    end
    
%     function [f,df] = backward_robust_dynamics_fun(obj,h,x0,x1,u,du,w)
%       [xdot,dxdot] = obj.plant.dynamics_w(0,x1,u+du,w);
%       f = x0 + h*xdot;
%       df = [xdot ... h
%         eye(obj.nX) ... x0
%         h*dxdot(:,2:1+obj.nX) ... x1
%         h*dxdot(:,1+obj.nX+(1:obj.nU)) ... u
%         h*dxdot(:,1+obj.nX+(1:obj.nU)) ... du
%         h*dxdot(:,1+obj.nX+obj.nU+(1:obj.nW))]; % w
%     end
    
%     function [f,df] = midpoint_robust_dynamics_fun(obj,h,x0,x1,u0,du0,u1,du1,w0,w1)
%       [xdot,dxdot] = obj.plant.dynamics_w(0,.5*(x0+x1),.5*(u0+du0+u1+du1),.5*(w0+w1));
%       f = x0 + h*xdot;
%       df = [xdot ... h
%         (eye(obj.nX) + .5*h*dxdot(:,1+(1:obj.nX))) ... x0
%         .5*h*dxdot(:,1+(1:obj.nX))... x1
%         .5*h*dxdot(:,1+obj.nX+(1:obj.nU)) ... u0
%         .5*h*dxdot(:,1+obj.nX+(1:obj.nU)) ... du0
%         .5*h*dxdot(:,1+obj.nX+(1:obj.nU)) ... u1
%         .5*h*dxdot(:,1+obj.nX+(1:obj.nU)) ... du1
%         .5*h*dxdot(:,1+obj.nX+obj.nU+(1:obj.nW)) ... w0
%         .5*h*dxdot(:,1+obj.nX+obj.nU+(1:obj.nW))]; % w1
%     end

    
    function [f,df] = linear_delta_u_constraint(obj,dx,du)
      f = du - obj.K*dx;

      df = [-obj.K, ... dx
            eye(obj.nU)]; % du 
    end
    
    function [f,df] = robust_bound_fun(obj,robust_cost,j,gamma,dx,du)
      w = obj.disturbances(:,j);
      
      [g,dg] = robust_cost(dx,du,w);
      f = gamma - g;
      
      df = [1, ... gamma
            -dg(1:obj.nX), ... dx
            -dg(obj.nX+(1:obj.nU))]; % du 
    end

    
    function [f,df] = forward_delta_x_constraint(obj,j,h,x0,dx0,x1,dx1j,u,du)
      w = obj.disturbances(:,j);
      
      dx0_max = zeros(obj.nX,1);
      for j=1:obj.M
         dx0j = dx0((j-1)*obj.nX+(1:obj.nX));
         if norm(dx0j)>=norm(dx0_max)
            dx0_max = dx0j;
            df_ddx0 = zeros(obj.nX,obj.nX*obj.M);
            df_ddx0(:,(j-1)*obj.nX+(1:obj.nX)) = eye(obj.nX);
         end        
      end
      
      [r,dr] = obj.forward_robust_dynamics_fun(h,x0,u,du,w);
    
      f = dx0_max + dx1j - r + x1; % delta_x(1) = max(delta_x(0)) + x0 + f(x0,u0+du0,w) - x(1)
      
      dr_dh = dr(:,1);
      dr_dx0 = dr(:,1+(1:obj.nX));
      dr_du = dr(:,1+obj.nX+(1:obj.nU));
      dr_ddu = dr(:,1+obj.nX+obj.nU+(1:obj.nU));
      
      df_dh = -dr_dh;
      df_dx0 = -dr_dx0;
      df_dx1 = eye(obj.nX);
      df_ddx1 = eye(obj.nX);
      df_du = -dr_du;
      df_ddu = -dr_ddu;
      
      df = [df_dh, df_dx0, df_ddx0 df_dx1, df_ddx1, df_du, df_ddu];
    end
    
    
%     function [f,df] = robust_backward_bound_fun(obj,robust_cost,j,gamma,h,x0,x1,u,du)
%       w = obj.disturbances(:,j);
%       
%       [r,dr] = obj.backward_robust_dynamics_fun(h,x0,x1,u,du,w);
%     
%       xerr = x1-r;
%       [g,dg] = robust_cost(xerr,du,w);
%       f = gamma - g;
%       
%       dr_dh = dr(:,1);
%       dr_dx0 = dr(:,1+(1:obj.nX));
%       dr_dx1 = dr(:,1+obj.nX+(1:obj.nX));
%       dr_du = dr(:,1+2*obj.nX+(1:obj.nU));
%       dr_ddu = dr(:,1+2*obj.nX+obj.nU+(1:obj.nU));
%       
%       dxerr_dh = -dr_dh;
%       dxerr_dx0 = -dr_dx0;
%       dxerr_dx1 = eye(obj.nX) - dr_dx1;
%       dxerr_du = -dr_du;
%       dxerr_ddu = -dr_ddu;
%       
%       df = [1, ... gamma
%             -dg((1:obj.nX))*dxerr_dh, ... h
%             -dg((1:obj.nX))*dxerr_dx0, ... x0
%             -dg((1:obj.nX))*dxerr_dx1, ... x1
%             -dg((1:obj.nX))*dxerr_du, ... u
%             -dg((1:obj.nX))*dxerr_ddu-dg(obj.nX + (1:obj.nU))]; %ddu
%     end
    
    
       
%     function [f,df] = robust_midpoint_bound_fun(obj,robust_cost,j,gamma,h,x0,x1,u0,du0,u1)
%       w = obj.disturbances(:,j);
%       
%       [r,dr] = obj.midpoint_robust_dynamics_fun(h,x0,x1,u0,du0,u1,0*u1,w,0*w);
%     
%       xerr = x1-r;
%       [g,dg] = robust_cost(xerr,du0,w);
%       f = gamma - g;
%           
%       dr_dh = dr(:,1);
%       dr_dx0 = dr(:,1+(1:obj.nX));
%       dr_dx1 = dr(:,1+obj.nX+(1:obj.nX));
%       dr_du0 = dr(:,1+2*obj.nX+(1:obj.nU));
%       dr_ddu0 = dr(:,1+2*obj.nX+obj.nU+(1:obj.nU));
%       dr_du1 = dr(:,1+2*obj.nX+2*obj.nU+(1:obj.nU));
%       
%       dxerr_dh = -dr_dh;
%       dxerr_dx0 = -dr_dx0;
%       dxerr_dx1 = eye(obj.nX) - dr_dx1;
%       dxerr_du0 = -dr_du0;
%       dxerr_ddu0 = -dr_ddu0;
%       dxerr_du1 = -dr_du1;
%       
%       df = [1, ... gamma
%             -dg((1:obj.nX))*dxerr_dh, ... h
%             -dg((1:obj.nX))*dxerr_dx0, ... x0
%             -dg((1:obj.nX))*dxerr_dx1, ... x1
%             -dg((1:obj.nX))*dxerr_du0, ... u0
%             -dg((1:obj.nX))*dxerr_ddu0-dg(obj.nX + (1:obj.nU)) ... du0
%             -dg((1:obj.nX))*dxerr_du1]; % du1
%     end
    
    
    function [f,df] = gamma_cost(obj,gamma)
      f = ones(1,obj.N-1)*gamma; 
      df = ones(1,obj.N-1);
    end

    function [f,df] = state_regularizer(obj,x0,x1)
      nq = obj.nX/2;
      Q = 50*diag([ones(1,nq),0.1*ones(1,nq)]);
      xerr = x1-x0;
      f = xerr'*Q*xerr; 
      df = [-2*xerr'*Q,2*xerr'*Q];
    end

    function [f,df] = delta_u_bound(obj,u,du)
      f = u+du;
      df = [eye(obj.nU), eye(obj.nU)];
    end
    
  end
end