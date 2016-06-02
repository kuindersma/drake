classdef RobustDirtranTrajectoryOptimization < DirtranTrajectoryOptimization
  properties
    M       % number of sample disturbance per timestep
    dx_inds  % n x N indices for delta state in disturbed trajectory
    du_inds  % m x N indices for delta inputs in disturbed trajectory
    w_inds  % d x N indices for disturbances
    gamma_inds  % N-1 x 1 indices for uppper bound on cost gain
    z_inds  % M x N indices for indicator variables on disturbances---identifies which
    % of the sampled disturbances incurs the highest cost gain at each knot point    
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
%       obj = obj.addZCost;
    
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

        inds = {obj.u_inds(:); obj.du_inds(:)};
        constraint = FunctionHandleConstraint(repmat(plant.umin,N,1),repmat(plant.umax,N,1),2*obj.nU*N,@obj.delta_u_bound);
        obj = obj.addConstraint(constraint, inds);  
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
      
      num_vars = nH + N*(2*nX+2*nU+nW) + nG + nZ;
      obj.h_inds = (1:nH)';
      obj.x_inds = reshape(nH + (1:nX*N),nX,N);
      obj.u_inds = reshape(nH + nX*N + (1:nU*N),nU,N);
      obj.dx_inds = reshape(nH + (nX+nU)*N +(1:nX*N),nX,N);
      obj.du_inds = reshape(nH + (2*nX+nU)*N +(1:nU*N),nU,N);
      obj.w_inds = reshape(nH + (2*nX+2*nU)*N +(1:nW*N),nW,N);
      obj.gamma_inds = (nH + (2*nX+2*nU+nW)*N + (1:nG))';
      obj.z_inds = (nH + (2*nX+2*nU+nW)*N + nG + (1:nZ))';
      
      x_names = cell(num_vars,1);
      for i = 1:N
        if(i<N)
          x_names{i} = sprintf('h[%d]',i);
          x_names{nH+(2*nX+2*nU+nW)*N+i} = sprintf('g[%d]',i);
          for j = 1:M
            x_names{nH+(2*nX+2*nU+nW)*N+nG+(i-1)*M+j} = sprintf('z%d[%d]',j,i);
          end
        end
        for j = 1:nX
          x_names{nH+(i-1)*nX+j}=sprintf('x%d[%d]',j,i);
          x_names{nH+(nX+nU)*N+(i-1)*nX+j} = sprintf('dx%d[%d]',j,i);
        end
        for j = 1:nU
          x_names{nH+nX*N+(i-1)*nU+j} = sprintf('u%d[%d]',j,i);
          x_names{nH+(2*nX+nU)*N+(i-1)*nU+j} = sprintf('du%d[%d]',j,i);
        end
        for j = 1:nW
          x_names{nH+(2*nX+2*nU)*N+(i-1)*nW+j} = sprintf('w%d[%d]',j,i);
        end
      end

      obj = obj.addDecisionVariable(num_vars,x_names);
    end
    
    function obj = addGammaCost(obj)
      running_cost = FunctionHandleObjective(obj.N-1, @obj.gamma_cost);
      obj = obj.addCost(running_cost,{obj.gamma_inds});
    end
    
    function obj = addZCost(obj)
      for i=1:obj.N-1
        running_cost = FunctionHandleObjective(obj.M, @obj.z_cost);
        inds_i = {obj.z_inds((i-1)*obj.M+(1:obj.M))};
        obj = obj.addCost(running_cost,inds_i);
      end
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
          
          switch obj.options.integration_method
            case DirtranTrajectoryOptimization.FORWARD_EULER
  
              n_vars = 2 + 3*nX + 2*nU;
              inds = {obj.gamma_inds(i); ...
                      obj.h_inds(i); ...
                      obj.x_inds(:,i); ...
                      obj.dx_inds(:,i); ...
                      obj.u_inds(:,i); ...
                      obj.du_inds(:,i); ...
                      obj.x_inds(:,i+1)};
                    
                    
              forward_robust_bound_fun_j = @(gamma,h,x0,dx0,u,du,x1) obj.forward_robust_bound_fun(robust_cost,j,gamma,h,x0,dx0,u,du,x1);

              constraint = FunctionHandleConstraint(0,inf,n_vars,forward_robust_bound_fun_j);
              obj = obj.addConstraint(constraint, inds);

            case DirtranTrajectoryOptimization.BACKWARD_EULER
              error('not implemented yet');

            case DirtranTrajectoryOptimization.MIDPOINT
              error('not implemented yet');
              
            otherwise
              error('Drake:DirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
          end
          
        end
        
        %---------------------------------------------------------------
        %---------------------------------------------------------------
        %---------------------------------------------------------------        
        % complementarity constraint: (gamma-ell)'z=0
        %---------------------------------------------------------------
        %---------------------------------------------------------------
        n_vars = 2 + 3*nX + 2*nU + M;
        inds = {obj.gamma_inds(i); ...
                obj.h_inds(i); ...
                obj.x_inds(:,i); ...
                obj.dx_inds(:,i); ...
                obj.u_inds(:,i); ...
                obj.du_inds(:,i); ...
                obj.x_inds(:,i+1); ...
                obj.z_inds((i-1)*M+(1:M))};
              
        complementarity_fun = @(gamma,h,x0,dx0,u,du,x1,zi) obj.complementarity_fun(robust_cost,gamma,h,x0,dx0,u,du,x1,zi);
              
        constraint = FunctionHandleConstraint(0,0,n_vars,complementarity_fun);
        obj = obj.addConstraint(constraint, inds);
 
%         % Sum constraint on z: \sum_j zij = 1
%         inds = {obj.z_inds((i-1)*M+(1:M))};
%         constraint = FunctionHandleConstraint(1,1,M,@obj.z_sum_constr);
%         obj = obj.addConstraint(constraint, inds);

        % norm constraints on zi
        inds = {obj.z_inds((i-1)*M+(1:M))};
        constraint = FunctionHandleConstraint(1,1,M,@obj.z_norm_constr);
        obj = obj.addConstraint(constraint, inds);
        
        % equality constraint on wi: wi - \sum_j zij = 0
        inds = {obj.w_inds(:,i);obj.z_inds((i-1)*M+(1:M))};
        constraint = FunctionHandleConstraint(zeros(nW,1),zeros(nW,1),nW+M,@obj.w_equality);
        obj = obj.addConstraint(constraint, inds);        
        
        
        
          
        %---------------------------------------------------------------
        %---------------------------------------------------------------
        %---------------------------------------------------------------
        % Add constraint on delta_x for all i=1:N
        %---------------------------------------------------------------
        %---------------------------------------------------------------
        switch obj.options.integration_method
          case DirtranTrajectoryOptimization.FORWARD_EULER
            % add delta_x1 constraint 
            n_vars = 1 + 4*nX + 2*nU + nW;
            inds = {obj.h_inds(i); ...
                    obj.x_inds(:,i); ...
                    obj.dx_inds(:,i); ...
                    obj.x_inds(:,i+1); ...
                    obj.dx_inds(:,i+1); ...
                    obj.u_inds(:,i); ...
                    obj.du_inds(:,i); ...
                    obj.w_inds(:,i)};
            constraint = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.forward_delta_x_constraint);
            obj = obj.addConstraint(constraint, inds);


          case DirtranTrajectoryOptimization.BACKWARD_EULER
            error('not implemented yet');

          case DirtranTrajectoryOptimization.MIDPOINT
            error('not implemented yet');

            % add delta_x1 constraint 
            n_vars = 1 + 4*nX + 3*nU + nW;
            inds = {obj.h_inds(i); ...
                    obj.x_inds(:,i); ...
                    obj.dx_inds(:,i); ...
                    obj.x_inds(:,i+1); ...
                    obj.dx_inds(:,i+1); ...
                    obj.u_inds(:,i); ...
                    obj.du_inds(:,i); ...
                    obj.u_inds(:,i+1); ...
                    obj.w_inds(:,i)};

            constraint = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.midpoint_delta_x_constraint);
            obj = obj.addConstraint(constraint, inds);


          otherwise
            error('Drake:DirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
        end
      end
    
          
      
      %---------------------------------------------------------------
      %---------------------------------------------------------------
      %---------------------------------------------------------------
      %-------- Constrain first dx=0        --------------------------
      %---------------------------------------------------------------
      obj = obj.addConstraint(ConstantConstraint(zeros(nX,1)),obj.dx_inds(1:nX,1)); % first dx=0
      obj = obj.addConstraint(BoundingBoxConstraint(zeros(nZ,1),ones(nZ,1)),obj.z_inds);
%       obj = obj.addConstraint(ConstantConstraint(zeros(nU,1)),obj.u_inds(:,N)); % last u=0
      obj = obj.addConstraint(ConstantConstraint(zeros(nU,1)),obj.du_inds(:,N)); % last du=0
    end
       
    function obj = addDeltaXEqualsZeroConstraint(obj)
      obj = obj.addConstraint(ConstantConstraint(zeros(obj.nX*obj.N,1)),obj.dx_inds(:)); 
    end
  
    function obj = addLinearControlConstraint(obj,K)
      obj.K = K;
      n_vars = obj.nX + obj.nU;
      for i=1:obj.N
        for j=1:obj.M
          inds = {obj.dx_inds((j-1)*obj.nX+(1:obj.nX),i); ...
                  obj.du_inds((j-1)*obj.nU+(1:obj.nU),i)};

          constraint = FunctionHandleConstraint(zeros(obj.nU,1),zeros(obj.nU,1),n_vars,@obj.linear_delta_u_constraint);
          obj = obj.addConstraint(constraint, inds);
        end
      end
    end    
    
    function [f,df] = complementarity_fun(obj,robust_cost,gamma,h,x0,dx0,u,du,x1,zi)
      v = zeros(obj.M,1);
      dv = [];
      for j=1:obj.M
        [g,dg] = forward_robust_bound_fun(obj,robust_cost,j,gamma,h,x0,dx0,u,du,x1);
        v(j) = g;
        dv = [dv;dg];
      end
      f = v'*zi;
      df = [zi'*dv, v'];
    end
    
    function [f,df] = linear_delta_u_constraint(obj,dx,du)
      f = du - obj.K*dx;

      df = [-obj.K, ... dx
            eye(obj.nU)]; % du 
    end
    
    function [f,df] = forward_robust_bound_fun(obj,robust_cost,j,gamma,h,x0,dx0,u,du,x1)
      w = obj.disturbances(:,j);
      
      [x1_w, dx1_w] = forward_robust_dynamics_fun(obj,h,x0+dx0,u,du,w);
     
      dx = x1_w - x1;
      [g,dg] = robust_cost(dx,du,w);
   
      f = gamma - g;
    
      df = [1, ... gamma
            -dg(1:obj.nX)*dx1_w(:,1), ... h
            -dg(1:obj.nX)*dx1_w(:,1+(1:obj.nX)), ... x0
            -dg(1:obj.nX)*dx1_w(:,1+(1:obj.nX)), ... dx0
            -dg(1:obj.nX)*dx1_w(:,1+obj.nX+(1:obj.nU)), ... u
            -dg(1:obj.nX)*dx1_w(:,1+obj.nX+(1:obj.nU)) - dg(obj.nX+(1:obj.nU)), ... du
            dg(1:obj.nX)]; %x1
    end
    
    function [f,df] = forward_delta_x_constraint(obj,h,x0,dx0,x1,dx1,u,du,w)

      [r,dr] = obj.forward_robust_dynamics_fun(h,x0+dx0,u,du,w);
      
      f = dx1 - r + x1;
      
      dr_dh = dr(:,1);
      dr_dx0 = dr(:,1+(1:obj.nX));
      dr_du = dr(:,1+obj.nX+(1:obj.nU));
      dr_ddu = dr(:,1+obj.nX+obj.nU+(1:obj.nU));
      dr_dw = dr(:,1+obj.nX+2*obj.nU+(1:obj.nW));
      
      df_dh = -dr_dh;
      df_dx0 = -dr_dx0;
      df_ddx0 = -dr_dx0;
      df_dx1 = eye(obj.nX);
      df_ddx1 = eye(obj.nX);
      df_du = -dr_du;
      df_ddu = -dr_ddu;
      df_dw = -dr_dw;
      
      df = [df_dh, df_dx0, df_ddx0, df_dx1, df_ddx1, df_du, df_ddu, df_dw];
    end
    
%     function [f,df] = midpoint_delta_x_constraint(obj,h,x0,dx0,x1,dx1j,u0,du0,u1,w)
% 
%       [r,dr] = obj.midpoint_robust_dynamics_fun(h,x0,dx0,x1,u0,du0,u1,w);
%       
%       f = dx1j - r + x1; 
%       
%       dr_dh = dr(:,1);
%       dr_dx0 = dr(:,1+(1:obj.nX));
%       dr_ddx0 = dr(:,1+obj.nX+(1:obj.nX));
%       dr_dx1 = dr(:,1+2*obj.nX+(1:obj.nX));
%       dr_du0 = dr(:,1+3*obj.nX+(1:obj.nU));
%       dr_ddu0 = dr(:,1+3*obj.nX+obj.nU+(1:obj.nU));
%       dr_du1 = dr(:,1+3*obj.nX+2*obj.nU+(1:obj.nU));
%       dr_dw = dr(:,1+3*obj.nX+3*obj.nU+(1:obj.nW));
%       
%       df_dh = -dr_dh;
%       df_dx0 = -dr_dx0;
%       df_ddx0 = -dr_ddx0;
%       df_dx1 = eye(obj.nX) - dr_dx1;
%       df_ddx1 = eye(obj.nX);
%       df_du0 = -dr_du0;
%       df_ddu0 = -dr_ddu0;
%       df_du1 = -dr_du1;
%       df_dw = -dr_dw;
%       
%       df = [df_dh, df_dx0, df_ddx0, df_dx1, df_ddx1, df_du0, df_ddu0, df_du1, df_dw];
%     end
    
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
    
%     function [f,df] = midpoint_robust_dynamics_fun(obj,h,x0,dx0,x1,u0,du0,u1,w)
%       [xdot,dxdot] = obj.plant.dynamics_w(0,.5*(x0+x1)+dx0,.5*(u0+u1)+du0,w);
%       f = x0 + dx0 + h*xdot;
%       df = [xdot ... h
%         (eye(obj.nX) + .5*h*dxdot(:,1+(1:obj.nX))) ... x0
%         (eye(obj.nX) + h*dxdot(:,1+(1:obj.nX))) ... dx0
%         .5*h*dxdot(:,1+(1:obj.nX))... x1
%         .5*h*dxdot(:,1+obj.nX+(1:obj.nU)) ... u0
%         h*dxdot(:,1+obj.nX+(1:obj.nU)) ... du0
%         .5*h*dxdot(:,1+obj.nX+(1:obj.nU)) ... u1
%         h*dxdot(:,1+obj.nX+obj.nU+(1:obj.nW))]; % w
%     end
    
    function [f,df] = gamma_cost(obj,gamma)
      f = ones(1,obj.N-1)*gamma; 
      df = ones(1,obj.N-1);
    end

    function [f,df] = delta_u_bound(obj,u,du)
      f = u+du;
      df = [eye(length(u)), eye(length(u))];
    end
    
    function [f,df] = w_equality(obj,wi,zi)
      f = wi - obj.disturbances*zi;
      df = [eye(obj.nW), -obj.disturbances];
    end

%     function [f,df] = z_diff_constr(obj,z)
%       mf = factorial(obj.M-1);
%       f = zeros(mf,1);
%       df = zeros(mf,obj.M);
%       
%       n=1;
%       for i=1:obj.M
%         for j=i+1:obj.M
%           err = z(i)-z(j);
%           f(n) = err'*err;
%           df(n,[i,j]) = [2*err,-2*err];
%           n=n+1;
%         end
%       end
%     end
    
    function [f,df] = z_sum_constr(~,zi)
      f = ones(1,length(zi))*zi;
      df = ones(1,length(zi));
    end
    
    function [f,df] = z_norm_constr(~,zi)
      f = zi'*zi;
      df = 2*zi';
    end
    
%     function [f,df] = z_cost(obj,zi)      
%       f = 0;
%       df = zeros(1,obj.M);
%       for j=1:obj.M
%         e = zeros(obj.M,1); e(j) = 1;
%         err=zi-e;
%         [c,dc] = smin(err'*err,0.5,0.1);
%         f = f+c;
%         df = df + dc*2*err';
%       end
%       f = f*500;
%       df = df*500;
%     end
  end
end