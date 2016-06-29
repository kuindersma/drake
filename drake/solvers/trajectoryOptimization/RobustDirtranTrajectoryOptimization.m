classdef RobustDirtranTrajectoryOptimization < DirtranTrajectoryOptimization
  properties
    dx_inds  % n x N indices for delta state in disturbed trajectory
    du_inds  % m x N indices for delta inputs in disturbed trajectory
    w_inds  % d x N indices for disturbances
    D % Disturbance ellipsoid matrix
    nX
    nU
    nW
    K
  end
  
  methods
    function obj = RobustDirtranTrajectoryOptimization(plant,N,D,duration,options)
      if nargin < 5
        options = struct();
      end
      if isscalar(duration), duration=[duration,duration]; end

      if ~isfield(options,'integration_method')
        options.integration_method = DirtranTrajectoryOptimization.FORWARD_EULER;
      end
      obj = obj@DirtranTrajectoryOptimization(plant,N,duration,options);
      obj.D = D;

      % add control inputs constraints
      if any(~isinf(plant.umin)) || any(~isinf(plant.umax))
        % u+du limit constraint
        inds = {obj.u_inds(:); obj.du_inds(:)};
        n = obj.nU*N;
        idx_u = zeros(n,2*n); idx_u(:,1:n) = eye(n);
        idx_du = zeros(n,2*n); idx_du(:,n+(1:n)) = eye(n);
        A = sparse(idx_u + idx_du); 
        constraint = LinearConstraint(repmat(plant.umin,N,1),repmat(plant.umax,N,1),A);
        constraint = constraint.setName('u+du limit');
        obj = obj.addConstraint(constraint, inds);
      end        
    end
    
    function obj = setupVariables(obj, N)
      nH = N-1;
      nx = obj.plant.getNumStates();
      nu = obj.plant.getNumInputs();
      nw = obj.plant.getNumDisturbances();
      
      num_vars = nH + N*(2*nx+2*nu+nw);
      obj.h_inds = (1:nH)';
      obj.x_inds = reshape(nH + (1:nx*N),nx,N);
      obj.u_inds = reshape(nH + nx*N + (1:nu*N),nu,N);
      obj.dx_inds = reshape(nH + (nx+nu)*N +(1:nx*N),nx,N);
      obj.du_inds = reshape(nH + (2*nx+nu)*N +(1:nu*N),nu,N);
      obj.w_inds = reshape(nH + (2*nx+2*nu)*N +(1:nw*N),nw,N);
      
      x_names = cell(num_vars,1);
      for i = 1:N
        if(i<N)
          x_names{i} = sprintf('h[%d]',i);
        end
        for j = 1:nx
          x_names{nH+(i-1)*nx+j}=sprintf('x%d[%d]',j,i);
          x_names{nH+(nx+nu)*N+(i-1)*nx+j} = sprintf('dx%d[%d]',j,i);
        end
        for j = 1:nu
          x_names{nH+nx*N+(i-1)*nu+j} = sprintf('u%d[%d]',j,i);
          x_names{nH+(2*nx+nu)*N+(i-1)*nu+j} = sprintf('du%d[%d]',j,i);
        end
        for j = 1:nw
          x_names{nH+(2*nx+2*nu)*N+(i-1)*nw+j} = sprintf('w%d[%d]',j,i);
        end
      end

      obj = obj.addDecisionVariable(num_vars,x_names);

      obj.nX = nx;
      obj.nU = nu;
      obj.nW = nw;
      obj.N = N;
    end
    
    function obj = addRobustConstraints(obj,robust_cost)
      nx = obj.nX;
      nu = obj.nU;
      nw = obj.nW;
      N = obj.N;
      
      for i=1:N-1
        
        switch obj.options.integration_method
          case DirtranTrajectoryOptimization.FORWARD_EULER
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Add constraint on delta_x for all i=1:N
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            n_vars = 1 + 4*nx + 2*nu + nw;
            inds = {obj.h_inds(i); ...
                    obj.x_inds(:,i); ...
                    obj.dx_inds(:,i); ...
                    obj.x_inds(:,i+1); ...
                    obj.dx_inds(:,i+1); ...
                    obj.u_inds(:,i); ...
                    obj.du_inds(:,i); ...
                    obj.w_inds(:,i)};
            constraint = FunctionHandleConstraint(zeros(nx,1),zeros(nx,1),n_vars,@obj.forward_delta_x_constraint);
            constraint = constraint.setName(sprintf('delta_x_%d',i));
            obj = obj.addConstraint(constraint, inds);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Equality constraint on w_i
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            n_vars = 1 + 2*nx + 2*nu + nw;
            inds = {obj.h_inds(i); ...
                    obj.x_inds(:,i); ...
                    obj.dx_inds(:,i); ...
                    obj.u_inds(:,i); ...
                    obj.du_inds(:,i); ...
                    obj.w_inds(:,i)};

            forward_w_equality_constraint_i = @(h,x0,dx0,u,du,w0) obj.forward_w_equality_constraint(robust_cost,h,x0,dx0,u,du,w0);

            constraint = FunctionHandleConstraint(zeros(nw,1),zeros(nw,1),n_vars,forward_w_equality_constraint_i);
            constraint = constraint.setName(sprintf('w_equality_%d',i));
            obj = obj.addConstraint(constraint, inds);      
            

          case DirtranTrajectoryOptimization.MIDPOINT
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % Add constraint on delta_x for all i=1:N
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             n_vars = 1 + 4*nx + 3*n + nw;
%             inds = {obj.h_inds(i); ...
%                     obj.x_inds(:,i); ...
%                     obj.dx_inds(:,i); ...
%                     obj.x_inds(:,i+1); ...
%                     obj.dx_inds(:,i+1); ...
%                     obj.u_inds(:,i); ... 
%                     obj.du_inds(:,i); ...
%                     obj.u_inds(:,i+1); ...
%                     obj.w_inds(:,i)};
% 
%             constraint = FunctionHandleConstraint(zeros(nx,1),zeros(nx,1),n_vars,@obj.midpoint_delta_x_constraint);
%             obj = obj.addConstraint(constraint, inds);
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % Bound on robust cost of w_i
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             n_vars = 2 + 3*nX + 3*nU + nW;
%             inds = {obj.gamma_inds(i); ...
%                     obj.h_inds(i); ...
%                     obj.x_inds(:,i); ...
%                     obj.dx_inds(:,i); ...
%                     obj.x_inds(:,i+1); ...
%                     obj.u_inds(:,i); ...
%                     obj.du_inds(:,i); ...
%                     obj.u_inds(:,i+1); ...
%                     obj.w_inds(:,i)};
% 
%             midpoint_w_bound_fun_i = @(gamma,h,x0,dx0,x1,u0,du0,u1,w0) obj.midpoint_w_bound_fun(robust_cost,gamma,h,x0,dx0,x1,u0,du0,u1,w0);
%             constraint = FunctionHandleConstraint(0,inf,n_vars,midpoint_w_bound_fun_i);
%             obj = obj.addConstraint(constraint, inds);            
            
          case DirtranTrajectoryOptimization.BACKWARD_EULER
            error('not implemented yet');
          
          otherwise
            error('Drake:DirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
        end
      end
      
          
      %---------------------------------------------------------------
      %---------------------------------------------------------------
      %-------- Constrain first dx=0        --------------------------
      %---------------------------------------------------------------
      obj = obj.addConstraint(ConstantConstraint(zeros(nx,1)),obj.dx_inds(:,1)); % first dx=0
      %obj = obj.addConstraint(ConstantConstraint(zeros(nu,1)),obj.u_inds(:,N)); % last u=0
      %obj = obj.addConstraint(ConstantConstraint(zeros(nu,1)),obj.du_inds(:,N)); % last du=0
    end
       
    function obj = addDeltaXEqualsZeroConstraint(obj)
      constraint = ConstantConstraint(zeros(obj.nX*obj.N,1));
      constraint = constraint.setName('dx_eq_0');      
      obj = obj.addConstraint(constraint,obj.dx_inds(:)); 
    end
  
    function obj = addLinearControlConstraint(obj,K)
      obj.K = K;
      for i=1:obj.N
        inds = {obj.dx_inds(:,i); ...
                obj.du_inds(:,i)};

        nx = obj.nX;
        nu = obj.nU;
        idx_dx = zeros(nx,nx+nu); idx_dx(:,1:nx) = eye(nx);
        idx_du = zeros(nu,nx+nu); idx_du(:,nx+(1:nu)) = eye(nu);
        A = idx_du - K*idx_dx; 
        constraint = LinearConstraint(zeros(nu,1),zeros(nu,1),A);
        constraint = constraint.setName(sprintf('linear_control_%d',i));
        obj = obj.addConstraint(constraint, inds);  
      end
    end  
   
    function [f,df] = forward_w_equality_constraint(obj,robust_cost,h,x0,dx0,u,du,w)
      option.grad_method = 'numerical';
      %[w_opt, dw] = geval(@(hq,xq,uq)solve_qcqp(obj, robust_cost, hq, xq, uq), h, x0+dx0, u+du, option);
      [w_opt, dw] = solve_qcqp(obj, robust_cost, h, x0+dx0, u+du);
     
      f = w - w_opt;

      df = [-dw(:,1), ... h
            -dw(:,1+(1:obj.nX)), ... x0
            -dw(:,1+(1:obj.nX)), ... dx0
            -dw(:,1+obj.nX+(1:obj.nU)), ... u
            -dw(:,1+obj.nX+(1:obj.nU)), ... % du
            eye(obj.nW)]; 
    end
    
    function [w, dw] = solve_qcqp(obj,robust_cost,h,x0,u0)
        %Setup QCQP
        
        [~,dx1,ddx1] = geval(@obj.forward_robust_dynamics_fun,h,x0,u0,zeros(obj.nW,1));
        [~,dJ,ddJ] = robust_cost(zeros(obj.nX,1),zeros(obj.nU,1),zeros(obj.nW,1));
        
        %Dynamics derivatives
        G = dx1(:,1+obj.nX+obj.nU+(1:obj.nW));
        T = reshape(full(ddx1),obj.nX,1+obj.nX+obj.nU+obj.nW,1+obj.nX+obj.nU+obj.nW); %3rd derivative
        dG = T(:,1+obj.nX+obj.nU+(1:obj.nW),:);
        
        %Cost derivatives for QP
        H = -G'*ddJ(1:obj.nX,1:obj.nX)*G;
        f = -G'*dJ(1:obj.nX)';
        
        %[w,lambda] = qcqp_mex(H,f,obj.D);
        [w,lambda] = qcqp(H,f,obj.D);
        
        %Evaluate derivatives
        dw = -H\tvMult(dG,2*ddJ(1:obj.nX,1:obj.nX)*G*w + dJ(1:obj.nX)',1);
        if lambda < 1e-4
            %Project onto the constraint surface
            n = obj.D*w;
            n = n/norm(n);
            P = (eye(obj.nW)-n*n');
            dw = P*dw;
        end   
        
        function y = tvMult(T,x,ind)
            if ind == 1
                y = zeros(size(T,2),size(T,3));
                for j = 1:size(T,2)
                    for k = 1:size(T,3)
                        y(j,k) = T(:,j,k)'*x;
                    end
                end
            elseif ind == 2
                y = zeros(size(T,1),size(T,3));
                for j = 1:size(T,1)
                    for k = 1:size(T,3)
                        y(j,k) = T(j,:,k)*x;
                    end
                end
            else %ind == 3
                y = zeros(size(T,1),size(T,2));
                for j = 1:size(T,1)
                    for k = 1:size(T,2)
                        y(j,k) = squeeze(T(j,k,:))'*x;
                    end
                end
            end
        end
        
    end
    
    function [f,df] = forward_delta_x_constraint(obj,h,x0,dx0,x1,dx1,u,du,w)

      [r,dr] = obj.forward_robust_dynamics_fun(h,x0+dx0,u+du,w);
      
      f = x1 + dx1 - r;
      
      dr_dh = dr(:,1);
      dr_dx0 = dr(:,1+(1:obj.nX));
      dr_du = dr(:,1+obj.nX+(1:obj.nU));
      dr_dw = dr(:,1+obj.nX+obj.nU+(1:obj.nW));
      
      df_dh = -dr_dh;
      df_dx0 = -dr_dx0;
      df_ddx0 = -dr_dx0;
      df_dx1 = eye(obj.nX);
      df_ddx1 = eye(obj.nX);
      df_du = -dr_du;
      df_ddu = -dr_du;
      df_dw = -dr_dw;
      
      df = [df_dh, df_dx0, df_ddx0, df_dx1, df_ddx1, df_du, df_ddu, df_dw];
    end
    
    function [f,df] = forward_robust_dynamics_fun(obj,h,x,u,w)
      [xdot,dxdot] = obj.plant.dynamics_w(x,u,w);
      f = x + h*xdot;
      df = [xdot ... h
        eye(obj.nX) + h*dxdot(:,1+(1:obj.nX)) ... x0
        h*dxdot(:,1+obj.nX+(1:obj.nU)) ... u
        h*dxdot(:,1+obj.nX+obj.nU+(1:obj.nW))]; % w
    end
    
%      
%       function [f,df] = midpoint_w_bound_fun(obj,robust_cost,gamma,h,x0,dx0,x1,u0,du0,u1,w0)
%       [r, dr] = midpoint_robust_dynamics_fun(obj,h,x0,dx0,x1,u0,du0,u1,w0);
%      
%       dx = r - x1;
%       [g,dg] = robust_cost(dx,du0,w0);
%    
%       f = g - gamma;
%     
%       df = [-1, ... gamma
%             dg(1:obj.nX)*dr(:,1), ... h
%             dg(1:obj.nX)*dr(:,1+(1:obj.nX)), ... x0
%             dg(1:obj.nX)*dr(:,1+obj.nX+(1:obj.nX)), ... dx0
%             dg(1:obj.nX)*(dr(:,1+2*obj.nX+(1:obj.nX)) - eye(obj.nX)), ... x1
%             dg(1:obj.nX)*dr(:,1+3*obj.nX+(1:obj.nU)), ... u0
%             dg(1:obj.nX)*dr(:,1+3*obj.nX+obj.nU+(1:obj.nU)) + dg(obj.nX+(1:obj.nU)), ... du0
%             dg(1:obj.nX)*dr(:,1+3*obj.nX+2*obj.nU+(1:obj.nU)), ... u1
%             dg(1:obj.nX)*dr(:,1+3*obj.nX+3*obj.nU+(1:obj.nW)) + dg(obj.nX+obj.nU+(1:obj.nW))]; % w0
%     end
%     
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
%     
%     
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
%     

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
    
%     function [f,df] = forward_w_bound_fun(obj,robust_cost,gamma,dx1,du0,w0)
%       [g,dg] = robust_cost(dx1,du0,w0);
%    
%       f = g - gamma;
%     
%       df = [-1, ... gamma
%             dg(1:obj.nX), ... dx1
%             dg(obj.nX+(1:obj.nU)), ... du0
%             dg(obj.nX+obj.nU+(1:obj.nW))]; %w0
%     end


%     function [f,df] = delta_u_w_constr(~,du,w)
%       f = du+w;
%       df = [eye(length(du)), eye(length(w))];
%     end



    

  end
end