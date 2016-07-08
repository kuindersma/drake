classdef RobustDirtranTrajectoryOptimization < DirtranTrajectoryOptimization
    
  properties
    dx_inds  % n x N indices for delta state in disturbed trajectory
    du_inds  % m x N indices for delta inputs in disturbed trajectory
    w_inds  % d x N indices for disturbances
    D % Disturbance ellipsoid matrix
    nX
    nU
    nW
    Q
    R
    Qf
  end
  
  methods
    function obj = RobustDirtranTrajectoryOptimization(plant,N,D,Q,R,Qf,duration,options)
      if nargin < 5
        options = struct();
      end
      if isscalar(duration), duration=[duration,duration]; end

      if ~isfield(options,'integration_method')
        options.integration_method = DirtranTrajectoryOptimization.MIDPOINT;
      end
      obj = obj@DirtranTrajectoryOptimization(plant,N,duration,options);
      obj.D = D;
      obj.Q = Q;
      obj.R = R;
      obj.Qf = Qf;

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
    
    function obj = addRobustCost(obj,robust_cost)
        nx = obj.nX;
        nu = obj.nU;
        nw = obj.nW;
        N = obj.N;
        
        dim = nx+nu+nw;
        for k = 1:N
            costObj = FunctionHandleObjective(dim,robust_cost,2);
            obj = obj.addCost(costObj, {obj.dx_inds(:,k); obj.du_inds(:,k); obj.w_inds(:,k)});
        end
    end
    
    function obj = addRobustConstraints(obj,robust_cost)
        nx = obj.nX;
        nu = obj.nU;
        nw = obj.nW;
        N = obj.N;
        
        switch obj.options.integration_method
            case DirtranTrajectoryOptimization.FORWARD_EULER
                for i=1:N-1
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
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Equality constraint on delta-u
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                lqrconst = FunctionHandleConstraint(zeros(N*nu,1), zeros(N*nu,1), obj.num_vars, @obj.lqr_constraint);
                lqrconst.grad_level = 0; %need to add derivatives
                lqrconst.grad_method = 'numerical';
                obj = obj.addConstraint(lqrconst);
                
            case DirtranTrajectoryOptimization.MIDPOINT
                for i=1:N-1
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
                    
                    constraint = FunctionHandleConstraint(zeros(nx,1),zeros(nx,1),n_vars,@obj.midpoint_delta_x_constraint);
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
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Equality constraint on delta-u
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                lqrconst = FunctionHandleConstraint(zeros(N*nu,1), zeros(N*nu,1), obj.num_vars, @obj.lqr_constraint);
                lqrconst.grad_level = 0; %need to add derivatives
                lqrconst.grad_method = 'numerical';
                obj = obj.addConstraint(lqrconst);
                
            case DirtranTrajectoryOptimization.BACKWARD_EULER
                error('not implemented yet');
                
            otherwise
                error('Drake:DirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
        end
      
      %Constrain first dx=0
      obj = obj.addConstraint(ConstantConstraint(zeros(nx,1)),obj.dx_inds(:,1)); % first dx=0
    end
       
    function obj = addDeltaXEqualsZeroConstraint(obj)
      constraint = ConstantConstraint(zeros(obj.nX*obj.N,1));
      constraint = constraint.setName('dx_eq_0');      
      obj = obj.addConstraint(constraint,obj.dx_inds(:)); 
    end
    
    function f = lqr_constraint(obj,z,sp)
        nx = obj.nX;
        nu = obj.nU;
        N = obj.N;
        
        %Get linearized dynamics along trajectory + derivatives
        Ak = zeros(nx,nx,N-1);
        Bk = zeros(nx,nu,N-1);
        switch obj.options.integration_method
            case DirtranTrajectoryOptimization.FORWARD_EULER
                for k = 1:(N-1)
                    [~,dx1] = obj.forward_robust_dynamics_fun(z(obj.h_inds(k)),z(obj.x_inds(:,k)),z(obj.u_inds(:,k)), zeros(1,obj.nW));
                    Ak(:,:,k) = dx1(:,1+(1:nx));
                    Bk(:,:,k) = dx1(:,1+nx+(1:nu));
                end
            case DirtranTrajectoryOptimization.MIDPOINT
                for k = 1:(N-1)
                    [~,dx1] = obj.forward_robust_dynamics_fun(z(obj.h_inds(k)),.5*(z(obj.x_inds(:,k))+z(obj.x_inds(:,k+1))),z(obj.u_inds(:,k)), zeros(1,obj.nW));
                    Ak(:,:,k) = dx1(:,1+(1:nx));
                    Bk(:,:,k) = dx1(:,1+nx+(1:nu));
                end
        end
        
        %Solve Riccati Equation
        S = obj.Qf;
        dx = z(obj.dx_inds(:));
        du = zeros(N*nu,1);
        for k = (N-1):-1:1
            K = (Bk(:,:,k)'*S*Bk(:,:,k)+obj.R)\(Bk(:,:,k)'*S*Ak(:,:,k));
            du(k) = -K*dx((k-1)*nx+(1:nx));
            S = obj.Q + K'*obj.R*K + (Ak(:,:,k) - Bk(:,:,k)*K)'*S*(Ak(:,:,k) - Bk(:,:,k)*K);
        end
        
        f = z(obj.du_inds(:)) - du;
        
    end
   
    function [f,df] = forward_w_equality_constraint(obj,robust_cost,h,x0,dx0,u,du,w)
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
        
        persistent wstar; %keep last value for warm starting
        
        %Setup QCQP
        [~,dx1,ddx1] = obj.forward_robust_dynamics_fun(h,x0,u0,zeros(obj.nW,1));
        [~,dJ,ddJ] = robust_cost(zeros(obj.nX,1),zeros(obj.nU,1),zeros(obj.nW,1));
        
        %Dynamics derivatives
        G = dx1(:,1+obj.nX+obj.nU+(1:obj.nW));
        T = reshape(full(ddx1),obj.nX,1+obj.nX+obj.nU+obj.nW,1+obj.nX+obj.nU+obj.nW); %2nd derivative
        dG = T(:,1+obj.nX+obj.nU+(1:obj.nW),:);
        
        %Cost derivatives for QP
        H = -G'*ddJ(1:obj.nX,1:obj.nX)*G;
        f = -G'*dJ(1:obj.nX)';
        
        %wstar = qcqp_mex(H,f,obj.D,wstar);
        w = qcqp(H,f,obj.D,wstar);
        wstar = w;
        
        %Evaluate derivatives
        dw = -H\tvMult(dG,2*ddJ(1:obj.nX,1:obj.nX)*G*w + dJ(1:obj.nX)',1);
        if abs(.5*w'*obj.D*w - 1) < .01
            %Project onto the constraint surface
            n = obj.D*w;
            n = n/norm(n);
            P = (eye(obj.nW)-n*n');
            dw = P*dw;
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
    
    function [f,df,d2f] = forward_robust_dynamics_fun(obj,h,x,u,w)
      if nargout == 1
          xdot = obj.plant.dynamics_w(0,x,u,w);
          f = x + h*xdot;
      elseif nargout == 2
        [xdot,dxdot] = obj.plant.dynamics_w(0,x,u,w);
        f = x + h*xdot;
        df = [xdot ... h
          eye(obj.nX) + h*dxdot(:,1+(1:obj.nX)) ... x0
          h*dxdot(:,1+obj.nX+(1:obj.nU)) ... u
          h*dxdot(:,1+obj.nX+obj.nU+(1:obj.nW))]; % w
      else %nargout == 3
          [xdot,dxdot,d2xdot] = obj.plant.dynamics_w(0,x,u,w);
          f = x + h*xdot;
          df = [xdot ... h
            eye(obj.nX) + h*dxdot(:,1+(1:obj.nX)) ... x0
            h*dxdot(:,1+obj.nX+(1:obj.nU)) ... u
            h*dxdot(:,1+obj.nX+obj.nU+(1:obj.nW))]; % w
          d2f = h*d2xdot;
          d2f(:,1:(1+obj.nX+obj.nU+obj.nW)) = dxdot;
          d2f(:,1:(1+obj.nX+obj.nU+obj.nW):end) = dxdot;
      end
    end

   
    function [f,df] = midpoint_delta_x_constraint(obj,h,x0,dx0,x1,dx1,u0,du0,w)

      [r,dr] = obj.midpoint_robust_dynamics_fun(h,x0,dx0,x1,dx1,u0,du0,w);
      
      f = dx1 - r + x1; 
      
      dr_dh = dr(:,1);
      dr_dx0 = dr(:,1+(1:obj.nX));
      dr_ddx0 = dr(:,1+obj.nX+(1:obj.nX));
      dr_dx1 = dr(:,1+2*obj.nX+(1:obj.nX));
      dr_ddx1 = dr(:,1+3*obj.nX+(1:obj.nX));
      dr_du0 = dr(:,1+4*obj.nX+(1:obj.nU));
      dr_ddu0 = dr(:,1+4*obj.nX+obj.nU+(1:obj.nU));
      dr_dw = dr(:,1+4*obj.nX+2*obj.nU+(1:obj.nW));
      
      df_dh = -dr_dh;
      df_dx0 = -dr_dx0;
      df_ddx0 = -dr_ddx0;
      df_dx1 = eye(obj.nX) - dr_dx1;
      df_ddx1 = eye(obj.nX) - dr_ddx1;
      df_du0 = -dr_du0;
      df_ddu0 = -dr_ddu0;
      df_dw = -dr_dw;
      
      df = [df_dh, df_dx0, df_ddx0, df_dx1, df_ddx1, df_du0, df_ddu0, df_dw];
    end

    function [f,df] = midpoint_robust_dynamics_fun(obj,h,x0,dx0,x1,dx1,u0,du0,w)
      %2nd order midpoint on dynamics, FOH on control and disturbance inputs
      if nargout == 1
          xdot = obj.plant.dynamics_w(0,.5*(x0+dx0+x1+dx1),u0+du0,w);
          f = x0 + dx0 + h*xdot;
      elseif nargout == 2
          [xdot,dxdot] = obj.plant.dynamics_w(0,.5*(x0+dx0+x1+dx1),u0+du0,w);
          f = x0 + dx0 + h*xdot;
          df = [xdot ... h
              (eye(obj.nX) + .5*h*dxdot(:,1+(1:obj.nX))) ... x0
              (eye(obj.nX) + .5*h*dxdot(:,1+(1:obj.nX))) ... dx0
              .5*h*dxdot(:,1+(1:obj.nX))... x1
              .5*h*dxdot(:,1+(1:obj.nX))... dx1
              h*dxdot(:,1+obj.nX+(1:obj.nU)) ... u0
              h*dxdot(:,1+obj.nX+(1:obj.nU)) ... du0
              h*dxdot(:,1+obj.nX+obj.nU+(1:obj.nW))]; % w
%       else %nargout == 3
%           [xdot,dxdot,d2xdot] = obj.plant.dynamics_w(0,.5*(x0+dx0+x1+dx1),u0+du0,w);
%           
%           f = x0 + dx0 + h*xdot;
%           
%           df = [xdot ... h
%               (eye(obj.nX) + .5*h*dxdot(:,1+(1:obj.nX))) ... x0
%               (eye(obj.nX) + .5*h*dxdot(:,1+(1:obj.nX))) ... dx0
%               .5*h*dxdot(:,1+(1:obj.nX))... x1
%               .5*h*dxdot(:,1+(1:obj.nX))... dx1
%               h*dxdot(:,1+obj.nX+(1:obj.nU)) ... u0
%               h*dxdot(:,1+obj.nX+(1:obj.nU)) ... du0
%               h*dxdot(:,1+obj.nX+obj.nU+(1:obj.nW))]; % w
%           
%           nVars = 1+4*obj.nX+2*obj.nU+obj.nW;
%           dxdotdx = dxdot(:,1+(1:obj.nX));
%           dxdotdu = dxdot(:,1+obj.nX+(1:obj.nU));
%           dxdotdw = dxdot(:,1+obj.nX+obj.nU+(1:obj.nW));
%           d2f = zeros(obj.nX,(1+4*obj.nX+2*obj.nU+obj.nW).^2);
%           d2f(:,1:nVars) = [zeros(obj.nX,1), ... h
%                             .5*dxdotdx, ... x0
%                             .5*dxdotdx, ... dx0
%                             .5*dxdotdx, ... x1
%                             .5*dxdotdx, ... dx1
%                             dxdotdu, ... u0
%                             dxdotdu, ... du0
%                             dxdotdw]; % w
      end
      
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

  end
end