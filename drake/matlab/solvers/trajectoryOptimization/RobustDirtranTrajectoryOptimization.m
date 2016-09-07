 classdef RobustDirtranTrajectoryOptimization < DirectTrajectoryOptimization
    %  For forward euler integratino:
    %    dynamics constraints are: x(k+1) = x(k) + h(k)*f(x(k),u(k))
    %    integrated cost is sum of g(h(k),x(k),u(k))
    %  For midpoint integration:
    %    dynamics constraints are: x(k+1) = x(k) + h(k)*f(.5*x(k)+.5*x(k+1),u(k))
    %    integrated cost is sum of g(h(k),.5*x(k)+.5*x(k+1),u(k))
    
  properties (Constant)
    FORWARD_EULER = 1;
    MIDPOINT = 3;  % DEFAULT
  end
  
  properties
    nX
    nU
    nW
    
    D % Disturbance ellipsoid matrix w'*D*w <= 1
    Dinv
    L
    
    Q % LQR state cost matrix
    R % LQR input cost matrix
    Qf% LQR terminal cost matrix
    
    Qr %Robust cost matrix
    Rr %Robust cost matrix
    Qrf %Robust cost matrix
  end
  
  methods
    function obj = RobustDirtranTrajectoryOptimization(plant,N,D,Q,R,Qf,duration,options)
      
      if nargin < 8
        options = struct();
      end
      if ~isfield(options,'integration_method')
        options.integration_method = RobustDirtranTrajectoryOptimization.MIDPOINT;
      end
      if isscalar(duration)
          duration=[duration,duration];
      end

      obj = obj@DirectTrajectoryOptimization(plant,N,duration,options);
      
      obj.nX = plant.getNumStates;
      obj.nU = plant.getNumInputs;
      obj.nW = plant.getNumDisturbances;
      obj.D = D;
      obj.Dinv = inv(D);
      obj.L = chol(obj.Dinv,'lower');
      obj.Q = Q;
      obj.R = R;
      obj.Qf = Qf;
      
    end
    
    function obj = setupVariables(obj, N)
      % Assumes that there are N-1 time steps
      % N corresponding state variables
      % and N-1 corresponding input variables
      %
      % Generates num_vars total number of decision variables
      %   h_inds (N-1) x 1 indices for timesteps h so that z(h_inds(i)) = h(i)
      %   x_inds N x n indices for state
      %   u_inds (N-1) x m indices for input
      %
      % @param N number of knot points
      nH = N-1;
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();

      num_vars = nH + N*nX + (N-1)*nU;
      obj.h_inds = (1:nH)';
      obj.x_inds = reshape(nH + (1:nX*N),nX,N);
      obj.u_inds = reshape(nH + nX*N + (1:nU*(N-1)),nU,N-1);

      obj.N = N;
      x_names = cell(num_vars,1);
      for i = 1:(N-1)
        x_names{i} = sprintf('h[%d]',i);
        for j = 1:nX
          x_names{nH+(i-1)*nX+j}=sprintf('x%d[%d]',j,i);
        end
        for j = 1:nU
          x_names{nH+nX*N+(i-1)*nU+j} = sprintf('u%d[%d]',j,i);
        end
      end
      for j = 1:nX
          x_names{nH+(N-1)*nX+j}=sprintf('x%d[%d]',j,N);
      end

      obj = obj.addDecisionVariable(num_vars,x_names);
      
      % Ensure that all h values are non-negative
      obj = obj.addConstraint(BoundingBoxConstraint(zeros(N-1,1),inf(N-1,1)),obj.h_inds);
      
      % create constraints for dynamics and add them
      obj = obj.addDynamicConstraints();
      
      % add control inputs as bounding box constraints
      if any(~isinf(obj.plant.umin)) || any(~isinf(obj.plant.umax))
          control_limit = BoundingBoxConstraint(repmat(obj.plant.umin,N-1,1),repmat(obj.plant.umax,N-1,1));
          obj = obj.addConstraint(control_limit,obj.u_inds(:));
      end
      
    end
    
    function z0 = getInitialVars(obj,t_init,traj_init)
        % evaluates the initial trajectories at the sampled times and
        % constructs the nominal z0.
        if isscalar(t_init)
            t_init = linspace(0,t_init,obj.N);
        elseif length(t_init) ~= obj.N
            error('The initial sample times must have the same length as property N')
        end
        z0 = zeros(obj.num_vars,1);
        z0(obj.h_inds) = diff(t_init);
        
        if nargin<3, traj_init = struct(); end
        
        nU = getNumInputs(obj.plant);
        if isfield(traj_init,'u')
            z0(obj.u_inds) = traj_init.u.eval(t_init);
        else
            z0(obj.u_inds) = 0.01*randn(nU,obj.N-1);
        end
        
        if isfield(traj_init,'x')
            z0(obj.x_inds) = traj_init.x.eval(t_init);
        else
            if nU>0
                if ~isfield(traj_init,'u')
                    traj_init.u = setOutputFrame(PPTrajectory(foh(t_init,reshape(z0(obj.u_inds),nU,obj.N-1))),getInputFrame(obj.plant));
                end
                
                % todo: if x0 and xf are equality constrained, then initialize with
                % a straight line from x0 to xf (this was the previous behavior)
                
                %simulate
                sys_ol = cascade(traj_init.u,obj.plant);
            else
                sys_ol = obj.plant;
            end
            
            if ~isfield(traj_init,'x0')
                [~,x_sim] = sys_ol.simulate([t_init(1) t_init(end)]);
            else
                [~,x_sim] = sys_ol.simulate([t_init(1) t_init(end)],traj_init.x0);
            end
            
            z0(obj.x_inds) = x_sim.eval(t_init);
        end
    end
    
    function obj = addDynamicConstraints(obj)
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
      N = obj.N;
      
      constraints = cell(N-1,1);
      dyn_inds = cell(N-1,1);
      
      switch obj.options.integration_method
        case RobustDirtranTrajectoryOptimization.FORWARD_EULER
          n_vars = 2*nX + nU + 1;
          cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.forward_constraint_fun);
        case RobustDirtranTrajectoryOptimization.MIDPOINT
          n_vars = 2*nX + nU + 1;
          cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.midpoint_constraint_fun);
        otherwise
          error('Drake:RobustDirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
      end
      
      for i=1:obj.N-1,
        switch obj.options.integration_method
          case RobustDirtranTrajectoryOptimization.FORWARD_EULER
            dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i)};
          case RobustDirtranTrajectoryOptimization.MIDPOINT
            dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i)};
          otherwise
            error('Drake:RobustDirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
        end
        cnstr = cnstr.setName(sprintf('dynamics_constr_%d_',i));
        constraints{i} = cnstr;
        
        obj = obj.addConstraint(constraints{i}, dyn_inds{i});
      end
    end
    
    function obj = addRunningCost(obj,running_cost_function,grad_level)
      % Adds an integrated cost to all time steps, which is
      % numerical implementation specific (thus abstract)
      % this cost is assumed to be time-invariant
      % @param running_cost_function a function handle
      %  of the form running_cost_function(dt,x,u)
      
      if nargin < 3
        grad_level = -1;
      end
      
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
      
      for i=1:obj.N-1,
        switch obj.options.integration_method
          case DirtranTrajectoryOptimization.FORWARD_EULER
            running_cost = FunctionHandleObjective(1+nX+nU, running_cost_function,grad_level);
            inds_i = {obj.h_inds(i);obj.x_inds(:,i);obj.u_inds(:,i)};
          case DirtranTrajectoryOptimization.BACKWARD_EULER
            running_cost = FunctionHandleObjective(1+nX+nU, running_cost_function,grad_level);
            inds_i = {obj.h_inds(i);obj.x_inds(:,i+1);obj.u_inds(:,i)};
          case DirtranTrajectoryOptimization.MIDPOINT
            running_cost = FunctionHandleObjective(1+2*nX+nU,...
              @(h,x0,x1,u0) obj.midpoint_running_fun(running_cost_function,h,x0,x1,u0),grad_level);
            inds_i = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i)};
          case DirtranTrajectoryOptimization.DT_SYSTEM
            running_cost = FunctionHandleObjective(1+nX+nU, running_cost_function,grad_level);
            inds_i = {obj.h_inds(i);obj.x_inds(:,i);obj.u_inds(:,i)};
          otherwise
            error('Drake:DirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
        end
        
        obj = obj.addCost(running_cost,inds_i);
      end
    end
    
    function obj = addRobustCost(obj,Qr,Rr,Qrf)
        nX = obj.nX;
        nU = obj.nU;
        N = obj.N;
        
        obj.Qr = Qr;
        obj.Rr = Rr;
        obj.Qrf = Qrf;
        
        dim = N-1 + N*nX + (N-1)*nU;
        cost = FunctionHandleObjective(dim,@obj.robust_cost_grad,1);
        cost.grad_method = 'user';
        obj = obj.addCost(cost, {reshape([obj.h_inds'; obj.x_inds(:,1:end-1); obj.u_inds],[],1); obj.x_inds(:,end)});
    end
    
    function obj = addRobustConstraint(obj)
        nX = obj.nX;
        nU = obj.nU;
        nw = obj.nW;
        N = obj.N;
        
        lb = repmat(obj.plant.umin,2*nw*(N-1),1);
        ub = repmat(obj.plant.umax,2*nw*(N-1),1);
        constraint = FunctionHandleConstraint(lb,ub,N-1 + N*nX + (N-1)*nU,@obj.robust_constraint_grad,1);
        constraint.grad_method = 'user';
        obj = obj.addConstraint(constraint, {reshape([obj.h_inds'; obj.x_inds(:,1:end-1); obj.u_inds],[],1); obj.x_inds(:,end)});
    end
    
    function [f,df] = forward_constraint_fun(obj,h,x0,x1,u)
      nX = obj.plant.getNumStates();
      [xdot,dxdot] = obj.plant.dynamics(0,x0,u);
      f = x1 - x0 - h*xdot;
      df = [-xdot (-eye(nX) - h*dxdot(:,2:1+nX)) eye(nX) -h*dxdot(:,nX+2:end)];
    end
    
    function [f,df] = midpoint_running_fun(obj,running_handle,h,x0,x1,u0)
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
      [f,dg] = running_handle(h,.5*(x0+x1),u0);
      
      df = [dg(:,1) .5*dg(:,2:1+nX) .5*dg(:,2:1+nX) dg(:,2+nX:1+nX+nU)];
    end
    
    function [f,df] = midpoint_constraint_fun(obj,h,x0,x1,u0)
      nX = obj.plant.getNumStates();
      [xdot,dxdot] = obj.plant.dynamics(0,.5*(x0+x1),u0);
      f = x1 - x0 - h*xdot;
      df = [-xdot (-eye(nX) - .5*h*dxdot(:,2:1+nX)) (eye(nX)- .5*h*dxdot(:,2:1+nX)) -h*dxdot(:,nX+2:end)];
    end
    
    function [c, dc] = robust_cost_grad(obj,y,xf)
        nX = obj.nX;
        nU = obj.nU;
        nw = obj.nW;
        N = obj.N;
        
        [K,A,B,G,dK,dA,dB,dG] = lqrController(obj,y,xf);
        
        c = 0;
        dc = zeros(1,(N-1)*(1+nX+nU)+nX);
        P = zeros(nX,nX);
        dP = zeros(nX*nX,(N-1)*(1+nX+nU)+nX);
        switch obj.options.integration_method
            case DirtranTrajectoryOptimization.FORWARD_EULER
                for k = 1:(N-1)
                    c = c + trace((obj.Qr + K(:,:,k)'*obj.Rr*K(:,:,k))*P);
                    
                    dcdP = vec((obj.Qr + K(:,:,k)'*obj.Rr*K(:,:,k)))';
                    dcdK = 2*vec(obj.Rr*K(:,:,k)*P)';
                    
                    dc = dc + dcdP*dP + dcdK*dK(:,:,k);
                    
                    dPdA = kron(eye(nX), A(:,:,k)*P)*comm(nX,nX) + kron(A(:,:,k)*P, eye(nX)) - kron(B(:,:,k)*K(:,:,k)*P, eye(nX)) - kron(eye(nX), B(:,:,k)*K(:,:,k)*P)*comm(nX,nX);
                    dPdB = -kron(eye(nX), A(:,:,k)*P*K(:,:,k)')*comm(nX,nU) - kron(A(:,:,k)*P*K(:,:,k)', eye(nX)) + kron(eye(nX), B(:,:,k)*K(:,:,k)*P*K(:,:,k)')*comm(nX,nU) + kron(B(:,:,k)*K(:,:,k)*P*K(:,:,k)', eye(nX));
                    dPdG = kron(eye(nX), G(:,:,k)*obj.Dinv)*comm(nX,nw) + kron(G(:,:,k)*obj.Dinv, eye(nX));
                    dPdK = -kron(B(:,:,k), A(:,:,k)*P)*comm(nU,nX) - kron(A(:,:,k)*P, B(:,:,k)) + kron(B(:,:,k)*K(:,:,k)*P, B(:,:,k)) + kron(B(:,:,k),B(:,:,k)*K(:,:,k)*P)*comm(nU,nX);
                    dPdP = kron(A(:,:,k)-B(:,:,k)*K(:,:,k), A(:,:,k)-B(:,:,k)*K(:,:,k));
                    
                    P = (A(:,:,k)-B(:,:,k)*K(:,:,k))*P*(A(:,:,k)-B(:,:,k)*K(:,:,k))' + G(:,:,k)*obj.Dinv*G(:,:,k)';
                    
                    dP = dPdP*dP + dPdK*dK(:,:,k);
                    dP(:,(k-1)*(1+nX+nU)+(1:1+nX+nU)) = dP(:,(k-1)*(1+nX+nU)+(1:1+nX+nU)) + dPdA*dA(:,:,k) + dPdB*dB(:,:,k) + dPdG*dG(:,:,k);
                end
                c = c + trace(obj.Qrf*P);
                dcdP = vec(obj.Qrf)';
                dc = dc + dcdP*dP;
            case DirtranTrajectoryOptimization.MIDPOINT
                for k = 1:(N-2)
                    c = c + trace((obj.Qr + K(:,:,k)'*obj.Rr*K(:,:,k))*P);
                    
                    dcdP = vec((obj.Qr + K(:,:,k)'*obj.Rr*K(:,:,k)))';
                    dcdK = 2*vec(obj.Rr*K(:,:,k)*P)';
                    
                    dc = dc + dcdP*dP + dcdK*dK(:,:,k);
                    
                    dPdA = kron(eye(nX), A(:,:,k)*P)*comm(nX,nX) + kron(A(:,:,k)*P, eye(nX)) - kron(B(:,:,k)*K(:,:,k)*P, eye(nX)) - kron(eye(nX), B(:,:,k)*K(:,:,k)*P)*comm(nX,nX);
                    dPdB = -kron(eye(nX), A(:,:,k)*P*K(:,:,k)')*comm(nX,nU) - kron(A(:,:,k)*P*K(:,:,k)', eye(nX)) + kron(eye(nX), B(:,:,k)*K(:,:,k)*P*K(:,:,k)')*comm(nX,nU) + kron(B(:,:,k)*K(:,:,k)*P*K(:,:,k)', eye(nX));
                    dPdG = kron(eye(nX), G(:,:,k)*obj.Dinv)*comm(nX,nw) + kron(G(:,:,k)*obj.Dinv, eye(nX));
                    dPdK = -kron(B(:,:,k), A(:,:,k)*P)*comm(nU,nX) - kron(A(:,:,k)*P, B(:,:,k)) + kron(B(:,:,k)*K(:,:,k)*P, B(:,:,k)) + kron(B(:,:,k),B(:,:,k)*K(:,:,k)*P)*comm(nU,nX);
                    dPdP = kron(A(:,:,k)-B(:,:,k)*K(:,:,k), A(:,:,k)-B(:,:,k)*K(:,:,k));
                    
                    P = (A(:,:,k)-B(:,:,k)*K(:,:,k))*P*(A(:,:,k)-B(:,:,k)*K(:,:,k))' + G(:,:,k)*obj.Dinv*G(:,:,k)';
                    
                    dP = dPdP*dP + dPdK*dK(:,:,k);
                    dP(:,(k-1)*(1+nX+nU)+(1:2*(1+nX+nU))) = dP(:,(k-1)*(1+nX+nU)+(1:2*(1+nX+nU))) + dPdA*dA(:,:,k) + dPdB*dB(:,:,k) + dPdG*dG(:,:,k);
                end
                k = N-1;
                
                c = c + trace((obj.Qr + K(:,:,k)'*obj.Rr*K(:,:,k))*P);
                
                dcdP = vec((obj.Qr + K(:,:,k)'*obj.Rr*K(:,:,k)))';
                dcdK = 2*vec(obj.Rr*K(:,:,k)*P)';
                
                dc = dc + dcdP*dP + dcdK*dK(:,:,k);
                
                dPdA = kron(eye(nX), A(:,:,k)*P)*comm(nX,nX) + kron(A(:,:,k)*P, eye(nX)) - kron(B(:,:,k)*K(:,:,k)*P, eye(nX)) - kron(eye(nX), B(:,:,k)*K(:,:,k)*P)*comm(nX,nX);
                dPdB = -kron(eye(nX), A(:,:,k)*P*K(:,:,k)')*comm(nX,nU) - kron(A(:,:,k)*P*K(:,:,k)', eye(nX)) + kron(eye(nX), B(:,:,k)*K(:,:,k)*P*K(:,:,k)')*comm(nX,nU) + kron(B(:,:,k)*K(:,:,k)*P*K(:,:,k)', eye(nX));
                dPdG = kron(eye(nX), G(:,:,k)*obj.Dinv)*comm(nX,nw) + kron(G(:,:,k)*obj.Dinv, eye(nX));
                dPdK = -kron(B(:,:,k), A(:,:,k)*P)*comm(nU,nX) - kron(A(:,:,k)*P, B(:,:,k)) + kron(B(:,:,k)*K(:,:,k)*P, B(:,:,k)) + kron(B(:,:,k),B(:,:,k)*K(:,:,k)*P)*comm(nU,nX);
                dPdP = kron(A(:,:,k)-B(:,:,k)*K(:,:,k), A(:,:,k)-B(:,:,k)*K(:,:,k));
                
                P = (A(:,:,k)-B(:,:,k)*K(:,:,k))*P*(A(:,:,k)-B(:,:,k)*K(:,:,k))' + G(:,:,k)*obj.Dinv*G(:,:,k)';
                
                dP = dPdP*dP + dPdK*dK(:,:,k);
                dP(:,(k-1)*(1+nX+nU)+(1:1+nX+nU)) = dP(:,(k-1)*(1+nX+nU)+(1:1+nX+nU)) + dPdA*dA(:,1:(1+nX+nU),k) + dPdB*dB(:,1:(1+nX+nU),k) + dPdG*dG(:,1:(1+nX+nU),k);
                dP(:,k*(1+nX+nU)+(1:nX)) = dP(:,k*(1+nX+nU)+(1:nX)) + dPdA*dA(:,(1+nX+nU+1)+(1:nX),k) + dPdB*dB(:,(1+nX+nU+1)+(1:nX),k) + dPdG*dG(:,(1+nX+nU+1)+(1:nX),k);
                
                c = c + trace(obj.Qrf*P);
                dcdP = vec(obj.Qrf)';
                dc = dc + dcdP*dP;
        end
    end

    function [c, dc] = robust_constraint_grad(obj,y,xf)
        nX = obj.nX;
        nU = obj.nU;
        nw = obj.nW;
        N = obj.N;
        
        [K,A,B,G,dK,dA,dB,dG] = lqrController(obj,y,xf);
        
        v = zeros((N-1)*nU*nw,1);
        dv = zeros((N-1)*nU*nw,(N-1)*(1+nX+nU)+nX);
        M = zeros(nX,nw);
        dM = zeros(nX*nw,(N-1)*(1+nX+nU)+nX);
        switch obj.options.integration_method
            case DirtranTrajectoryOptimization.FORWARD_EULER
                for k = 1:(N-2)
                    v((k-1)*(nU*nw)+(1:nU*nw)) = vec(K(:,:,k)*M);
                    
                    dvdK = kron(M', eye(nU));
                    dvdM = kron(eye(nw), K(:,:,k));
                    
                    dv((k-1)*(nU*nw)+(1:nU*nw),:) = dvdK*dK(:,:,k) + dvdM*dM;
                    
                    dMdA = kron(M', eye(nX));
                    dMdB = -kron((K(:,:,k)*M)', eye(nX));
                    dMdG = kron(obj.L', eye(nX));
                    dMdK = -kron(M', B(:,:,k));
                    dMdM = kron(eye(nw), A(:,:,k)-B(:,:,k)*K(:,:,k));
                    
                    dM = dMdM*dM + dMdK*dK(:,:,k);
                    dM(:,(k-1)*(1+nX+nU)+(1:(1+nX+nU))) = dM(:,(k-1)*(1+nX+nU)+(1:(1+nX+nU))) + dMdA*dA(:,:,k) + dMdB*dB(:,:,k) + dMdG*dG(:,:,k);
                    
                    M = (A(:,:,k)-B(:,:,k)*K(:,:,k))*M + G(:,:,k)*obj.L;
                end
                k = N-1;
                v((k-1)*(nU*nw)+(1:nU*nw)) = vec(K(:,:,k)*M);
                
                dvdK = kron(M', eye(nU));
                dvdM = kron(eye(nw), K(:,:,k));
                
                dv((k-1)*(nU*nw)+(1:nU*nw),:) = dvdK*dK(:,:,k) + dvdM*dM;
                
            case DirtranTrajectoryOptimization.MIDPOINT
                for k = 1:(N-2)
                    v((k-1)*(nU*nw)+(1:nU*nw)) = vec(K(:,:,k)*M);
                    
                    dvdK = kron(M', eye(nU));
                    dvdM = kron(eye(nw), K(:,:,k));
                    
                    dv((k-1)*(nU*nw)+(1:nU*nw),:) = dvdK*dK(:,:,k) + dvdM*dM;
                    
                    dMdA = kron(M', eye(nX));
                    dMdB = -kron((K(:,:,k)*M)', eye(nX));
                    dMdG = kron(obj.L', eye(nX));
                    dMdK = -kron(M', B(:,:,k));
                    dMdM = kron(eye(nw), A(:,:,k)-B(:,:,k)*K(:,:,k));
                    
                    dM = dMdM*dM + dMdK*dK(:,:,k);
                    dM(:,(k-1)*(1+nX+nU)+(1:2*(1+nX+nU))) = dM(:,(k-1)*(1+nX+nU)+(1:2*(1+nX+nU))) + dMdA*dA(:,:,k) + dMdB*dB(:,:,k) + dMdG*dG(:,:,k);
                    
                    M = (A(:,:,k)-B(:,:,k)*K(:,:,k))*M + G(:,:,k)*obj.L;
                end
                k = N-1;
                v((k-1)*(nU*nw)+(1:nU*nw)) = vec(K(:,:,k)*M);
                
                dvdK = kron(M', eye(nU));
                dvdM = kron(eye(nw), K(:,:,k));
                
                dv((k-1)*(nU*nw)+(1:nU*nw),:) = dvdK*dK(:,:,k) + dvdM*dM;
        end
        
        uinds = 1+nX+(0:N-2)'*(1+nX+nU)+kron(ones(N-1,1), (1:nU)');
        u = y(uinds);
        uc = kron(ones(nw,1), u);
        c = [uc+v(:); uc-v(:)];
        
        du = sparse(1:(N-1)*nU, uinds, ones((N-1)*nU,1),(N-1)*nU,(N-1)*(1+nX+nU)+nX);
        dc = [du+dv; du-dv];
    end
    
    function [K,A,B,G,dK,dA,dB,dG] = lqrController(obj,y,xf)
        nX = obj.nX;
        nU = obj.nU;
        nw = obj.nW;
        N = obj.N;
        
        if nargout < 5 %don't need derivatives
            %Get linearized dynamics along trajectory + derivatives
            A = zeros(nX,nX,N-1);
            B = zeros(nX,nU,N-1);
            G = zeros(nX,nw,N-1);
            switch obj.options.integration_method
                case DirtranTrajectoryOptimization.FORWARD_EULER
                    for k = 1:(N-1)
                        [~,dx1] = obj.robust_dynamics(y((k-1)*(1+nX+nU)+1),y((k-1)*(1+nX+nU)+1+(1:nX)),y((k-1)*(1+nX+nU)+1+nX+(1:nU)),zeros(1,nw));
                        A(:,:,k) = dx1(:,1+(1:nX));
                        B(:,:,k) = dx1(:,1+nX+(1:nU));
                        G(:,:,k) = dx1(:,1+nX+nU+(1:nw));
                    end
                case DirtranTrajectoryOptimization.MIDPOINT
                    for k = 1:(N-2)
                        [~,dx1] = obj.robust_dynamics(y((k-1)*(1+nX+nU)+1),.5*(y((k-1)*(1+nX+nU)+1+(1:nX))+y((k)*(1+nX+nU)+1+(1:nX))),y((k-1)*(1+nX+nU)+1+nX+(1:nU)),zeros(1,nw));
                        A(:,:,k) = dx1(:,1+(1:nX));
                        B(:,:,k) = dx1(:,1+nX+(1:nU));
                        G(:,:,k) = dx1(:,1+nX+nU+(1:nw));
                    end
                    k = N-1;
                    [~,dx] = obj.robust_dynamics(y((k-1)*(1+nX+nU)+1),.5*(y((k-1)*(1+nX+nU)+1+(1:nX))+xf),y((k-1)*(1+nX+nU)+1+nX+(1:nU)),zeros(1,nw));
                    A(:,:,k) = dx(:,1+(1:nX));
                    B(:,:,k) = dx(:,1+nX+(1:nU));
                    G(:,:,k) = dx(:,1+nX+nU+(1:nw));
            end
            
            %Solve Riccati Equation
            S = obj.Qf;
            K = zeros(nU,nX,N-1);
            for k = (N-1):-1:1
                K(:,:,k) = (B(:,:,k).'*S*B(:,:,k)+obj.R)\(B(:,:,k).'*S*A(:,:,k));
                S = obj.Q + K(:,:,k).'*obj.R*K(:,:,k) + (A(:,:,k) - B(:,:,k)*K(:,:,k)).'*S*(A(:,:,k) - B(:,:,k)*K(:,:,k));
            end
            
        else %need derivatives
            %Get dynamics derivatives along trajectory
            A = zeros(nX,nX,N-1);
            B = zeros(nX,nU,N-1);
            G = zeros(nX,nw,N-1);
            switch obj.options.integration_method
                case DirtranTrajectoryOptimization.FORWARD_EULER
                    dA = zeros(nX*nX,1+nX+nU,N-1);
                    dB = zeros(nX*nU,1+nX+nU,N-1);
                    dG = zeros(nX*nw,1+nX+nU,N-1);
                    for k = 1:(N-1)
                        [~,dx,d2x] = obj.robust_dynamics(y((k-1)*(1+nX+nU)+1),y((k-1)*(1+nX+nU)+1+(1:nX)),y((k-1)*(1+nX+nU)+1+nX+(1:nU)),zeros(1,nw));
                        A(:,:,k) = dx(:,1+(1:nX));
                        B(:,:,k) = dx(:,1+nX+(1:nU));
                        G(:,:,k) = dx(:,1+nX+nU+(1:nw));
                        dvec = reshape(d2x,nX*(1+nX+nU+nw),1+nX+nU+nw);
                        dA(:,:,k) = dvec(nX+(1:nX*nX),1:(1+nX+nU));
                        dB(:,:,k) = dvec((1+nX)*nX+(1:nX*nU),1:(1+nX+nU));
                        dG(:,:,k) = dvec((1+nX+nU)*nX+(1:nX*nw),1:(1+nX+nU));
                    end
                    
                    %Solve Riccati Equation
                    S = obj.Qf;
                    dS = zeros(nX*nX,(N-1)*(1+nX+nU)+nX);
                    K = zeros(nU,nX,N-1);
                    dK = zeros(nU*nX,(N-1)*(1+nX+nU)+nX,N-1);
                    for k = (N-1):-1:1
                        K(:,:,k) = (B(:,:,k).'*S*B(:,:,k)+obj.R)\(B(:,:,k).'*S*A(:,:,k));
                        dKdA = kron(eye(nX),(B(:,:,k)'*S*B(:,:,k)+obj.R)\B(:,:,k)'*S);
                        dKdB = kron(A(:,:,k)'*S, inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*comm(nX,nU) - kron(A(:,:,k)'*S*B(:,:,k), eye(nU))*kron(inv(B(:,:,k)'*S*B(:,:,k)+obj.R)', inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*(kron(eye(nU), B(:,:,k)'*S) + kron(B(:,:,k)'*S, eye(nU))*comm(nX,nU));
                        dKdS = kron(A(:,:,k)', (B(:,:,k)'*S*B(:,:,k)+obj.R)\B(:,:,k)') - kron(A(:,:,k)'*S*B(:,:,k), eye(nU))*kron(inv(B(:,:,k)'*S*B(:,:,k)+obj.R)', inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*kron(B(:,:,k)', B(:,:,k)');
                        dK(:,:,k) = dKdS*dS;
                        dK(:,(k-1)*(1+nX+nU)+(1:(1+nX+nU)),k) = dK(:,(k-1)*(1+nX+nU)+(1:(1+nX+nU)),k) + dKdA*dA(:,:,k) + dKdB*dB(:,:,k);
                        
                        dSdA = kron(eye(nX), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S) + kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nX))*comm(nX,nX);
                        dSdB = -kron(eye(nX), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S)*kron(K(:,:,k)', eye(nX)) - kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nX))*kron(eye(nX), K(:,:,k)')*comm(nX,nU);
                        dSdK = kron(eye(nX), K(:,:,k)'*obj.R) + kron(K(:,:,k)'*obj.R, eye(nX))*comm(nU,nX) - kron(eye(nX), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S)*kron(eye(nX), B(:,:,k)) - kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nX))*kron(B(:,:,k), eye(nX))*comm(nU,nX);
                        dSdS = kron((A(:,:,k)-B(:,:,k)*K(:,:,k))', (A(:,:,k)-B(:,:,k)*K(:,:,k))');
                        
                        S = obj.Q + K(:,:,k).'*obj.R*K(:,:,k) + (A(:,:,k) - B(:,:,k)*K(:,:,k)).'*S*(A(:,:,k) - B(:,:,k)*K(:,:,k));
                        dS = dSdS*dS + dSdK*dK(:,:,k);
                        dS(:,(k-1)*(1+nX+nU)+(1:(1+nX+nU))) = dS(:,(k-1)*(1+nX+nU)+(1:(1+nX+nU))) + dSdA*dA(:,:,k) + dSdB*dB(:,:,k);
                    end
                case DirtranTrajectoryOptimization.MIDPOINT
                    dA = zeros(nX*nX,2*(1+nX+nU),N-1);
                    dB = zeros(nX*nU,2*(1+nX+nU),N-1);
                    dG = zeros(nX*nw,2*(1+nX+nU),N-1);
                    for k = 1:(N-2)
                        [~,dx,d2x] = obj.robust_dynamics(y((k-1)*(1+nX+nU)+1),.5*(y((k-1)*(1+nX+nU)+1+(1:nX))+y((k)*(1+nX+nU)+1+(1:nX))),y((k-1)*(1+nX+nU)+1+nX+(1:nU)),zeros(1,nw));
                        A(:,:,k) = dx(:,1+(1:nX));
                        B(:,:,k) = dx(:,1+nX+(1:nU));
                        G(:,:,k) = dx(:,1+nX+nU+(1:nw));
                        dvec = reshape(d2x,nX*(1+nX+nU+nw),1+nX+nU+nw);
                        dA(:,:,k) = [dvec(nX+(1:nX*nX),1), .5*dvec(nX+(1:nX*nX),1+(1:nX)), dvec(nX+(1:nX*nX),1+nX+(1:nU)), zeros(nX*nX,1), .5*dvec(nX+(1:nX*nX),1+(1:nX)), zeros(nX*nX,nU)];
                        dB(:,:,k) = [dvec((1+nX)*nX+(1:nX*nU),1), .5*dvec((1+nX)*nX+(1:nX*nU),1+(1:nX)), dvec((1+nX)*nX+(1:nX*nU),1+nX+(1:nU)), zeros(nX*nU,1), .5*dvec((1+nX)*nX+(1:nX*nU),1+(1:nX)), zeros(nX*nU,nU)];
                        dG(:,:,k) = [dvec((1+nX+nU)*nX+(1:nX*nw),1), .5*dvec((1+nX+nU)*nX+(1:nX*nw),1+(1:nX)), dvec((1+nX+nU)*nX+(1:nX*nw),1+nX+(1:nU)), zeros(nX*nw,1), .5*dvec((1+nX+nU)*nX+(1:nX*nw),1+(1:nX)), zeros(nX*nw,nU)];
                    end
                    k = N-1;
                    [~,dx,d2x] = obj.robust_dynamics(y((k-1)*(1+nX+nU)+1),.5*(y((k-1)*(1+nX+nU)+1+(1:nX))+xf),y((k-1)*(1+nX+nU)+1+nX+(1:nU)),zeros(1,nw));
                    A(:,:,k) = dx(:,1+(1:nX));
                    B(:,:,k) = dx(:,1+nX+(1:nU));
                    G(:,:,k) = dx(:,1+nX+nU+(1:nw));
                    dvec = reshape(d2x,nX*(1+nX+nU+nw),1+nX+nU+nw);
                    dA(:,:,k) = [dvec(nX+(1:nX*nX),1), .5*dvec(nX+(1:nX*nX),1+(1:nX)), dvec(nX+(1:nX*nX),1+nX+(1:nU)), zeros(nX*nX,1), .5*dvec(nX+(1:nX*nX),1+(1:nX)), zeros(nX*nX,nU)];
                    dB(:,:,k) = [dvec((1+nX)*nX+(1:nX*nU),1), .5*dvec((1+nX)*nX+(1:nX*nU),1+(1:nX)), dvec((1+nX)*nX+(1:nX*nU),1+nX+(1:nU)), zeros(nX*nU,1), .5*dvec((1+nX)*nX+(1:nX*nU),1+(1:nX)), zeros(nX*nU,nU)];
                    dG(:,:,k) = [dvec((1+nX+nU)*nX+(1:nX*nw),1), .5*dvec((1+nX+nU)*nX+(1:nX*nw),1+(1:nX)), dvec((1+nX+nU)*nX+(1:nX*nw),1+nX+(1:nU)), zeros(nX*nw,1), .5*dvec((1+nX+nU)*nX+(1:nX*nw),1+(1:nX)), zeros(nX*nw,nU)];
                    
                    %Solve Riccati Equation
                    S = obj.Qf;
                    dS = zeros(nX*nX,(N-1)*(1+nX+nU)+nX);
                    K = zeros(nU,nX,N-1);
                    dK = zeros(nU*nX,(N-1)*(1+nX+nU)+nX,N-1);
                    
                    k = N-1;
                    K(:,:,k) = (B(:,:,k).'*S*B(:,:,k)+obj.R)\(B(:,:,k).'*S*A(:,:,k));
                    dKdA = kron(eye(nX),(B(:,:,k)'*S*B(:,:,k)+obj.R)\B(:,:,k)'*S);
                    dKdB = kron(A(:,:,k)'*S, inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*comm(nX,nU) - kron(A(:,:,k)'*S*B(:,:,k), eye(nU))*kron(inv(B(:,:,k)'*S*B(:,:,k)+obj.R)', inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*(kron(eye(nU), B(:,:,k)'*S) + kron(B(:,:,k)'*S, eye(nU))*comm(nX,nU));
                    dKdS = kron(A(:,:,k)', (B(:,:,k)'*S*B(:,:,k)+obj.R)\B(:,:,k)') - kron(A(:,:,k)'*S*B(:,:,k), eye(nU))*kron(inv(B(:,:,k)'*S*B(:,:,k)+obj.R)', inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*kron(B(:,:,k)', B(:,:,k)');
                    dK(:,:,k) = dKdS*dS;
                    dK(:,(k-1)*(1+nX+nU)+(1:(1+nX+nU)),k) = dKdA*dA(:,1:(1+nX+nU),k) + dKdB*dB(:,1:(1+nX+nU),k);
                    dK(:,k*(1+nX+nU)+(1:nX),k) = dKdA*dA(:,(1+nX+nU+1)+(1:nX),k) + dKdB*dB(:,(1+nX+nU+1)+(1:nX),k);
                    dSdA = kron(eye(nX), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S) + kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nX))*comm(nX,nX);
                    dSdB = -kron(eye(nX), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S)*kron(K(:,:,k)', eye(nX)) - kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nX))*kron(eye(nX), K(:,:,k)')*comm(nX,nU);
                    dSdK = kron(eye(nX), K(:,:,k)'*obj.R) + kron(K(:,:,k)'*obj.R, eye(nX))*comm(nU,nX) - kron(eye(nX), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S)*kron(eye(nX), B(:,:,k)) - kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nX))*kron(B(:,:,k), eye(nX))*comm(nU,nX);
                    dSdS = kron((A(:,:,k)-B(:,:,k)*K(:,:,k))', (A(:,:,k)-B(:,:,k)*K(:,:,k))');
                    S = obj.Q + K(:,:,k).'*obj.R*K(:,:,k) + (A(:,:,k) - B(:,:,k)*K(:,:,k)).'*S*(A(:,:,k) - B(:,:,k)*K(:,:,k));
                    dS = dSdS*dS + dSdK*dK(:,:,k);
                    dS(:,(k-1)*(1+nX+nU)+(1:(1+nX+nU))) = dS(:,(k-1)*(1+nX+nU)+(1:(1+nX+nU))) + dSdA*dA(:,1:(1+nX+nU),k) + dSdB*dB(:,1:(1+nX+nU),k);
                    dS(:,k*(1+nX+nU)+(1:nX)) = dS(:,k*(1+nX+nU)+(1:nX))+ dSdA*dA(:,(1+nX+nU+1)+(1:nX),k) + dSdB*dB(:,(1+nX+nU+1)+(1:nX),k);
                    for k = (N-2):-1:1
                        K(:,:,k) = (B(:,:,k).'*S*B(:,:,k)+obj.R)\(B(:,:,k).'*S*A(:,:,k));
                        dKdA = kron(eye(nX),(B(:,:,k)'*S*B(:,:,k)+obj.R)\B(:,:,k)'*S);
                        dKdB = kron(A(:,:,k)'*S, inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*comm(nX,nU) - kron(A(:,:,k)'*S*B(:,:,k), eye(nU))*kron(inv(B(:,:,k)'*S*B(:,:,k)+obj.R)', inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*(kron(eye(nU), B(:,:,k)'*S) + kron(B(:,:,k)'*S, eye(nU))*comm(nX,nU));
                        dKdS = kron(A(:,:,k)', (B(:,:,k)'*S*B(:,:,k)+obj.R)\B(:,:,k)') - kron(A(:,:,k)'*S*B(:,:,k), eye(nU))*kron(inv(B(:,:,k)'*S*B(:,:,k)+obj.R)', inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*kron(B(:,:,k)', B(:,:,k)');
                        dK(:,:,k) = dKdS*dS;
                        dK(:,(k-1)*(1+nX+nU)+(1:2*(1+nX+nU)),k) = dK(:,(k-1)*(1+nX+nU)+(1:2*(1+nX+nU)),k) + dKdA*dA(:,:,k) + dKdB*dB(:,:,k);
                        
                        dSdA = kron(eye(nX), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S) + kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nX))*comm(nX,nX);
                        dSdB = -kron(eye(nX), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S)*kron(K(:,:,k)', eye(nX)) - kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nX))*kron(eye(nX), K(:,:,k)')*comm(nX,nU);
                        dSdK = kron(eye(nX), K(:,:,k)'*obj.R) + kron(K(:,:,k)'*obj.R, eye(nX))*comm(nU,nX) - kron(eye(nX), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S)*kron(eye(nX), B(:,:,k)) - kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nX))*kron(B(:,:,k), eye(nX))*comm(nU,nX);
                        dSdS = kron((A(:,:,k)-B(:,:,k)*K(:,:,k))', (A(:,:,k)-B(:,:,k)*K(:,:,k))');
                        
                        S = obj.Q + K(:,:,k).'*obj.R*K(:,:,k) + (A(:,:,k) - B(:,:,k)*K(:,:,k)).'*S*(A(:,:,k) - B(:,:,k)*K(:,:,k));
                        dS = dSdS*dS + dSdK*dK(:,:,k);
                        dS(:,(k-1)*(1+nX+nU)+(1:2*(1+nX+nU))) = dS(:,(k-1)*(1+nX+nU)+(1:2*(1+nX+nU))) + dSdA*dA(:,:,k) + dSdB*dB(:,:,k);
                    end
            end
        end
    end
    
    function [f,df,d2f] = robust_dynamics(obj,h,x,u,w)
      % Euler integration of continuous dynamics
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
    
    function utraj = reconstructInputTrajectory(obj,z)
      % zero-order holds
      t = [0; cumsum(z(obj.h_inds(1:end-1)))];
      if size(obj.u_inds,1)>0
        u = reshape(z(obj.u_inds),[],obj.N-1);
        utraj = PPTrajectory(zoh(t,u));
        utraj = utraj.setOutputFrame(obj.plant.getInputFrame);
      else
        utraj=[];
      end
    end

    function xtraj = reconstructStateTrajectory(obj,z)
      % first-order holds
      t = [0; cumsum(z(obj.h_inds))];
      x = reshape(z(obj.x_inds),[],obj.N);
      xtraj = PPTrajectory(foh(t,x));
      xtraj = xtraj.setOutputFrame(obj.plant.getStateFrame);
    end

%     function c = robust_cost(obj,y,xf)
%         nX = obj.nX;
%         nU = obj.nU;
%         N = obj.N;
%         
%         [K,A,B,G] = lqrController(obj,y,xf);
%         
%         c = 0;
%         P = zeros(nX,nX);
%         for k = 1:(N-1)
%             c = c + trace((obj.Qr + K(:,:,k)'*obj.Rr*K(:,:,k))*P);
%             P = (A(:,:,k)-B(:,:,k)*K(:,:,k))*P*(A(:,:,k)-B(:,:,k)*K(:,:,k))' + G(:,:,k)*obj.Dinv*G(:,:,k)';
%         end
%         c = c + trace(obj.Qrf*P);
%     end
%     
%     function [c, dc] = robust_cost_fd(obj,y,xf)
%         nX = obj.nX;
%         nU = obj.nU;
%         N = obj.N;
%         
%         [K,A,B,G] = lqrController(obj,y,xf);
%         
%         c = 0;
%         P = zeros(nX,nX);
%         for k = 1:(N-1)
%             c = c + trace((obj.Qr + K(:,:,k)'*obj.Rr*K(:,:,k))*P);
%             P = (A(:,:,k)-B(:,:,k)*K(:,:,k))*P*(A(:,:,k)-B(:,:,k)*K(:,:,k))' + G(:,:,k)*obj.Dinv*G(:,:,k)';
%         end
%         c = c + trace(obj.Qrf*P);
%         
%         delta = 5e-7;
%         dc = zeros(1,length(y)+length(xf));
%         dy = zeros(size(y));
%         for k = 1:length(y)
%             dy(k) = delta;
%             dc(k) = (robust_cost(obj,y+dy,xf) - robust_cost(obj,y-dy,xf))/(2*delta);
%             dy(k) = 0;
%         end
%         dxf = zeros(size(xf));
%         for k = 1:length(xf)
%             dxf(k) = delta;
%             dc(length(y)+k) = (robust_cost(obj,y,xf+dxf) - robust_cost(obj,y,xf-dxf))/(2*delta);
%             dxf(k) = 0;
%         end
%     end

%     function [c, dc] = robust_constraint_fd(obj,y,xf)
%         %nX = obj.nX;
%         nU = obj.nU;
%         N = obj.N;
%         delta = 1e-5;
%         
%         c = robust_constraint(obj,y,xf);
%         
%         dc = zeros(2*(N-1)*nU,length(y)+length(xf));
%         dy = zeros(size(y));
%         for k = 1:length(y)
%             dy(k) = delta;
%             dc(:,k) = (robust_constraint(obj,y+dy,xf) - robust_constraint(obj,y-dy,xf))./(2*delta);
%             dy(k) = 0;
%         end
%         dxf = zeros(size(xf));
%         for k = 1:length(xf)
%             dxf(k) = delta;
%             dc(:,length(y)+k) = (robust_constraint(obj,y,xf+dxf) - robust_constraint(obj,y,xf-dxf))./(2*delta);
%             dxf(k) = 0;
%         end
%     end
%     
%     function c = robust_constraint(obj,y,xf)
%         nX = obj.nX;
%         nU = obj.nU;
%         nw = obj.nW;
%         N = obj.N;
%         
%         [K,A,B,G] = lqrController(obj,y,xf);
%         
%         v = zeros((N-1)*nU,nw);
%         M = zeros(nX,nw);
%         for k = 1:(obj.N-1)
%             v((k-1)*nU+(1:nU),nw) = K(:,:,k)*M;
%             M = (A(:,:,k)-B(:,:,k)*K(:,:,k))*M + G(:,:,k)*obj.L;
%         end
%         
%         u = y(1+nX+(0:N-2)'*(1+nX+nU)+kron(ones(N-1,1), (1:nU)'));
%         uc = kron(ones(nw,1), u);
%         c = [uc+v(:); uc-v(:)];
%     end
    
  end
end