classdef RobustDirtranTrajectoryOptimization < DirtranTrajectoryOptimization
    
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
  end
  
  methods
    function obj = RobustDirtranTrajectoryOptimization(plant,N,D,Q,R,Qf,duration,options)
      if nargin < 8
        options = struct();
      end
      if isscalar(duration), duration=[duration,duration]; end

      if ~isfield(options,'integration_method')
        options.integration_method = DirtranTrajectoryOptimization.MIDPOINT;
      end
      obj = obj@DirtranTrajectoryOptimization(plant,N,duration,options);
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
    
    function obj = addRobustCost(obj,Qr,Rr)
        nx = obj.nX;
        nu = obj.nU;
        N = obj.N;
        
        obj.Qr = Qr;
        obj.Rr = Rr;
        
        dim = N-1+N*(nx+nu);
        cost = FunctionHandleObjective(dim,@obj.robust_cost,0);
        cost.grad_method = 'numerical';
        obj = obj.addCost(cost, {obj.h_inds(:); obj.x_inds(:); obj.u_inds(:)});
    end
    
    function obj = addRobustConstraint(obj)
        nx = obj.nX;
        nu = obj.nU;
        N = obj.N;
        
        lb = repmat(obj.plant.umin,2*N,1);
        ub = repmat(obj.plant.umax,2*N,1);
        constraint = FunctionHandleConstraint(lb,ub,N-1+N*(nx+nu),@obj.robust_constraint,0);
        constraint.grad_method = 'numerical';
        obj = obj.addConstraint(constraint, {obj.h_inds(:); obj.x_inds(:); obj.u_inds(:)});
    end
    
    function c = robust_cost(obj,h,x,u)
        nx = obj.nX;
        N = obj.N;
        
        [K,G,A,B] = lqrController(obj,h,x,u);
        
        c = 0;
        P = zeros(nx,nx);
        for k = 1:(N-1)
            U = chol(obj.Qr + K(:,:,k)'*obj.Rr*K(:,:,k));
            c = c + h(k)*trace(U*P*U');
            P = (A(:,:,k)-B(:,:,k)*K(:,:,k))*P*(A(:,:,k)-B(:,:,k)*K(:,:,k))' + G(:,:,k)*obj.Dinv*G(:,:,k)';
        end
    end
    
    function c = robust_constraint(obj,h,x,u)
        nx = obj.nX;
        nu = obj.nU;
        nw = obj.nW;
        N = obj.N;
        
        [K,G,A,B] = lqrController(obj,h,x,u);
        
        m = zeros(N*nu,1);
        M = zeros(nx,nw);
        for k = 1:(obj.N-1)
            m((k-1)*nu+(1:nu)) = max(K(:,:,k)*M,[],2);
            M = (A(:,:,k)-B(:,:,k)*K(:,:,k))*M + G(:,:,k)*obj.L;
        end
        
        c = [u+m; u-m];
    end
    
    function [K,G,A,B] = lqrController(obj,h,x,u)
        nx = obj.nX;
        nu = obj.nU;
        nw = obj.nW;
        N = obj.N;
        
        %Get linearized dynamics along trajectory + derivatives
        A = zeros(nx,nx,N-1);
        B = zeros(nx,nu,N-1);
        G = zeros(nx,nw,N-1);
        switch obj.options.integration_method
            case DirtranTrajectoryOptimization.FORWARD_EULER
                for k = 1:(N-1)
                    [~,dx1] = obj.forward_robust_dynamics(h(k),x((k-1)*nx+(1:nx)),u((k-1)*nu+(1:nu)), zeros(1,nw));
                    A(:,:,k) = dx1(:,1+(1:nx));
                    B(:,:,k) = dx1(:,1+nx+(1:nu));
                    G(:,:,k) = dx1(:,1+nx+nu+(1:nw));
                end
            case DirtranTrajectoryOptimization.MIDPOINT
                for k = 1:(N-1)
                    [~,dx1] = obj.forward_robust_dynamics(h(k),.5*(x((k-1)*nx+(1:nx))+x(k*nx+(1:nx))),u((k-1)*nu+(1:nu)), zeros(1,nw));
                    A(:,:,k) = dx1(:,1+(1:nx));
                    B(:,:,k) = dx1(:,1+nx+(1:nu));
                    G(:,:,k) = dx1(:,1+nx+nu+(1:nw));
                end
        end
        
        %Solve Riccati Equation
        S = zeros(nx,nx,N);
        S(:,:,N) = obj.Qf;
        K = zeros(nu,nx,N);
        for k = (N-1):-1:1
            K(:,:,k) = (B(:,:,k).'*S(:,:,k+1)*B(:,:,k)+obj.R)\(B(:,:,k).'*S(:,:,k+1)*A(:,:,k));
            S(:,:,k) = obj.Q + K(:,:,k).'*obj.R*K(:,:,k) + (A(:,:,k) - B(:,:,k)*K(:,:,k)).'*S(:,:,k+1)*(A(:,:,k) - B(:,:,k)*K(:,:,k));
        end
    end
    
    function [f,df,d2f] = forward_robust_dynamics(obj,h,x,u,w)
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

    function [f,df] = midpoint_robust_dynamics(obj,h,x0,dx0,x1,dx1,u0,du0,w)
      %2nd order midpoint on dynamics, ZOH on control and disturbance inputs
      if nargout == 1
          xdot = obj.plant.dynamics_w(0,.5*(x0+dx0+x1+dx1),u0+du0,w);
          f = x0 + dx0 + h*xdot;
      else %nargout == 2
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
      end
    end

  end
end