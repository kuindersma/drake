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
    Qrf %Robust cost matrix
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
    
    function obj = addRobustCost(obj,Qr,Rr,Qrf)
        nx = obj.nX;
        nu = obj.nU;
        N = obj.N;
        
        obj.Qr = Qr;
        obj.Rr = Rr;
        obj.Qrf = Qrf;
        
        dim = N-1 + N*nx + (N-1)*nu;
        cost = FunctionHandleObjective(dim,@obj.robust_cost_grad,1);
        cost.grad_method = 'user';
        obj = obj.addCost(cost, {reshape([obj.h_inds'; obj.x_inds(:,1:end-1); obj.u_inds],[],1); obj.x_inds(:,end)});
    end
    
    function obj = addRobustConstraint(obj)
        nx = obj.nX;
        nu = obj.nU;
        nw = obj.nW;
        N = obj.N;
        
        lb = repmat(obj.plant.umin,2*nw*(N-1),1);
        ub = repmat(obj.plant.umax,2*nw*(N-1),1);
        constraint = FunctionHandleConstraint(lb,ub,N-1 + N*nx + (N-1)*nu,@obj.robust_constraint_fd,1);
        constraint.grad_method = 'user';
        obj = obj.addConstraint(constraint, {reshape([obj.h_inds'; obj.x_inds(:,1:end-1); obj.u_inds],[],1); obj.x_inds(:,end)});
    end
    
    function c = robust_cost(obj,y,xf)
        nx = obj.nX;
        nu = obj.nU;
        N = obj.N;
        
        [K,A,B,G] = lqrController(obj,y,xf);
        
        c = 0;
        P = zeros(nx,nx);
        for k = 1:(N-1)
            c = c + trace((obj.Qr + K(:,:,k)'*obj.Rr*K(:,:,k))*P);
            P = (A(:,:,k)-B(:,:,k)*K(:,:,k))*P*(A(:,:,k)-B(:,:,k)*K(:,:,k))' + G(:,:,k)*obj.Dinv*G(:,:,k)';
        end
        c = c + trace(obj.Qrf*P);
    end
    
    function [c, dc] = robust_cost_fd(obj,y,xf)
        nx = obj.nX;
        nu = obj.nU;
        N = obj.N;
        
        [K,A,B,G] = lqrController(obj,y,xf);
        
        c = 0;
        P = zeros(nx,nx);
        for k = 1:(N-1)
            c = c + trace((obj.Qr + K(:,:,k)'*obj.Rr*K(:,:,k))*P);
            P = (A(:,:,k)-B(:,:,k)*K(:,:,k))*P*(A(:,:,k)-B(:,:,k)*K(:,:,k))' + G(:,:,k)*obj.Dinv*G(:,:,k)';
        end
        c = c + trace(obj.Qrf*P);
        
        delta = 5e-7;
        dc = zeros(1,length(y)+length(xf));
        dy = zeros(size(y));
        for k = 1:length(y)
            dy(k) = delta;
            dc(k) = (robust_cost(obj,y+dy,xf) - robust_cost(obj,y-dy,xf))/(2*delta);
            dy(k) = 0;
        end
        dxf = zeros(size(xf));
        for k = 1:length(xf)
            dxf(k) = delta;
            dc(length(y)+k) = (robust_cost(obj,y,xf+dxf) - robust_cost(obj,y,xf-dxf))/(2*delta);
            dxf(k) = 0;
        end
    end
    
    function [c, dc] = robust_cost_grad(obj,y,xf)
        nx = obj.nX;
        nu = obj.nU;
        nw = obj.nW;
        N = obj.N;
        
        [K,A,B,G,dK,dA,dB,dG] = lqrController(obj,y,xf);
        
        c = 0;
        dc = zeros(1,(N-1)*(1+nx+nu)+nx);
        P = zeros(nx,nx);
        dP = zeros(nx*nx,(N-1)*(1+nx+nu)+nx);
        switch obj.options.integration_method
            case DirtranTrajectoryOptimization.FORWARD_EULER
                for k = 1:(N-1)
                    c = c + trace((obj.Qr + K(:,:,k)'*obj.Rr*K(:,:,k))*P);
                    
                    dcdP = vec((obj.Qr + K(:,:,k)'*obj.Rr*K(:,:,k)))';
                    dcdK = 2*vec(obj.Rr*K(:,:,k)*P)';
                    
                    dc = dc + dcdP*dP + dcdK*dK(:,:,k);
                    
                    dPdA = kron(eye(nx), A(:,:,k)*P)*comm(nx,nx) + kron(A(:,:,k)*P, eye(nx)) - kron(B(:,:,k)*K(:,:,k)*P, eye(nx)) - kron(eye(nx), B(:,:,k)*K(:,:,k)*P)*comm(nx,nx);
                    dPdB = -kron(eye(nx), A(:,:,k)*P*K(:,:,k)')*comm(nx,nu) - kron(A(:,:,k)*P*K(:,:,k)', eye(nx)) + kron(eye(nx), B(:,:,k)*K(:,:,k)*P*K(:,:,k)')*comm(nx,nu) + kron(B(:,:,k)*K(:,:,k)*P*K(:,:,k)', eye(nx));
                    dPdG = kron(eye(nx), G(:,:,k)*obj.Dinv)*comm(nx,nw) + kron(G(:,:,k)*obj.Dinv, eye(nx));
                    dPdK = -kron(B(:,:,k), A(:,:,k)*P)*comm(nu,nx) - kron(A(:,:,k)*P, B(:,:,k)) + kron(B(:,:,k)*K(:,:,k)*P, B(:,:,k)) + kron(B(:,:,k),B(:,:,k)*K(:,:,k)*P)*comm(nu,nx);
                    dPdP = kron(A(:,:,k)-B(:,:,k)*K(:,:,k), A(:,:,k)-B(:,:,k)*K(:,:,k));
                    
                    P = (A(:,:,k)-B(:,:,k)*K(:,:,k))*P*(A(:,:,k)-B(:,:,k)*K(:,:,k))' + G(:,:,k)*obj.Dinv*G(:,:,k)';
                    
                    dP = dPdP*dP + dPdK*dK(:,:,k);
                    dP(:,(k-1)*(1+nx+nu)+(1:1+nx+nu)) = dP(:,(k-1)*(1+nx+nu)+(1:1+nx+nu)) + dPdA*dA(:,:,k) + dPdB*dB(:,:,k) + dPdG*dG(:,:,k);
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
                    
                    dPdA = kron(eye(nx), A(:,:,k)*P)*comm(nx,nx) + kron(A(:,:,k)*P, eye(nx)) - kron(B(:,:,k)*K(:,:,k)*P, eye(nx)) - kron(eye(nx), B(:,:,k)*K(:,:,k)*P)*comm(nx,nx);
                    dPdB = -kron(eye(nx), A(:,:,k)*P*K(:,:,k)')*comm(nx,nu) - kron(A(:,:,k)*P*K(:,:,k)', eye(nx)) + kron(eye(nx), B(:,:,k)*K(:,:,k)*P*K(:,:,k)')*comm(nx,nu) + kron(B(:,:,k)*K(:,:,k)*P*K(:,:,k)', eye(nx));
                    dPdG = kron(eye(nx), G(:,:,k)*obj.Dinv)*comm(nx,nw) + kron(G(:,:,k)*obj.Dinv, eye(nx));
                    dPdK = -kron(B(:,:,k), A(:,:,k)*P)*comm(nu,nx) - kron(A(:,:,k)*P, B(:,:,k)) + kron(B(:,:,k)*K(:,:,k)*P, B(:,:,k)) + kron(B(:,:,k),B(:,:,k)*K(:,:,k)*P)*comm(nu,nx);
                    dPdP = kron(A(:,:,k)-B(:,:,k)*K(:,:,k), A(:,:,k)-B(:,:,k)*K(:,:,k));
                    
                    P = (A(:,:,k)-B(:,:,k)*K(:,:,k))*P*(A(:,:,k)-B(:,:,k)*K(:,:,k))' + G(:,:,k)*obj.Dinv*G(:,:,k)';
                    
                    dP = dPdP*dP + dPdK*dK(:,:,k);
                    dP(:,(k-1)*(1+nx+nu)+(1:2*(1+nx+nu))) = dP(:,(k-1)*(1+nx+nu)+(1:2*(1+nx+nu))) + dPdA*dA(:,:,k) + dPdB*dB(:,:,k) + dPdG*dG(:,:,k);
                end
                k = N-1;
                
                c = c + trace((obj.Qr + K(:,:,k)'*obj.Rr*K(:,:,k))*P);
                
                dcdP = vec((obj.Qr + K(:,:,k)'*obj.Rr*K(:,:,k)))';
                dcdK = 2*vec(obj.Rr*K(:,:,k)*P)';
                
                dc = dc + dcdP*dP + dcdK*dK(:,:,k);
                
                dPdA = kron(eye(nx), A(:,:,k)*P)*comm(nx,nx) + kron(A(:,:,k)*P, eye(nx)) - kron(B(:,:,k)*K(:,:,k)*P, eye(nx)) - kron(eye(nx), B(:,:,k)*K(:,:,k)*P)*comm(nx,nx);
                dPdB = -kron(eye(nx), A(:,:,k)*P*K(:,:,k)')*comm(nx,nu) - kron(A(:,:,k)*P*K(:,:,k)', eye(nx)) + kron(eye(nx), B(:,:,k)*K(:,:,k)*P*K(:,:,k)')*comm(nx,nu) + kron(B(:,:,k)*K(:,:,k)*P*K(:,:,k)', eye(nx));
                dPdG = kron(eye(nx), G(:,:,k)*obj.Dinv)*comm(nx,nw) + kron(G(:,:,k)*obj.Dinv, eye(nx));
                dPdK = -kron(B(:,:,k), A(:,:,k)*P)*comm(nu,nx) - kron(A(:,:,k)*P, B(:,:,k)) + kron(B(:,:,k)*K(:,:,k)*P, B(:,:,k)) + kron(B(:,:,k),B(:,:,k)*K(:,:,k)*P)*comm(nu,nx);
                dPdP = kron(A(:,:,k)-B(:,:,k)*K(:,:,k), A(:,:,k)-B(:,:,k)*K(:,:,k));
                
                P = (A(:,:,k)-B(:,:,k)*K(:,:,k))*P*(A(:,:,k)-B(:,:,k)*K(:,:,k))' + G(:,:,k)*obj.Dinv*G(:,:,k)';
                
                dP = dPdP*dP + dPdK*dK(:,:,k);
                dP(:,(k-1)*(1+nx+nu)+(1:1+nx+nu)) = dP(:,(k-1)*(1+nx+nu)+(1:1+nx+nu)) + dPdA*dA(:,1:(1+nx+nu),k) + dPdB*dB(:,1:(1+nx+nu),k) + dPdG*dG(:,1:(1+nx+nu),k);
                dP(:,k*(1+nx+nu)+(1:nx)) = dP(:,k*(1+nx+nu)+(1:nx)) + dPdA*dA(:,(1+nx+nu+1)+(1:nx),k) + dPdB*dB(:,(1+nx+nu+1)+(1:nx),k) + dPdG*dG(:,(1+nx+nu+1)+(1:nx),k);
                
                c = c + trace(obj.Qrf*P);
                dcdP = vec(obj.Qrf)';
                dc = dc + dcdP*dP;
        end
    end
    
    function [c, dc] = robust_constraint_fd(obj,y,xf)
        %nx = obj.nX;
        nu = obj.nU;
        N = obj.N;
        delta = 1e-5;
        
        c = robust_constraint(obj,y,xf);
        
        dc = zeros(2*(N-1)*nu,length(y)+length(xf));
        dy = zeros(size(y));
        for k = 1:length(y)
            dy(k) = delta;
            dc(:,k) = (robust_constraint(obj,y+dy,xf) - robust_constraint(obj,y-dy,xf))./(2*delta);
            dy(k) = 0;
        end
        dxf = zeros(size(xf));
        for k = 1:length(xf)
            dxf(k) = delta;
            dc(:,length(y)+k) = (robust_constraint(obj,y,xf+dxf) - robust_constraint(obj,y,xf-dxf))./(2*delta);
            dxf(k) = 0;
        end
    end
    
    function c = robust_constraint(obj,y,xf)
        nx = obj.nX;
        nu = obj.nU;
        nw = obj.nW;
        N = obj.N;
        
        [K,A,B,G] = lqrController(obj,y,xf);
        
        v = zeros((N-1)*nu,nw);
        M = zeros(nx,nw);
        for k = 1:(obj.N-1)
            v((k-1)*nu+(1:nu),nw) = K(:,:,k)*M;
            M = (A(:,:,k)-B(:,:,k)*K(:,:,k))*M + G(:,:,k)*obj.L;
        end
        
        u = y(1+nx+(0:N-2)'*(1+nx+nu)+kron(ones(N-1,1), (1:nu)'));
        uc = kron(ones(nw,1), u);
        c = [uc+v(:); uc-v(:)];
    end

%     function [c, dc] = robust_constraint_grad(obj,y,xf)
%         nx = obj.nX;
%         nu = obj.nU;
%         nw = obj.nW;
%         N = obj.N;
%         
%         [K,A,B,G,dK,dA,dB,dG] = lqrController(obj,y,xf);
%         
%         v = zeros((N-1)*nu*nw,1);
%         dv = zeros((N-1)*nu*nw,(N-1)*(1+nx+nu)+nx);
%         M = zeros(nx,nw);
%         dM = zeros(nx*nw,(N-1)*(1+nx+nu)+nx);
%         for k = 1:(obj.N-1)
%             
%             v((k-1)*(nu*nw)+(1:nu*nw)) = vec(K(:,:,k)*M);
%             
%             dvdK = kron(M', eye(nu));
%             dvdM = kron(eye(nw), K(:,:,k));
%             
%             dv((k-1)*(nu*nw)+(1:nu*nw),:) = dvdK*dK(:,:,k) + dvdM*dM;
%             
%             dMdA = kron(M', eye(nx));
%             dMdB = -kron((K(:,:,k)*M)', eye(nx));
%             dMdG = kron(obj.L', eye(nx));
%             dMdK = -kron(M', B(:,:,k));
%             dMdM = kron(eye(nw), A(:,:,k)-B(:,:,k)*K(:,:,k));
%             
%             dM = dMdM*dM + dMdK*dK(:,:,k);
%             dM(:,(k-1)*(1+nx+nu)+(1:(1+nx+nu))) = dM(:,(k-1)*(1+nx+nu)+(1:(1+nx+nu))) + dMdA*dA(:,:,k) + dMdB*dB(:,:,k) + dMdG*dG(:,:,k);
%             
%             M = (A(:,:,k)-B(:,:,k)*K(:,:,k))*M + G(:,:,k)*obj.L;
%         end
%         
%         uinds = 1+nx+(0:N-2)'*(1+nx+nu)+kron(ones(N-1,1), (1:nu)');
%         u = y(uinds);
%         uc = kron(ones(nw,1), u);
%         c = [uc+v(:); uc-v(:)];
%         
%         du = sparse(1:(N-1)*nu, uinds, ones((N-1)*nu,1),(N-1)*nu,(N-1)*(1+nx+nu)+nx);
%         dc = [du+dv; du-dv];
%     end
    
    function [K,A,B,G,dK,dA,dB,dG] = lqrController(obj,y,xf)
        nx = obj.nX;
        nu = obj.nU;
        nw = obj.nW;
        N = obj.N;
        
        if nargout < 5 %don't need derivatives
            %Get linearized dynamics along trajectory + derivatives
            A = zeros(nx,nx,N-1);
            B = zeros(nx,nu,N-1);
            G = zeros(nx,nw,N-1);
            switch obj.options.integration_method
                case DirtranTrajectoryOptimization.FORWARD_EULER
                    for k = 1:(N-1)
                        [~,dx1] = obj.robust_dynamics(y((k-1)*(1+nx+nu)+1),y((k-1)*(1+nx+nu)+1+(1:nx)),y((k-1)*(1+nx+nu)+1+nx+(1:nu)),zeros(1,nw));
                        A(:,:,k) = dx1(:,1+(1:nx));
                        B(:,:,k) = dx1(:,1+nx+(1:nu));
                        G(:,:,k) = dx1(:,1+nx+nu+(1:nw));
                    end
                case DirtranTrajectoryOptimization.MIDPOINT
                    for k = 1:(N-2)
                        [~,dx1] = obj.robust_dynamics(y((k-1)*(1+nx+nu)+1),.5*(y((k-1)*(1+nx+nu)+1+(1:nx))+y((k)*(1+nx+nu)+1+(1:nx))),y((k-1)*(1+nx+nu)+1+nx+(1:nu)),zeros(1,nw));
                        A(:,:,k) = dx1(:,1+(1:nx));
                        B(:,:,k) = dx1(:,1+nx+(1:nu));
                        G(:,:,k) = dx1(:,1+nx+nu+(1:nw));
                    end
                    k = N-1;
                    [~,dx] = obj.robust_dynamics(y((k-1)*(1+nx+nu)+1),.5*(y((k-1)*(1+nx+nu)+1+(1:nx))+xf),y((k-1)*(1+nx+nu)+1+nx+(1:nu)),zeros(1,nw));
                    A(:,:,k) = dx(:,1+(1:nx));
                    B(:,:,k) = dx(:,1+nx+(1:nu));
                    G(:,:,k) = dx(:,1+nx+nu+(1:nw));
            end
            
            %Solve Riccati Equation
            S = obj.Qf;
            K = zeros(nu,nx,N-1);
            for k = (N-1):-1:1
                K(:,:,k) = (B(:,:,k).'*S*B(:,:,k)+obj.R)\(B(:,:,k).'*S*A(:,:,k));
                S = obj.Q + K(:,:,k).'*obj.R*K(:,:,k) + (A(:,:,k) - B(:,:,k)*K(:,:,k)).'*S*(A(:,:,k) - B(:,:,k)*K(:,:,k));
            end
            
        else %need derivatives
            %Get dynamics derivatives along trajectory
            A = zeros(nx,nx,N-1);
            B = zeros(nx,nu,N-1);
            G = zeros(nx,nw,N-1);
            switch obj.options.integration_method
                case DirtranTrajectoryOptimization.FORWARD_EULER
                    dA = zeros(nx*nx,1+nx+nu,N-1);
                    dB = zeros(nx*nu,1+nx+nu,N-1);
                    dG = zeros(nx*nw,1+nx+nu,N-1);
                    for k = 1:(N-1)
                        [~,dx,d2x] = obj.robust_dynamics(y((k-1)*(1+nx+nu)+1),y((k-1)*(1+nx+nu)+1+(1:nx)),y((k-1)*(1+nx+nu)+1+nx+(1:nu)),zeros(1,nw));
                        A(:,:,k) = dx(:,1+(1:nx));
                        B(:,:,k) = dx(:,1+nx+(1:nu));
                        G(:,:,k) = dx(:,1+nx+nu+(1:nw));
                        dvec = reshape(d2x,nx*(1+nx+nu+nw),1+nx+nu+nw);
                        dA(:,:,k) = dvec(nx+(1:nx*nx),1:(1+nx+nu));
                        dB(:,:,k) = dvec((1+nx)*nx+(1:nx*nu),1:(1+nx+nu));
                        dG(:,:,k) = dvec((1+nx+nu)*nx+(1:nx*nw),1:(1+nx+nu));
                    end
                    
                    %Solve Riccati Equation
                    S = obj.Qf;
                    dS = zeros(nx*nx,(N-1)*(1+nx+nu)+nx);
                    K = zeros(nu,nx,N-1);
                    dK = zeros(nu*nx,(N-1)*(1+nx+nu)+nx,N-1);
                    for k = (N-1):-1:1
                        K(:,:,k) = (B(:,:,k).'*S*B(:,:,k)+obj.R)\(B(:,:,k).'*S*A(:,:,k));
                        dKdA = kron(eye(nx),(B(:,:,k)'*S*B(:,:,k)+obj.R)\B(:,:,k)'*S);
                        dKdB = kron(A(:,:,k)'*S, inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*comm(nx,nu) - kron(A(:,:,k)'*S*B(:,:,k), eye(nu))*kron(inv(B(:,:,k)'*S*B(:,:,k)+obj.R)', inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*(kron(eye(nu), B(:,:,k)'*S) + kron(B(:,:,k)'*S, eye(nu))*comm(nx,nu));
                        dKdS = kron(A(:,:,k)', (B(:,:,k)'*S*B(:,:,k)+obj.R)\B(:,:,k)') - kron(A(:,:,k)'*S*B(:,:,k), eye(nu))*kron(inv(B(:,:,k)'*S*B(:,:,k)+obj.R)', inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*kron(B(:,:,k)', B(:,:,k)');
                        dK(:,:,k) = dKdS*dS;
                        dK(:,(k-1)*(1+nx+nu)+(1:(1+nx+nu)),k) = dK(:,(k-1)*(1+nx+nu)+(1:(1+nx+nu)),k) + dKdA*dA(:,:,k) + dKdB*dB(:,:,k);
                        
                        dSdA = kron(eye(nx), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S) + kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nx))*comm(nx,nx);
                        dSdB = -kron(eye(nx), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S)*kron(K(:,:,k)', eye(nx)) - kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nx))*kron(eye(nx), K(:,:,k)')*comm(nx,nu);
                        dSdK = kron(eye(nx), K(:,:,k)'*obj.R) + kron(K(:,:,k)'*obj.R, eye(nx))*comm(nu,nx) - kron(eye(nx), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S)*kron(eye(nx), B(:,:,k)) - kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nx))*kron(B(:,:,k), eye(nx))*comm(nu,nx);
                        dSdS = kron((A(:,:,k)-B(:,:,k)*K(:,:,k))', (A(:,:,k)-B(:,:,k)*K(:,:,k))');
                        
                        S = obj.Q + K(:,:,k).'*obj.R*K(:,:,k) + (A(:,:,k) - B(:,:,k)*K(:,:,k)).'*S*(A(:,:,k) - B(:,:,k)*K(:,:,k));
                        dS = dSdS*dS + dSdK*dK(:,:,k);
                        dS(:,(k-1)*(1+nx+nu)+(1:(1+nx+nu))) = dS(:,(k-1)*(1+nx+nu)+(1:(1+nx+nu))) + dSdA*dA(:,:,k) + dSdB*dB(:,:,k);
                    end
                case DirtranTrajectoryOptimization.MIDPOINT
                    dA = zeros(nx*nx,2*(1+nx+nu),N-1);
                    dB = zeros(nx*nu,2*(1+nx+nu),N-1);
                    dG = zeros(nx*nw,2*(1+nx+nu),N-1);
                    for k = 1:(N-2)
                        [~,dx,d2x] = obj.robust_dynamics(y((k-1)*(1+nx+nu)+1),.5*(y((k-1)*(1+nx+nu)+1+(1:nx))+y((k)*(1+nx+nu)+1+(1:nx))),y((k-1)*(1+nx+nu)+1+nx+(1:nu)),zeros(1,nw));
                        A(:,:,k) = dx(:,1+(1:nx));
                        B(:,:,k) = dx(:,1+nx+(1:nu));
                        G(:,:,k) = dx(:,1+nx+nu+(1:nw));
                        dvec = reshape(d2x,nx*(1+nx+nu+nw),1+nx+nu+nw);
                        dA(:,:,k) = [dvec(nx+(1:nx*nx),1), .5*dvec(nx+(1:nx*nx),1+(1:nx)), dvec(nx+(1:nx*nx),1+nx+(1:nu)), zeros(nx*nx,1), .5*dvec(nx+(1:nx*nx),1+(1:nx)), zeros(nx*nx,nu)];
                        dB(:,:,k) = [dvec((1+nx)*nx+(1:nx*nu),1), .5*dvec((1+nx)*nx+(1:nx*nu),1+(1:nx)), dvec((1+nx)*nx+(1:nx*nu),1+nx+(1:nu)), zeros(nx*nu,1), .5*dvec((1+nx)*nx+(1:nx*nu),1+(1:nx)), zeros(nx*nu,nu)];
                        dG(:,:,k) = [dvec((1+nx+nu)*nx+(1:nx*nw),1), .5*dvec((1+nx+nu)*nx+(1:nx*nw),1+(1:nx)), dvec((1+nx+nu)*nx+(1:nx*nw),1+nx+(1:nu)), zeros(nx*nw,1), .5*dvec((1+nx+nu)*nx+(1:nx*nw),1+(1:nx)), zeros(nx*nw,nu)];
                    end
                    k = N-1;
                    [~,dx,d2x] = obj.robust_dynamics(y((k-1)*(1+nx+nu)+1),.5*(y((k-1)*(1+nx+nu)+1+(1:nx))+xf),y((k-1)*(1+nx+nu)+1+nx+(1:nu)),zeros(1,nw));
                    A(:,:,k) = dx(:,1+(1:nx));
                    B(:,:,k) = dx(:,1+nx+(1:nu));
                    G(:,:,k) = dx(:,1+nx+nu+(1:nw));
                    dvec = reshape(d2x,nx*(1+nx+nu+nw),1+nx+nu+nw);
                    dA(:,:,k) = [dvec(nx+(1:nx*nx),1), .5*dvec(nx+(1:nx*nx),1+(1:nx)), dvec(nx+(1:nx*nx),1+nx+(1:nu)), zeros(nx*nx,1), .5*dvec(nx+(1:nx*nx),1+(1:nx)), zeros(nx*nx,nu)];
                    dB(:,:,k) = [dvec((1+nx)*nx+(1:nx*nu),1), .5*dvec((1+nx)*nx+(1:nx*nu),1+(1:nx)), dvec((1+nx)*nx+(1:nx*nu),1+nx+(1:nu)), zeros(nx*nu,1), .5*dvec((1+nx)*nx+(1:nx*nu),1+(1:nx)), zeros(nx*nu,nu)];
                    dG(:,:,k) = [dvec((1+nx+nu)*nx+(1:nx*nw),1), .5*dvec((1+nx+nu)*nx+(1:nx*nw),1+(1:nx)), dvec((1+nx+nu)*nx+(1:nx*nw),1+nx+(1:nu)), zeros(nx*nw,1), .5*dvec((1+nx+nu)*nx+(1:nx*nw),1+(1:nx)), zeros(nx*nw,nu)];
                    
                    %Solve Riccati Equation
                    S = obj.Qf;
                    dS = zeros(nx*nx,(N-1)*(1+nx+nu)+nx);
                    K = zeros(nu,nx,N-1);
                    dK = zeros(nu*nx,(N-1)*(1+nx+nu)+nx,N-1);
                    
                    k = N-1;
                    K(:,:,k) = (B(:,:,k).'*S*B(:,:,k)+obj.R)\(B(:,:,k).'*S*A(:,:,k));
                    dKdA = kron(eye(nx),(B(:,:,k)'*S*B(:,:,k)+obj.R)\B(:,:,k)'*S);
                    dKdB = kron(A(:,:,k)'*S, inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*comm(nx,nu) - kron(A(:,:,k)'*S*B(:,:,k), eye(nu))*kron(inv(B(:,:,k)'*S*B(:,:,k)+obj.R)', inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*(kron(eye(nu), B(:,:,k)'*S) + kron(B(:,:,k)'*S, eye(nu))*comm(nx,nu));
                    dKdS = kron(A(:,:,k)', (B(:,:,k)'*S*B(:,:,k)+obj.R)\B(:,:,k)') - kron(A(:,:,k)'*S*B(:,:,k), eye(nu))*kron(inv(B(:,:,k)'*S*B(:,:,k)+obj.R)', inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*kron(B(:,:,k)', B(:,:,k)');
                    dK(:,:,k) = dKdS*dS;
                    dK(:,(k-1)*(1+nx+nu)+(1:(1+nx+nu)),k) = dKdA*dA(:,1:(1+nx+nu),k) + dKdB*dB(:,1:(1+nx+nu),k);
                    dK(:,k*(1+nx+nu)+(1:nx),k) = dKdA*dA(:,(1+nx+nu+1)+(1:nx),k) + dKdB*dB(:,(1+nx+nu+1)+(1:nx),k);
                    dSdA = kron(eye(nx), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S) + kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nx))*comm(nx,nx);
                    dSdB = -kron(eye(nx), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S)*kron(K(:,:,k)', eye(nx)) - kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nx))*kron(eye(nx), K(:,:,k)')*comm(nx,nu);
                    dSdK = kron(eye(nx), K(:,:,k)'*obj.R) + kron(K(:,:,k)'*obj.R, eye(nx))*comm(nu,nx) - kron(eye(nx), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S)*kron(eye(nx), B(:,:,k)) - kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nx))*kron(B(:,:,k), eye(nx))*comm(nu,nx);
                    dSdS = kron((A(:,:,k)-B(:,:,k)*K(:,:,k))', (A(:,:,k)-B(:,:,k)*K(:,:,k))');
                    S = obj.Q + K(:,:,k).'*obj.R*K(:,:,k) + (A(:,:,k) - B(:,:,k)*K(:,:,k)).'*S*(A(:,:,k) - B(:,:,k)*K(:,:,k));
                    dS = dSdS*dS + dSdK*dK(:,:,k);
                    dS(:,(k-1)*(1+nx+nu)+(1:(1+nx+nu))) = dS(:,(k-1)*(1+nx+nu)+(1:(1+nx+nu))) + dSdA*dA(:,1:(1+nx+nu),k) + dSdB*dB(:,1:(1+nx+nu),k);
                    dS(:,k*(1+nx+nu)+(1:nx)) = dS(:,k*(1+nx+nu)+(1:nx))+ dSdA*dA(:,(1+nx+nu+1)+(1:nx),k) + dSdB*dB(:,(1+nx+nu+1)+(1:nx),k);
                    for k = (N-2):-1:1
                        K(:,:,k) = (B(:,:,k).'*S*B(:,:,k)+obj.R)\(B(:,:,k).'*S*A(:,:,k));
                        dKdA = kron(eye(nx),(B(:,:,k)'*S*B(:,:,k)+obj.R)\B(:,:,k)'*S);
                        dKdB = kron(A(:,:,k)'*S, inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*comm(nx,nu) - kron(A(:,:,k)'*S*B(:,:,k), eye(nu))*kron(inv(B(:,:,k)'*S*B(:,:,k)+obj.R)', inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*(kron(eye(nu), B(:,:,k)'*S) + kron(B(:,:,k)'*S, eye(nu))*comm(nx,nu));
                        dKdS = kron(A(:,:,k)', (B(:,:,k)'*S*B(:,:,k)+obj.R)\B(:,:,k)') - kron(A(:,:,k)'*S*B(:,:,k), eye(nu))*kron(inv(B(:,:,k)'*S*B(:,:,k)+obj.R)', inv(B(:,:,k)'*S*B(:,:,k)+obj.R))*kron(B(:,:,k)', B(:,:,k)');
                        dK(:,:,k) = dKdS*dS;
                        dK(:,(k-1)*(1+nx+nu)+(1:2*(1+nx+nu)),k) = dK(:,(k-1)*(1+nx+nu)+(1:2*(1+nx+nu)),k) + dKdA*dA(:,:,k) + dKdB*dB(:,:,k);
                        
                        dSdA = kron(eye(nx), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S) + kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nx))*comm(nx,nx);
                        dSdB = -kron(eye(nx), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S)*kron(K(:,:,k)', eye(nx)) - kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nx))*kron(eye(nx), K(:,:,k)')*comm(nx,nu);
                        dSdK = kron(eye(nx), K(:,:,k)'*obj.R) + kron(K(:,:,k)'*obj.R, eye(nx))*comm(nu,nx) - kron(eye(nx), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*S)*kron(eye(nx), B(:,:,k)) - kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*S, eye(nx))*kron(B(:,:,k), eye(nx))*comm(nu,nx);
                        dSdS = kron((A(:,:,k)-B(:,:,k)*K(:,:,k))', (A(:,:,k)-B(:,:,k)*K(:,:,k))');
                        
                        S = obj.Q + K(:,:,k).'*obj.R*K(:,:,k) + (A(:,:,k) - B(:,:,k)*K(:,:,k)).'*S*(A(:,:,k) - B(:,:,k)*K(:,:,k));
                        dS = dSdS*dS + dSdK*dK(:,:,k);
                        dS(:,(k-1)*(1+nx+nu)+(1:2*(1+nx+nu))) = dS(:,(k-1)*(1+nx+nu)+(1:2*(1+nx+nu))) + dSdA*dA(:,:,k) + dSdB*dB(:,:,k);
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

  end
end