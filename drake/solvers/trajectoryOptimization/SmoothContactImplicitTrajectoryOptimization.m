classdef SmoothContactImplicitTrajectoryOptimization < DirtranTrajectoryOptimization
  properties
    l_inds  
    alpha_inds  
    nX
    nU
    nC
    nL
    nQ
  end
  
  methods
    function obj = SmoothContactImplicitTrajectoryOptimization(plant,N,duration,options)
      if nargin < 5
        options = struct();
      end
      if isscalar(duration), duration=[duration,duration]; end

      obj = obj@DirtranTrajectoryOptimization(plant,N,duration,options);
      
      % Ensure that all lambda values and Lagrange multipliers are non-negative
      obj = obj.addConstraint(BoundingBoxConstraint(zeros(N*obj.nL,1),inf(N*obj.nL,1)),obj.l_inds);
      obj = obj.addConstraint(BoundingBoxConstraint(zeros(N*obj.nL,1),inf(N*obj.nL,1)),obj.alpha_inds);
    end
    
    function obj = setupVariables(obj, N)
      nh = N-1;
      nx = obj.plant.getNumStates();
      nu = obj.plant.getNumInputs();
      nc = obj.plant.getNumContactPairs();
      if obj.plant.twoD
        num_d = 2;  
      else
        num_d = 4;
      end
      nl = nc * num_d;
      
      num_vars = nh + N*(nx+nu+2*nl);
      obj.h_inds = (1:nh)';
      obj.x_inds = reshape(nh + (1:nx*N),nx,N);
      obj.u_inds = reshape(nh + nx*N + (1:nu*N),nu,N);
      obj.l_inds = reshape(nh + (nx+nu)*N +(1:nl*N),nl,N);
      obj.alpha_inds = reshape(nh + (nx+nu+nl)*N +(1:nl*N),nl,N); %lagrange multipliers
      
      x_names = cell(num_vars,1);
      for i = 1:N
        if(i<N)
          x_names{i} = sprintf('h[%d]',i);
        end
        for j = 1:nx
          x_names{nh+(i-1)*nx+j}=sprintf('x%d[%d]',j,i);
        end
        for j = 1:nu
          x_names{nh+nx*N+(i-1)*nu+j} = sprintf('u%d[%d]',j,i);
        end
        for j = 1:nl
          x_names{nh+(nx+nu)*N+(i-1)*nl+j} = sprintf('l%d[%d]',j,i);
          x_names{nh+(nx+nu+nl)*N+(i-1)*nl+j} = sprintf('alpha%d[%d]',j,i);
        end
      end

      obj = obj.addDecisionVariable(num_vars,x_names);
      
      obj.nX = nx;
      obj.nU = nu;
      obj.nC = nc;
      obj.nQ = obj.plant.getNumPositions();
      obj.N = N;
      obj.nL = nl;
      
    end
    
    function obj = addDynamicConstraints(obj)
      nx = obj.nX;
      nu = obj.nU;
      nl = obj.nL;
      N = obj.N;
      
      for i=1:N-1
        % dynamics constraints
        n_vars = 1 + 2*nx + nu + nl;
        inds = {obj.h_inds(i); ...
                obj.x_inds(:,i); ...
                obj.x_inds(:,i+1); ...
                obj.u_inds(:,i); ...
                obj.l_inds(:,i)};
        constraint = FunctionHandleConstraint(zeros(nx,1),zeros(nx,1),n_vars,@obj.forward_constraint_fun);
        constraint = constraint.setName(sprintf('dynamics_%d',i));
        obj = obj.addConstraint(constraint, inds);
            
        
        % contact force dynamics constraints
        n_vars = 1 + nx + nu + 2*nl;
        inds = {obj.h_inds(i); ...
                obj.x_inds(:,i); ...
                obj.u_inds(:,i); ...
                obj.l_inds(:,i); ...
                obj.alpha_inds(:,i)};
        constraint = FunctionHandleConstraint(zeros(nl,1),zeros(nl,1),n_vars,@obj.lambda_constraint_fun);
        constraint = constraint.setName(sprintf('contact_dynamics_%d',i));
        obj = obj.addConstraint(constraint, inds);
        
        % non-penetration constraints
        nc = obj.nC;
        inds = {obj.x_inds(:,i)};
        constraint = FunctionHandleConstraint(zeros(nc,1),inf(nc,1),nx,@obj.phi_bound);
        constraint = constraint.setName(sprintf('phi_bound_%d',i));
        obj = obj.addConstraint(constraint, inds);
      end
      
      % lagrange multiplier equality constraints
      inds = {obj.l_inds(:),obj.alpha_inds(:)};
      constraint = FunctionHandleConstraint(zeros(nl*N,1),zeros(nl*N,1),2*nl*N,@obj.lagrange_eq);
      constraint = constraint.setName(sprintf('lagrange_eq_%d',i));
      obj = obj.addConstraint(constraint, inds);
      
    end
 
    
    function [f,df] = forward_constraint_fun(obj,h,x0,x1,u,l)
      [xn,dxn] = forward_dynamics_fun(obj,h,x0,u,l);
      
      f = xn-x1;
      dfdh = dxn(:,1);
      dfdx0 = dxn(:,1+(1:obj.nX));
      dfdx1 = -eye(obj.nX);
      dfdu = dxn(:,1+obj.nX+(1:obj.nU));
      dfdl = dxn(:,1+obj.nX+obj.nU+(1:obj.nL));
      df = [dfdh, dfdx0, dfdx1, dfdu, dfdl];
    end
    
    function [f,df] = forward_dynamics_fun(obj,h,x,u,l)
      p = obj.plant;
      nq = obj.nQ;
      nx = obj.nX;
      nu = obj.nU;
      q = x(1:nq);
      qd = x(nq+(1:nq));
      nl = obj.nL;
      nc = obj.nC;
      
      kinsol = doKinematics(p, q, [], struct('compute_gradients', true));

      [~,normal,d,Apts,Bpts,Aidx,Bidx,mu] = contactConstraints(p,kinsol,p.multiple_contacts);

      % TODO: clean up
      nk = length(d);  
      V = cell(1,2*nk);
      muI = sparse(diag(mu));
      norm_mat = sparse(diag(1./sqrt(1 + mu.^2)));
      for k=1:nk,
        V{k} = (normal + d{k}*muI)*norm_mat;
        V{nk+k} = (normal - d{k}*muI)*norm_mat;
      end

      JA = [];
      dJA = [];
      for i=1:length(Aidx)
        [~,J_,dJ_] = forwardKin(p,kinsol,Aidx(i),Apts(:,i));
        JA = [JA; J_];
        dJA = [dJA; dJ_];
      end

      JB = [];
      dJB = [];
      for i=1:length(Bidx)
        [~,J_,dJ_] = forwardKin(p,kinsol,Bidx(i),Bpts(:,i));
        JB = [JB; J_];
        dJB = [dJB; dJ_];
      end
      
      J = JA-JB;
      dJ = dJA-dJB;
      dJtranspose = reshape(permute(reshape(dJ,size(J,1),size(J,2),[]),[2,1,3]),numel(J),[]);      
      
      V = horzcat(V{:});
      I = eye(nl);
      V_cell = cell(1,nc);
      for i=1:nc
        idx_lambda = 1:nc:nl;
        V_cell{i} = V*I(idx_lambda,:)'; % basis vectors for ith contact
      end
      V = blkdiag(V_cell{:});
      
      [H,C,B,dH,dC] = manipulatorDynamics(p,q,qd);
      Hinv = inv(H);
      dHdq = dH(:,1:nq);
      
      tau = B*u-C;
      qdn = qd + Hinv*(tau*h + J'*V*l);
      qn = q + qd*h;  
      
      f = [qn;qdn];
      dtau_dq = -dC(:,1:nq);
      dtau_dqd = -dC(:,nq+(1:nq));
      dtau_du = B;
      dfdh = [qd; Hinv*tau];
      dfdq = [eye(nq); -Hinv*matGradMult(dHdq,Hinv*(tau*h + J'*V*l)) + Hinv*(dtau_dq*h + matGradMult(dJtranspose,V*l))];
      dfdqd = [eye(nq)*h; eye(nq) + Hinv*(dtau_dqd*h)];
      dfdu = [zeros(nq,nu); Hinv*dtau_du*h];
      dfdl = [zeros(nq,nl); Hinv*J'*V];
      
      df = [dfdh, dfdq, dfdqd, dfdu, dfdl];
    end
    
    function [f,df] = phi_bound(obj,x)
      nq = obj.plant.getNumPositions;
      q = x(1:nq);
      kinsol = doKinematics(obj.plant, q);
      [f,~,~,~,~,~,~,~,dfdq] = obj.plant.contactConstraints(kinsol,obj.plant.multiple_contacts);
      df = [dfdq,0*dfdq];
    end
    
    function [f,df] = lagrange_eq(~,l,alpha)
      f = l.*alpha;
      df = [diag(alpha),diag(l)];
    end
    
    function [f,df] = lambda_constraint_fun(obj,h,x,u,l,alpha)
      p = obj.plant;
      nq = obj.nQ;
      nx = obj.nX;
      nu = obj.nU;
      q = x(1:nq);
      qd = x(nq+(1:nq));
      nl = obj.nL;
      nc = obj.nC;
      dim = 3;
      
      kinsol = doKinematics(p, q, [], struct('compute_gradients', true));

      [phiC,normal,d,Apts,Bpts,Aidx,Bidx,mu,n] = contactConstraints(p,kinsol,p.multiple_contacts);

      % TODO: clean up
      nk = length(d);  
      V = cell(1,2*nk);
      muI = sparse(diag(mu));
      norm_mat = sparse(diag(1./sqrt(1 + mu.^2)));
      for k=1:nk,
        V{k} = (normal + d{k}*muI)*norm_mat;
        V{nk+k} = (normal - d{k}*muI)*norm_mat;
      end

      JA = [];
      dJA = [];
      for i=1:length(Aidx)
        [~,J_,dJ_] = forwardKin(p,kinsol,Aidx(i),Apts(:,i));
        JA = [JA; J_];
        dJA = [dJA; dJ_];
      end

      JB = [];
      dJB = [];
      for i=1:length(Bidx)
        [~,J_,dJ_] = forwardKin(p,kinsol,Bidx(i),Bpts(:,i));
        JB = [JB; J_];
        dJB = [dJB; dJ_];
      end
      
      J = JA-JB;
      dJ = dJA-dJB;
      dJtranspose = reshape(permute(reshape(dJ,size(J,1),size(J,2),[]),[2,1,3]),numel(J),[]);      
      
      V = horzcat(V{:});
      I = eye(nl);
      V_cell = cell(1,nc);
      for i=1:nc
        idx_lambda = 1:nc:nl;
        V_cell{i} = V*I(idx_lambda,:)'; % basis vectors for ith contact
      end
      V = blkdiag(V_cell{:});
      
      [H,C,B,dH,dC] = manipulatorDynamics(p,q,qd);
      Hinv = inv(H);
      dHdq = dH(:,1:nq);
      
      % contact smoothing matrix
      R_min = 1e-3; 
      R_max = 1e4;
      k = 10;
      y = 2*(phiC-p.contact_threshold)./(p.phi_max - p.contact_threshold) - 1; % scale 
      dydq = 2/(p.phi_max - p.contact_threshold) * n;
      r = R_min + R_max./(1+exp(-k*y));
      drdy = 1./(1+exp(-k*y)).^2 .* R_max*k .* exp(-k*y);
      drdq = repmat(drdy,1,nq).*dydq;

      r = repmat(r,1,dim)';
      r = r(:);

      tmp = repmat(reshape(drdq',1,numel(drdq)),dim,1);
      drdq_cell = cell(1,nc*dim);
      cnt = 1;
      for i=1:nc
        for j=1:dim
          drdq_cell{cnt} = tmp(j,(i-1)*nq+(1:nq));
          cnt = cnt+1;
        end
      end
      
      R = diag(r);
      dRdq = blkdiag(drdq_cell{:});
      
      tau = B*u-C;
      f = V'*R*V*l + V'*J*(qd + Hinv*(tau*h + J'*V*l)) - alpha;
      
      dtau_dq = -dC(:,1:nq);
      dtau_dqd = -dC(:,nq+(1:nq));
      dtau_du = B;
      
      dfdh = V'*J*Hinv*tau;
      dfdq = V'*(matGradMult(dRdq',V*l)' + matGradMult(dJ',qd + Hinv*(tau*h + J'*V*l))' + J*(-Hinv*matGradMult(dHdq,Hinv*(tau*h + J'*V*l)) + Hinv*(dtau_dq*h + matGradMult(dJtranspose,V*l))));
      dfdqd = V'*(J + J*Hinv*(dtau_dqd*h));
      dfdu = V'*J*Hinv*dtau_du*h;
      dfdl = V'*(R + J*Hinv*J')*V;
      dfdalpha = -eye(nl);
      
      df = [dfdh, dfdq, dfdqd, dfdu, dfdl,dfdalpha];
    end
    
   
      
  end

end