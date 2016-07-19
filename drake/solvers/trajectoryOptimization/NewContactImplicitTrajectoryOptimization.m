classdef NewContactImplicitTrajectoryOptimization < DirtranTrajectoryOptimization
  properties
    l_inds  
    alpha_inds  
    beta_inds
    nX
    nU
    nC
    nL
    nQ
    linc_mode = 1;
    linc_slack = 0; % slack bound for complementarity constraints
  end
  
  methods
    function obj = NewContactImplicitTrajectoryOptimization(plant,N,duration,options)
      if nargin < 4
        options = struct();
      end
      if isscalar(duration), duration=[duration,duration]; end

      obj = obj@DirtranTrajectoryOptimization(plant,N,duration,options);
      
      if isfield(options,'linc_mode')
        obj.linc_mode = options.linc_mode;
      end
      if isfield(options,'linc_slack')
        obj.linc_slack = options.linc_slack;
      end
     
      obj = obj.addContactConstraints();
      obj = obj.addComplementarityConstraints();
  
      % Ensure that all lambda values and Lagrange multipliers are non-negative
      obj = obj.addConstraint(BoundingBoxConstraint(zeros(N*obj.nL,1),inf(N*obj.nL,1)),obj.beta_inds);
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
      nb = nl;
      
      num_vars = nh + N*(nx+nu+2*nl+nb);
      obj.h_inds = (1:nh)';
      obj.x_inds = reshape(nh + (1:nx*N),nx,N);
      obj.u_inds = reshape(nh + nx*N + (1:nu*N),nu,N);
      obj.l_inds = reshape(nh + (nx+nu)*N +(1:nl*N),nl,N);
      obj.alpha_inds = reshape(nh + (nx+nu+nl)*N +(1:nl*N),nl,N); 
      obj.beta_inds = reshape(nh + (nx+nu+2*nl)*N +(1:nb*N),nb,N); 
      
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
        for j = 1:nb
          x_names{nh+(nx+nu+2*nl)*N+(i-1)*nb+j} = sprintf('beta%d[%d]',j,i);
        end
      end

      obj = obj.addDecisionVariable(num_vars,x_names);
      
      obj.nX = nx;
      obj.nU = nu;
      obj.nC = nc;
      obj.nQ = obj.plant.getNumPositions();
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
      end
    end 
        
    function obj = addContactConstraints(obj)
      nx = obj.nX;
      nu = obj.nU;
      nl = obj.nL;
      N = obj.N;
      nc = obj.nC;
     
      for i=1:N-1
        % contact force dynamics constraints
        n_vars = 1 + nx + nu + 3*nl;
        inds = {obj.h_inds(i); ...
                obj.x_inds(:,i); ...
                obj.u_inds(:,i); ...
                obj.l_inds(:,i); ...
                obj.alpha_inds(:,i); ... 
                obj.beta_inds(:,i)};
        constraint = FunctionHandleConstraint(zeros(nl,1),zeros(nl,1),n_vars,@obj.lambda_constraint_fun);
        constraint = constraint.setName(sprintf('contact_dynamics_%d',i));
        obj = obj.addConstraint(constraint, inds);        
      end
        
    end 
    
    function obj = addComplementarityConstraints(obj)
      nx = obj.nX;
      nl = obj.nL;
      N = obj.N;
      nc = obj.nC;
      % lambda >= 0, alpha >= 0, alpha'*lambda = 0
      constraint = LinearComplementarityConstraint(zeros(nl*N),zeros(nl*N,1),eye(nl*N),obj.linc_mode,obj.linc_slack);
      obj = obj.addConstraint(constraint, [obj.l_inds(:);obj.alpha_inds(:)]);

      % phi >= 0, phi'*lambda = 0
      for i=1:N
        % non-penetration constraints
        inds = {obj.x_inds(:,i)};
        constraint = FunctionHandleConstraint(zeros(nc,1),inf(nc,1),nx,@obj.phi_bound);
        constraint = constraint.setName(sprintf('phi_bound_%d',i));
        obj = obj.addConstraint(constraint, inds);        

        if i<N
          % lagrange multiplier equality constraint on vmin
          n_vars = nx + nl;
          inds = {obj.x_inds(:,i); ...
                  obj.l_inds(:,i)};     
          constraint = FunctionHandleConstraint(0,obj.linc_slack,n_vars,@obj.phi_comp);
          constraint = constraint.setName(sprintf('phi_comp_%d',i));
          obj = obj.addConstraint(constraint, inds);
        end
      end    
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
      q = x(1:nq);
      qd = x(nq+(1:nq));
      
      kinsol = doKinematics(p, q, [], struct('compute_gradients', true));

      [~,normal,d,Apts,Bpts,Aidx,Bidx,mu] = contactConstraints(p,kinsol,p.multiple_contacts);
      V = computeV(obj,normal,mu,d);
      [J,dJ] = computeJ(obj,kinsol,Aidx,Apts,Bidx,Bpts);

      dJtVl = zeros(nq,nq);
      for i=1:nq
        dJtVl(:,i) = dJ(:,(i-1)*nq+(1:nq))'*V*l;
      end
      
      [H,C,B,dH,dC] = manipulatorDynamics(p,q,qd);
      Hinv = inv(H);
      dHdq = dH(:,1:nq);
       
      tau = B*u-C;
      qdn = qd + Hinv*(tau + J'*V*l)*h;
      qn = q + qdn*h;  
      
      f = [qn;qdn];
      dtau_dq = -dC(:,1:nq);
      dtau_dqd = -dC(:,nq+(1:nq));
      dtau_du = B;
      
      dfdh = [qdn + Hinv*(tau + J'*V*l)*h; Hinv*(tau + J'*V*l)];
      
      dqdn_dq = -Hinv*matGradMult(dHdq,Hinv*(tau + J'*V*l)*h) + Hinv*(dtau_dq + dJtVl)*h;
      dqndq = eye(nq) + h*dqdn_dq;
      dfdq = [dqndq; dqdn_dq];
      
      dqdn_dqd = eye(nq) + Hinv*dtau_dqd*h;
      dfdqd = [dqdn_dqd*h; dqdn_dqd];
      
      dfdu = [Hinv*dtau_du*h^2; Hinv*dtau_du*h];
      dfdl = [Hinv*J'*V*h^2; Hinv*J'*V*h];
      
      df = [dfdh, dfdq, dfdqd, dfdu, dfdl];
    end
            
    function [f,df] = lambda_constraint_fun(obj,h,x,u,l,alpha,beta)
      p = obj.plant;
      nq = obj.nQ;
      q = x(1:nq);
      qd = x(nq+(1:nq));
      nl = obj.nL;
      nc = obj.nC;
      if obj.plant.twoD
        num_d = 2;  
      else
        num_d = 4;
      end
      
      kinsol = doKinematics(p, q, [], struct('compute_gradients', true));

      [phiC,normal,d,Apts,Bpts,Aidx,Bidx,mu,n] = contactConstraints(p,kinsol,p.multiple_contacts);

      V = computeV(obj,normal,mu,d);
      [J,dJ] = computeJ(obj,kinsol,Aidx,Apts,Bidx,Bpts);
      
      [H,C,B,dH,dC] = manipulatorDynamics(p,q,qd);
      Hinv = inv(H);
      dHdq = dH(:,1:nq);
     
      phi_ = repmat(phiC',num_d,1);
      phi_ = phi_(:);
      
      Abar = V'*J*Hinv*J'*V*h;
      tau = B*u-C;
      cbar = V'*(J*qd + J*Hinv*tau*h);
      
      f = Abar*l + cbar - alpha - phi_.*beta.*cbar;
      
      dJtVl = zeros(nq,nq);
      for i=1:nq
        dJtVl(:,i) = dJ(:,(i-1)*nq+(1:nq))'*V*l;
      end

      dphi_ = reshape(repmat(reshape(n,1,numel(n)),num_d,1),nc*num_d,numel(n)/nc);
      
      dtau_dq = -dC(:,1:nq);
      dtau_dqd = -dC(:,nq+(1:nq));
      dtau_du = B;
      
      dAbarl_dh = V'*J*Hinv*J'*V*l;
      dcbar_dh = V'*J*Hinv*tau;
      
      dAbarl_dq = V'*(matGradMult(dJ',Hinv*J'*V*l*h)' + J*(-Hinv*matGradMult(dHdq,Hinv*J'*V*l*h) + Hinv*dJtVl*h));
      dcbar_dq = V'*(matGradMult(dJ',qd + Hinv*tau*h)' + J*(-Hinv*matGradMult(dHdq,Hinv*tau*h) + Hinv*dtau_dq*h));

      dcbar_dqd = V'*(J + J*Hinv*dtau_dqd*h);
      dcbar_du = V'*J*Hinv*dtau_du*h;
      
      dfdh = dAbarl_dh + dcbar_dh - phi_.*beta.*dcbar_dh;
      dfdq = dAbarl_dq + dcbar_dq - (dphi_'*diag(beta.*cbar))' - diag(phi_.*beta)*dcbar_dq;
      dfdqd = dcbar_dqd - diag(phi_.*beta)*dcbar_dqd;
      dfdu = dcbar_du - diag(phi_.*beta)*dcbar_du;
      dfdl = Abar;
      dfdalpha = -eye(nl);
      dfdbeta = -diag(phi_.*cbar);

      df = [dfdh, dfdq, dfdqd, dfdu, dfdl, dfdalpha, dfdbeta];
    end
    
    function [f,df] = phi_bound(obj,x)
      nq = obj.plant.getNumPositions;
      q = x(1:nq);
      kinsol = doKinematics(obj.plant, q);
      [f,~,~,~,~,~,~,~,dfdq] = obj.plant.contactConstraints(kinsol,obj.plant.multiple_contacts);
      df = [dfdq,0*dfdq];
    end

    function [f,df] = phi_comp(obj,x,l)
      nq = obj.plant.getNumPositions;
      nc = obj.nC;
      q = x(1:nq);
      kinsol = doKinematics(obj.plant, q);
      [phi,normal,d,~,~,~,~,mu,dphi] = obj.plant.contactConstraints(kinsol,obj.plant.multiple_contacts);

      V = computeV(obj,normal,mu,d);
      ncell = cell(1,nc);
      for i=1:nc
        ncell{i} = normal(:,i)';
      end
      normal_mat = blkdiag(ncell{:});
      
      f = phi'*normal_mat*V*l;
      df = [(dphi'*normal_mat*V*l)', zeros(1,nq), phi'*normal_mat*V];
    end
    
    function V = computeV(obj,normal,mu,d)
      nc = obj.nC;
      nl = obj.nL;

      % TODO: clean up!
      nk = length(d);  
      V = cell(1,2*nk);
      muI = sparse(diag(mu));
      norm_mat = sparse(diag(1./sqrt(1 + mu.^2)));
      for k=1:nk,
        V{k} = (normal + d{k}*muI)*norm_mat;
        V{nk+k} = (normal - d{k}*muI)*norm_mat;
      end
      
      V = horzcat(V{:});
      I = eye(nl);
      V_cell = cell(1,nc);
      for i=1:nc
        idx_lambda = 1:nc:nl;
        V_cell{i} = V*I(idx_lambda,:)'; % basis vectors for ith contact
      end
      V = blkdiag(V_cell{:});
      
    end
    
    function [J,dJ] = computeJ(obj,kinsol,Aidx,Apts,Bidx,Bpts)
      JA = [];
      dJA = [];
      for i=1:length(Aidx)
        [~,J_,dJ_] = forwardKin(obj.plant,kinsol,Aidx(i),Apts(:,i));
        JA = [JA; J_];
        dJA = [dJA; dJ_];
      end

      JB = [];
      dJB = [];
      for i=1:length(Bidx)
        [~,J_,dJ_] = forwardKin(obj.plant,kinsol,Bidx(i),Bpts(:,i));
        JB = [JB; J_];
        dJB = [dJB; dJ_];
      end
      
      J = JA-JB;
      dJ = dJA-dJB;
    end
   
    function z0 = getInitialVars(obj,t_init,traj_init)
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
        z0(obj.u_inds) = 0.01*randn(nU,obj.N);
      end

      if isfield(traj_init,'x')
        z0(obj.x_inds) = traj_init.x.eval(t_init);
      else
        if nU>0
          if ~isfield(traj_init,'u')
            traj_init.u = setOutputFrame(PPTrajectory(foh(t_init,reshape(z0(obj.u_inds),nU,obj.N))),getInputFrame(obj.plant));
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
      
      if isfield(traj_init,'l')
        z0(obj.l_inds) = traj_init.l;
      end

      if isfield(traj_init,'alpha')
        z0(obj.alpha_inds) = traj_init.alpha;
      end
    
      if isfield(traj_init,'beta')
        z0(obj.beta_inds) = traj_init.beta;
      end
    end
    
  end

end