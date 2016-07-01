classdef SmoothContactImplicitTrajectoryOptimization < DirtranTrajectoryOptimization
  properties
    l_inds  
    alpha_inds  
    nX
    nU
    nC
    nL
    nQ
    linc_mode = 1;
    linc_slack = 0;
  end
  
  methods
    function obj = SmoothContactImplicitTrajectoryOptimization(plant,N,duration,options)
      if nargin < 4
        options = struct();
      end
      if isscalar(duration), duration=[duration,duration]; end

      obj = obj@DirtranTrajectoryOptimization(plant,N,duration,options);
      
%       % Ensure that all lambda values and Lagrange multipliers are non-negative
%       obj = obj.addConstraint(BoundingBoxConstraint(zeros(N*obj.nL,1),inf(N*obj.nL,1)),obj.l_inds);
%       obj = obj.addConstraint(BoundingBoxConstraint(zeros(N*obj.nL,1),inf(N*obj.nL,1)),obj.alpha_inds);
 
      if ~isfield(options,'linc_mode')
        options.linc_mode = 1;
      end
      if ~isfield(options,'linc_slack')
        options.linc_slack = 0;
      end
      obj.linc_mode = options.linc_mode;
      obj.linc_slack = options.linc_slack;
     
      obj = obj.addComplementarityConstraint();
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
      obj.nL = nl; 
    end
    
    function obj = addDynamicConstraints(obj)
      nx = obj.nX;
      nu = obj.nU;
      nl = obj.nL;
      N = obj.N;
      nc = obj.nC;
     
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
        
      end
      
      for i=1:N
        % non-penetration constraints
        inds = {obj.x_inds(:,i)};
        constraint = FunctionHandleConstraint(zeros(nc,1),inf(nc,1),nx,@obj.phi_bound);
        constraint = constraint.setName(sprintf('phi_bound_%d',i));
        obj = obj.addConstraint(constraint, inds);        
      end
  
    end 
    
    function obj = addComplementarityConstraint(obj)
      nl = obj.nL;
      N = obj.N;
      constraint = LinearComplementarityConstraint(zeros(nl*N),zeros(nl*N,1),eye(nl*N),obj.linc_mode,obj.linc_slack);
      obj = obj.addConstraint(constraint, [obj.l_inds(:);obj.alpha_inds(:)]);
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
      nu = obj.nU;
      q = x(1:nq);
      qd = x(nq+(1:nq));
      nl = obj.nL;
      nc = obj.nC;
      dim = 3;
      
      options.kinematics_cache_ptr_to_use = [];
      options.compute_gradients = true;
      kinsol = doKinematics(p, q, [], options);

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
       
      tau = (B*u-C);
      qdn = qd + Hinv*(tau*h + J'*V*l);
      qn = q + qd*h;  
      
      f = [qn;qdn];
      dtau_dq = -dC(:,1:nq);
      dtau_dqd = -dC(:,nq+(1:nq));
      dtau_du = B;
      dfdh = [qd; Hinv*tau];
      dfdq = [eye(nq); -Hinv*matGradMult(dHdq,Hinv*(tau*h + J'*V*l)) + Hinv*(dtau_dq*h + dJtVl)];
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
            
    function [f,df] = lambda_constraint_fun(obj,h,x,u,l,alpha)
      p = obj.plant;
      nq = obj.nQ;
      q = x(1:nq);
      qd = x(nq+(1:nq));
      nl = obj.nL;
      
      kinsol = doKinematics(p, q, [], struct('compute_gradients', true));

      [phiC,normal,d,Apts,Bpts,Aidx,Bidx,mu,n] = contactConstraints(p,kinsol,p.multiple_contacts);

      V = computeV(obj,normal,mu,d);
      [J,dJ] = computeJ(obj,kinsol,Aidx,Apts,Bidx,Bpts);
      
      [H,C,B,dH,dC] = manipulatorDynamics(p,q,qd);
      Hinv = inv(H);
      dHdq = dH(:,1:nq);
     
      [R,dR] = computeSmoothingTerm(obj,phiC,n);

      tau = B*u-C;
      
      f = V'*R*V*l + V'*J*(qd + Hinv*(tau*h + J'*V*l)) - alpha;
      
      dJtVl = zeros(nq,nq);
      for i=1:nq
        dJtVl(:,i) = dJ(:,(i-1)*nq+(1:nq))'*V*l;
      end
      
      dtau_dq = -dC(:,1:nq);
      dtau_dqd = -dC(:,nq+(1:nq));
      dtau_du = B;
      
      dfdh = V'*J*Hinv*tau;
      dfdq = V'*(matGradMult(dR',V*l)' + matGradMult(dJ',qd + Hinv*(tau*h + J'*V*l))' + J*(-Hinv*matGradMult(dHdq,Hinv*(tau*h + J'*V*l)) + Hinv*(dtau_dq*h + dJtVl)));
      dfdqd = V'*(J + J*Hinv*(dtau_dqd*h));
      dfdu = V'*J*Hinv*dtau_du*h;
      dfdl = V'*(R + J*Hinv*J')*V;
      dfdalpha = -eye(nl);

      df = [dfdh, dfdq, dfdqd, dfdu, dfdl, dfdalpha];
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
   
    function [R,dR] = computeSmoothingTerm(obj,phi,n)
      p = obj.plant;
      nc = obj.nC;
      nq = obj.nQ;
      dim = 3;
      
      R_min = p.contact_smoothing_params.R_min; 
      R_max = p.contact_smoothing_params.R_max;
      phi_min = p.contact_smoothing_params.contact_threshold;
      phi_max = p.contact_smoothing_params.phi_max;
      
      k = p.contact_smoothing_params.k;
      
      [r,dr_dphi] = sig(phi,k,phi_min,phi_max,R_min,R_max);
      
      % contact smoothing matrix
      drdq = repmat(dr_dphi,1,nq).*n;

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
      dR = blkdiag(drdq_cell{:});
      
      if any(isnan(dR))
        error('NaN in dR');
      end
    end
    


  end

end