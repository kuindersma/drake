classdef NewContactImplicitTrajectoryOptimization < DirtranTrajectoryOptimization
  properties
    l_inds  
    nX
    nU
    nC
    nL
    nQ
    w
  end
  
  methods
    function obj = NewContactImplicitTrajectoryOptimization(plant,N,duration,options)
      if nargin < 4
        options = struct();
      end
      if isscalar(duration), duration=[duration,duration]; end
      
      if ~isfield(options,'w')
        options.w = 1;
      end
      
      obj = obj@DirtranTrajectoryOptimization(plant,N,duration,options);

      obj.w = options.w;
      
      for i=1:N
        inds = {obj.h_inds(1); ...
                obj.x_inds(:,i); ...
                obj.u_inds(:,i); ...
                obj.l_inds(:,i)};       

        cost = FunctionHandleObjective(1+obj.nX+obj.nU+obj.nL,@obj.contact_cost_fun,1);
        obj = obj.addCost(cost,inds);
      end
      
      obj = obj.addConstraint(BoundingBoxConstraint(zeros(N*obj.nL,1),inf(N*obj.nL,1)),obj.l_inds);
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
      
      num_vars = nh + N*(nx+nu+nl);
      obj.h_inds = (1:nh)';
      obj.x_inds = reshape(nh + (1:nx*N),nx,N);
      obj.u_inds = reshape(nh + nx*N + (1:nu*N),nu,N);
      obj.l_inds = reshape(nh + (nx+nu)*N +(1:nl*N),nl,N); 
      
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
      nc = obj.nC;
      N = obj.N;
     
      for i=1:N-1
        % dynamics constraints
        n_vars = 1 + 2*nx + nu + nl;
        inds = {obj.h_inds(i); ...
                obj.x_inds(:,i); ...
                obj.x_inds(:,i+1); ...
                obj.u_inds(:,i+1); ...
                obj.l_inds(:,i+1)};
        constraint = FunctionHandleConstraint(zeros(nx,1),zeros(nx,1),n_vars,@obj.forward_constraint_fun);
        constraint = constraint.setName(sprintf('dynamics_%d',i));
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
             
    function [f,df] = forward_constraint_fun(obj,h,x0,x1,u,l)
      [xn,dxn] = forward_dynamics_fun(obj,h,x0,x1,u,l);
      
      f = xn-x1;
      dfdh = dxn(:,1);
      dfdx0 = dxn(:,1+(1:obj.nX));
      dfdx1 = dxn(:,1+obj.nX+(1:obj.nX))-eye(obj.nX);
      dfdu = dxn(:,1+obj.nX*2+(1:obj.nU));
      dfdl = dxn(:,1+obj.nX*2+obj.nU+(1:obj.nL));
      df = [dfdh, dfdx0, dfdx1, dfdu, dfdl];
    end
    
    function [f,df] = forward_dynamics_fun(obj,h,x0,x1,u,l)
      p = obj.plant;
      nq = obj.nQ;
      q0 = x0(1:nq);
      qd0 = x0(nq+(1:nq));
      q1 = x1(1:nq);
      qd1 = x1(nq+(1:nq));
      
      kinsol = doKinematics(p, q1, [], struct('compute_gradients', true));

      [~,normal,d,Apts,Bpts,Aidx,Bidx,mu] = contactConstraints(p,kinsol,p.multiple_contacts);
      V = computeV(obj,normal,mu,d);
      [J,dJ] = computeJ(obj,kinsol,Aidx,Apts,Bidx,Bpts);

      dJtVl = zeros(nq,nq);
      for i=1:nq
        dJtVl(:,i) = dJ(:,(i-1)*nq+(1:nq))'*V*l;
      end
      
      [H,C,B,dH,dC] = manipulatorDynamics(p,q1,qd1);
      Hinv = inv(H);
      dHdq1 = dH(:,1:nq);
       
      tau = B*u-C;
      qdn = qd0 + Hinv*(tau + J'*V*l)*h;
      qn = q0 + qdn*h;  
      
      f = [qn;qdn];
      
      dtau_dq1 = -dC(:,1:nq);
      dtau_dqd1 = -dC(:,nq+(1:nq));
      dtau_du = B;
      
      dqdn_dq0 = zeros(nq);
      dqn_dq0 = eye(nq);
      df_dq0 = [dqn_dq0; dqdn_dq0];
      
      dqn_dqd0 = eye(nq)*h;
      dqdn_dqd0 = eye(nq);
      df_dqd0 = [dqn_dqd0; dqdn_dqd0];
      
      dqdn_dq1 = -Hinv*matGradMult(dHdq1,Hinv*(tau + J'*V*l)*h) + Hinv*(dtau_dq1 + dJtVl)*h;
      dqn_dq1 = h*dqdn_dq1;      
      df_dq1 = [dqn_dq1; dqdn_dq1];
      
      dqdn_dqd1 = Hinv*dtau_dqd1*h;
      dqn_dqd1 = dqdn_dqd1*h;
      df_dqd1 = [dqn_dqd1; dqdn_dqd1];

      %----------
      dfdh = [qdn + Hinv*(tau + J'*V*l)*h; Hinv*(tau + J'*V*l)];
      dfdx0 = [df_dq0,df_dqd0];
      dfdx1 = [df_dq1,df_dqd1];
      dfdu = [Hinv*dtau_du*h^2; Hinv*dtau_du*h];
      dfdl = [Hinv*J'*V*h^2; Hinv*J'*V*h];
      %----------      
      
      df = [dfdh, dfdx0, dfdx1, dfdu, dfdl];
    end
            
    function [f,df] = contact_cost_fun(obj,h,x,u,l)
      p = obj.plant;
      nq = obj.nQ;
      q = x(1:nq);
      qd = x(nq+(1:nq));
      nc = obj.nC;
      
      kinsol = doKinematics(p, q, [], struct('compute_gradients', true));

      [phiC,normal,d,Apts,Bpts,Aidx,Bidx,mu,n] = contactConstraints(p,kinsol,p.multiple_contacts);

      V = computeV(obj,normal,mu,d);
      [J,dJ] = computeJ(obj,kinsol,Aidx,Apts,Bidx,Bpts);
      
      [H,C,B,dH,dC] = manipulatorDynamics(p,q,qd);
      Hinv = inv(H);
      dHdq = dH(:,1:nq);
      
      Abar = V'*J*Hinv*J'*V*h;
      tau = B*u-C;
      cbar = V'*(J*qd + J*Hinv*tau*h);
      [R,dR] = computeRegLinear(p,phiC,nc,n);
      
%       R = 0*R;
%       dR = 0*dR;
      
      f = 0.5*l'*(Abar + V'*R*V)*l + l'*cbar;
      
      dJtVl = zeros(nq,nq);
      for i=1:nq
        dJtVl(:,i) = dJ(:,(i-1)*nq+(1:nq))'*V*l;
      end

      dtau_dq = -dC(:,1:nq);
      dtau_dqd = -dC(:,nq+(1:nq));
      dtau_du = B;
      
      dAbarl_dh = V'*J*Hinv*J'*V*l;
      dcbar_dh = V'*J*Hinv*tau;
       
      dAbarl_dq = V'*(matGradMult(dJ',Hinv*J'*V*l*h)' + J*(-Hinv*matGradMult(dHdq,Hinv*J'*V*l*h) + Hinv*dJtVl*h));
      dcbar_dq = V'*(matGradMult(dJ',qd + Hinv*tau*h)' + J*(-Hinv*matGradMult(dHdq,Hinv*tau*h) + Hinv*dtau_dq*h));
 
      dcbar_dqd = V'*(J + J*Hinv*dtau_dqd*h);
      dcbar_du = V'*J*Hinv*dtau_du*h;
       
      dfdh = 0.5*l'*dAbarl_dh + l'*dcbar_dh;
      dfdq = 0.5*l'*(dAbarl_dq + V'*matGradMult(dR',V*l)') + l'*dcbar_dq;
      dfdqd = l'*dcbar_dqd;
      dfdu = l'*dcbar_du;
      dfdl = l'*(Abar + V'*R*V) + cbar';
      
      f = f*obj.w;
      df = obj.w*[dfdh, dfdq, dfdqd, dfdu, dfdl];
    end
    
    function [f,df] = phi_bound(obj,x)
      nq = obj.plant.getNumPositions;
      q = x(1:nq);
      kinsol = doKinematics(obj.plant, q);
      [f,~,~,~,~,~,~,~,dfdq] = obj.plant.contactConstraints(kinsol,obj.plant.multiple_contacts);
      df = [dfdq,0*dfdq];
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

      nu = getNumInputs(obj.plant);
      if isfield(traj_init,'u')
        z0(obj.u_inds) = traj_init.u.eval(t_init);
      else
        z0(obj.u_inds) = 0.01*randn(nu,obj.N);
      end

      if isfield(traj_init,'x')
        z0(obj.x_inds) = traj_init.x.eval(t_init);
      else
        if nu>0
          if ~isfield(traj_init,'u')
            traj_init.u = setOutputFrame(PPTrajectory(foh(t_init,reshape(z0(obj.u_inds),nu,obj.N))),getInputFrame(obj.plant));
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
    end
    
  end

end