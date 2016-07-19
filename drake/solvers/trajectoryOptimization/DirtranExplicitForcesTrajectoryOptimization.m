classdef DirtranExplicitForcesTrajectoryOptimization < DirtranTrajectoryOptimization
  properties
    l_inds  
    nX
    nU
    nC
    nL
    nQ
    right_foot_pos
    left_foot_pos
    left_hand_goal_pos
  end
  
  methods
    function obj = DirtranExplicitForcesTrajectoryOptimization(plant,N,duration,rf,lf,lh,options)
      if nargin < 4
        options = struct();
      end
      if isscalar(duration), duration=[duration,duration]; end

      obj = obj@DirtranTrajectoryOptimization(plant,N,duration,options);
      
%       % Ensure that all lambda values are non-negative
      obj = obj.addConstraint(BoundingBoxConstraint(zeros(N*obj.nL,1),inf(N*obj.nL,1)),obj.l_inds);
      
      [jlmin,jlmax] = plant.getJointLimits;
      xmin = [jlmin;-inf(obj.nQ,1)];
      xmax = [jlmax;inf(obj.nQ,1)];
      
      obj = obj.addConstraint(BoundingBoxConstraint(repmat(xmin,N,1),repmat(xmax,N,1)),obj.x_inds);
      
      obj.right_foot_pos = rf;
      obj.left_foot_pos = lf;
      obj.left_hand_goal_pos = lh;
%       obj = obj.addStateConstraints;
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

    function obj = addStateConstraints(obj)
      nx = obj.nX;
      N = obj.N;
      dim = 3;
      nfeet = 2;

%       for i=1:N
%         constraint = FunctionHandleConstraint(zeros(dim*2*nfeet,1),zeros(dim*2*nfeet,1),nx,@obj.foot_position_constraint);
%         constraint = constraint.setName(sprintf('foot_position_constraint_%d',i));
%         obj = obj.addConstraint(constraint, {obj.x_inds(:,i)});
%       end
      
%       constraint = FunctionHandleConstraint(-1e-4*ones(dim,1),1e-4*ones(dim,1),nx,@obj.l_hand_constraint);
%       constraint = constraint.setName('left_hand_constraint');
%       obj = obj.addConstraint(constraint, {obj.x_inds(:,N)});
           
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
      qdn = qd + Hinv*(tau*h + J'*V*l);
      qn = q + qdn*h;  
      
      f = [qn;qdn];
      dtau_dq = -dC(:,1:nq);
      dtau_dqd = -dC(:,nq+(1:nq));
      dtau_du = B;
      dfdh = [qdn + h*Hinv*tau; Hinv*tau];
      
      dqdn_dq = -Hinv*matGradMult(dHdq,Hinv*(tau*h + J'*V*l)) + Hinv*(dtau_dq*h + dJtVl);
      dqndq = eye(nq) + h*dqdn_dq;
      dfdq = [dqndq; dqdn_dq];
      
      dqdn_dqd = eye(nq) + Hinv*(dtau_dqd*h);
      dfdqd = [dqdn_dqd*h; dqdn_dqd];
      
      dfdu = [Hinv*dtau_du*h^2; Hinv*dtau_du*h];
      dfdl = [Hinv*J'*V*h; Hinv*J'*V];
      
      df = [dfdh, dfdq, dfdqd, dfdu, dfdl];
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
    end
    
%     function [f,df] = foot_motion_constraint(obj,x)
%       p = obj.plant;
%       nq = p.getNumPositions;
%       q = x(1:nq);
%       qd = x(nq+(1:nq));
%       kinsol = doKinematics(p, q, [], struct('compute_gradients', true));
%       [~,Jr,dJr] = forwardKin(p,kinsol,p.foot_body_id.right,[0;0;0],1);
%       [~,Jl,dJl] = forwardKin(p,kinsol,p.foot_body_id.left,[0;0;0],1);
%       
%       fl = Jl*qd;
%       fr = Jr*qd;
%       f = [fl;fr];
% 
%       dfldq = matGradMult(dJl',qd)';
%       dfldqd = Jl;
%       dfrdq = matGradMult(dJr',qd)';
%       dfrdqd = Jr;
%       df = [[dfldq;dfrdq],[dfldqd;dfrdqd]];
%     end
%     

    function [f,df] = foot_position_constraint(obj,x)
      p = obj.plant;
      nq = p.getNumPositions;
      q = x(1:nq);
      kinsol = doKinematics(p, q);
      [pr,Jr] = forwardKin(p,kinsol,p.foot_body_id.right,[0;0;0],1);
      [pl,Jl] = forwardKin(p,kinsol,p.foot_body_id.left,[0;0;0],1);

      fl = pl-obj.left_foot_pos;
      fr = pr-obj.right_foot_pos;
      f = [fl;fr];

      dfldq = Jl;
      dfrdq = Jr;
      df = [[dfldq;dfrdq],zeros(length(f),nq)];
    end

%     function [f,df] = phi_bound(obj,x)
%       nq = obj.plant.getNumPositions;
%       q = x(1:nq);
%       kinsol = doKinematics(obj.plant, q);
%       [f,~,~,~,~,~,~,~,dfdq] = obj.plant.contactConstraints(kinsol,obj.plant.multiple_contacts);
%       df = [dfdq,0*dfdq];
%     end
    
    function [f,df] = l_hand_constraint(obj,x)
      p = obj.plant;
      nq = p.getNumPositions;
      q = x(1:nq);
      kinsol = doKinematics(p, q);
      [pl,Jl] = forwardKin(p,kinsol,findLinkId(p,'l_hand'),[0;0;0]);

      f = pl-obj.left_hand_goal_pos;
      df = [Jl,zeros(length(f),nq)];
    end
    
    
  end

end