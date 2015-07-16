function [c,Ktraj,Straj,Ptraj,Btraj,tvec,Straj_full,Ftraj,xtraj,utraj,mode_data] = hybridconstrainedtvlqr(obj,xtraj,utraj,ltraj,Q,R,Qf,options)
%HYBRIDCONSTRAINEDTVLQR
% TVLQR for a hybrid model with constraints
% @input obj The plant object
% @input xtraj Nominal state trajectory
% @input utraj Nominal input trajectory
% @input ltraj Nominal constraint force trajectory (for extracting modes)
%        The indices of the constraints which are ordered
%        [joint_limits;contact_points] where each contact point
%        is 2x1 (tangential, normal)
% @input Q
% @input R
% @input Qf
% @input options.tspan time span, defaults to using the xtraj tspan
% Performs time varying lqr on a reduced order model
% @input options.force_threshold Threshold to determine if the force
%                                is active, default 1e-3
% @input options.periodic Attemps to find a periodic solution (default 0)
% @input options.periodic_jump Jump condition where x(0) = F*x(end)
% @return c Cell of controller objects, for each hybrid mode
% @return Ktraj Cell of control gains
% @return Straj Cell of quadratic cost terms
% @return Ptraj Cell of Orthonormal bases of controllable states
% @return Btraj Cell of B(t), where xdot = A(t)x + B(t)u
% @return tvec M+1 (for M modes) vector of times in the mode sequence
%              Where the system is in mode i for t in (tvec(i),tvec(i+1))
if nargin < 8
  options = struct();
end

if isfield(options,'tspan')
  typecheck(options.tspan,'double');
  tspan = options.tspan;
else
  tspan = xtraj.tspan;
end

if isfield(options,'force_threshold')
  typecheck(options.force_threshold,'double');
else
  options.force_threshold = 1e-3;
end

if ~isfield(options,'periodic')
  options.periodic = 0;
end

if ~isfield(options,'periodic_jump')
  options.periodic_jump = eye(xtraj.dim);
end

if ~isfield(options,'use_joint_limits')
  options.use_joint_limits = false;
end

if ~isfield(options,'use_zoh_qd')
  options.use_zoh_qd = false;
end

if ~isfield(options,'use_zoh_u')
  options.use_zoh_u = false;
end

if ~isfield(options,'shift_times')
  options.shift_times = false;
end

%extract hybrid sequence
t_t = xtraj.pp.breaks;
t_t = unique([tspan(1) t_t(t_t >= tspan(1) & t_t <= tspan(2)) tspan(2)]);
x = xtraj.eval(t_t);
u = utraj.eval(t_t);
l = ltraj.eval(t_t);

if options.use_zoh_qd
  qtraj = PPTrajectory(foh(t_t,x(1:obj.getNumPositions,:)));
  qdtraj = PPTrajectory(zoh(t_t,[x(1+obj.getNumPositions:end,2:end) zeros(obj.getNumVelocities,1)]));
  xtraj = [qtraj;qdtraj];
end

if options.use_zoh_u
  utraj = PPTrajectory(zoh(t_t,u));
end

if options.use_joint_limits
  nJL = obj.getNumJointLimitConstraints;
else
  nJL = 0;
end

nC = obj.getNumContactPairs;


jl_indices = l(1:nJL,:) > options.force_threshold;
contact_indices = l([1:4:end],:) > options.force_threshold;

if options.shift_times
  mask = [ones(size(contact_indices,1),1) contact_indices(:,1:end-1)];
  contact_indices = contact_indices.*mask;
end

mode = [jl_indices;contact_indices];

% used to get constraint indices for contact points
mode_expand = zeros(nJL+2*nC,nJL+nC);
mode_expand(1:nJL,1:nJL) = eye(nJL);
mode_expand(nJL+1:end,nJL+1:end) = kron(eye(nC),[1;1]);

constraint_bool = mode_expand*mode;
constraint_list = (1:nJL+2*nC)';

mode_num = sum(repmat(2.^(0:size(mode,1)-1)',1,size(mode,2)).*mode);


for i=1:size(mode,2)-1,
  [phi,normal,d,xA,xB,idxA,idxB,mu,n,D] = contactConstraints(obj,x(1:obj.getNumPositions,i+1));
  J = zeros(2*length(phi),obj.getNumPositions);
  constraint_ind = constraint_list(constraint_bool(:,i) == 1);
  J(1:2:end,:) = D{1};
  J(2:2:end,:) = n;
  if length(constraint_ind) > 0
    J_sub = J(constraint_ind(1),:);
    for k=2:length(constraint_ind),
      J_sub = [J_sub;J(constraint_ind(k),:)];
      %       if rank(J_sub) ~= size(J_sub,1),
      if cond(J_sub) > 1e2
        J_sub = J_sub(1:end-1,:);
        constraint_bool(constraint_ind(k),i) = 0;
      end
    end
  end
end

mode_ind = 1;
mode_start = 1;
tvec = tspan(1);
for i=1:size(mode,2)-1,
  if i == size(mode,2)-1 || (mode_num(:,i+1) ~= mode_num(:,i))
    mode_end = i + 1;
    
    mode_data{mode_ind}.tspan = [t_t(mode_start) t_t(mode_end)];
    mode_data{mode_ind}.constraint_ind = constraint_list(constraint_bool(:,i) == 1);
    tvec = [tvec t_t(mode_end)];
    
    if i ~= size(mode,2) - 1
      next_constraint_ind = constraint_list(constraint_bool(:,i+1) == 1);
      jumpfn = @(xm) jump(obj,xm,next_constraint_ind,options);
      [~,F] = geval(jumpfn,x(:,i+1),struct('grad_method','numerical'));
      mode_data{mode_ind}.F = full(F);
    end
    
    mode_ind = mode_ind + 1;
    mode_start = i+1;
  end
end

display(sprintf('Identified %d modes, with timings:',length(mode_data)));
for i=1:length(mode_data),
  display(sprintf('t %f to %f',mode_data{i}.tspan(1),mode_data{i}.tspan(2)));
end

if options.periodic
  jmax = 3;
else
  jmax = 1;
end
  
Qfi = Qf;
for j = 1:jmax;
  for i = length(mode_data):-1:1,
    ts = mode_data{i}.tspan;
    
    % construct constrained plant
    if ~isempty(jl_indices)
      error('Joint limits not implemented yet')
    end
    k = 1;
    constrained_plant = obj;
    [phi,normal,d,xA,xB,idxA,idxB,mu] = obj.contactConstraints(x(1:obj.getNumPositions,1),false,struct('terrain_only',true));
    % swap A and B
    tmp = xA;
    xA = xB;
    xB = tmp;
    tmp = idxA;
    idxA = idxB;
    idxB = tmp;
    while k <= length(mode_data{i}.constraint_ind)
      if false
        if mod(mode_data{i}.constraint_ind(k),2)~=0
          ind = (mode_data{i}.constraint_ind(k)+1)/2;
        else
          ind = (mode_data{i}.constraint_ind(k))/2;
        end
        % Use both first and second (relative position)
        position_fun = drakeFunction.kinematic.RelativePosition(constrained_plant,idxB(ind),idxA(ind),obj.T_2D_to_3D'*xB(:,ind));
        position_constraint = DrakeFunctionConstraint(obj.T_2D_to_3D'*xA(:,ind),obj.T_2D_to_3D'*xA(:,ind),position_fun);
        position_constraint.grad_level = 2;
        constrained_plant = constrained_plant.addPositionEqualityConstraint(position_constraint);
        k = k+2;
      else
        
        if mod(mode_data{i}.constraint_ind(k),2)~=0
          ind = (mode_data{i}.constraint_ind(k)+1)/2;
          if k < length(mode_data{i}.constraint_ind) && (mode_data{i}.constraint_ind(k+1) - mode_data{i}.constraint_ind(k)) == 1
            % Use both first and second (relative position)
            position_fun = drakeFunction.kinematic.RelativePosition(constrained_plant,idxB(ind),idxA(ind),obj.T_2D_to_3D'*xB(:,ind));
            position_constraint = DrakeFunctionConstraint(obj.T_2D_to_3D'*xA(:,ind),obj.T_2D_to_3D'*xA(:,ind),position_fun);
            position_constraint.grad_level = 2;
            constrained_plant = constrained_plant.addPositionEqualityConstraint(position_constraint);
            k = k+2;
          else
            % Odd -> only use first element
            posA = obj.T_2D_to_3D(:,1)'*xA(:,ind);
            position_fun = drakeFunction.kinematic.RelativePosition(constrained_plant,idxB(ind),idxA(ind),obj.T_2D_to_3D'*xB(:,ind),1);
            position_constraint = DrakeFunctionConstraint(posA,posA,position_fun);
            position_constraint.grad_level = 2;
            constrained_plant = constrained_plant.addPositionEqualityConstraint(position_constraint);
            k = k+1;
          end
        else
          ind = (mode_data{i}.constraint_ind(k))/2;
          % Even -> only use second element
          posA = obj.T_2D_to_3D(:,2)'*xA(:,ind);
          position_fun = drakeFunction.kinematic.RelativePosition(constrained_plant,idxB(ind),idxA(ind),obj.T_2D_to_3D'*xB(:,ind),2);
          position_constraint = DrakeFunctionConstraint(posA,posA,position_fun);
          position_constraint.grad_level = 2;
          constrained_plant = constrained_plant.addPositionEqualityConstraint(position_constraint);
          k = k+1;
        end
      end
    end
    mode_data{i}.constraint_ind;
    
    [c{i},Ktraj{i},Straj{i},Ptraj{i},Btraj{i},Ftraj{i},Straj_full{i}] = constrainedtvlqr(constrained_plant,xtraj,utraj,Q,R,Qfi,struct('tspan',ts));
    P0 = Ptraj{i}.eval(ts(1));
    S0 = Straj{i}.eval(ts(1));
    Qfi = pinv(P0*P0')*P0*S0*P0'*pinv(P0*P0')';  %extract least squares full rank solution
    
    if i > 1 || options.periodic
      if i~=1
        Qfi = mode_data{i-1}.F'*Qfi*mode_data{i-1}.F;
      else
        Qfi = options.periodic_jump'*Qfi*options.periodic_jump;
      end
      Qfi = (Qfi+Qfi')/2 + 1e-2*eye(size(Qfi,1));
      min(eig(Qfi))
    end
  end
end
end

function xp = jump(obj,xm,constraint_ind,options)
%Computes xp, the post impact state of jumping to this mode
q=xm(1:obj.getNumPositions()); qd=xm((obj.getNumPositions()+1):end);
H = manipulatorDynamics(obj,q,qd);
Hinv = inv(H);
if (size(constraint_ind) > 0)
  
  if options.use_joint_limits
    [~,J_limit] = obj.jointLimitConstraints(q);
  else
    J_limit = zeros(0,length(q));
  end
  
  
  [phi_n,normal,d,xA,xB,idxA,idxB,mu,n,D] = contactConstraints(obj,q);
  
  J = zeros(2*length(phi_n), obj.getNumPositions());
  J(1:2:end,:) = D{1};
  J(2:2:end,:) = n;
  
  J_full = [J_limit;J];
  
  J_sub = J_full(constraint_ind,:);
  
  qdp = qd-Hinv*J_sub'*inv(J_sub*Hinv*J_sub')*J_sub*qd;
  xp = [q;qdp];
else
  xp = xm;
end
end
