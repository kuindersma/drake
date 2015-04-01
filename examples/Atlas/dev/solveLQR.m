function solveLQR(p,xtraj,utraj,ltraj,Q,R,Qf)

if nargin < 5
  Q = diag([100*ones(p.getNumPositions,1);10*ones(p.getNumVelocities,1)]);
  
%   Q(1:3,1:3) = Q(1:3,1:3)*10;
  
  R = 0.01*eye(getNumInputs(p));
  Qf = 1*Q;
end

options.use_zoh_qd = true;
options.use_zoh_u = true;

[c,Ktraj,Straj,Ptraj,Btraj,tvec,Straj_full,Ftraj] = hybridconstrainedtvlqr(p,xtraj,utraj,ltraj,Q,R,Qf,options);

keyboard;
save('data/atlas_passive_ankle_lqr.mat','xtraj','utraj','ltraj','c','Ktraj','Straj','Ptraj','Btraj','tvec','Straj_full','Ftraj','Q','R','Qf');

end

