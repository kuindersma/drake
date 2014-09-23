function runAtlasLQR(use_qp)
if ~checkDependency('gurobi')
  warning('Must have gurobi installed to run this example');
  return;
end

if nargin < 1
  use_qp = false;
end
% silence some warnings
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints')
warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits')
warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits')

options.ignore_effort_limits = false;
options.floating = false;
r = RigidBodyManipulator('urdf/atlas_minimal_contact_glued_foot.urdf',options);
% r = r.changeRootLink('r_foot',zeros(3,1),zeros(3,1));
% r.resolveConstraints(r.getInitialState());
v=r.constructVisualizer(); 

q0 = [0.09;0.4;0.65;0.23;0;0;0;0;0;0;-1.1;1.5;0.4;0;0;0.5;-0.43;1.3;-0.3;0;0;0;1.1;1.5;-0.4;0;0;0];
x0 = [q0;0*q0];

[H,C,B] = r.manipulatorDynamics(q0,0*q0);
u0 = B\C;

% u0 = ConstantTrajectory(u0);
% u0 = setOutputFrame(u0,r.getInputFrame);
% sys = cascade(u0,r);
% traj = simulate(sys,[0 1],x0);
% v.playback(traj,struct('slider',true));
  
nq = r.getNumPositions;
nv = r.getNumVelocities;
nx = r.getNumStates;
nu = r.getNumInputs;

Q = diag([10*ones(nq,1); 1*ones(nv,1)]);
R = 0.001*eye(nu);

[c,V] = tilqr(r,x0,u0,Q,R,options);

[A,B2] = linearize(r,0,x0,u0);



if use_qp
  controller_data = SharedDataHandle(struct(...
    'Q',Q,...
    'R',R,...
    'S',V.S,...
    'x0',x0,...
    'B',B,...
    'c',c,...
    'u0',u0));
  c = SimpleTIQPController(r,controller_data);
end

sys = feedback(r,c);

x0(2) = x0(2) + 0; 
x0(nq+2) = 0.55; 

tic;
traj = simulate(sys,[0 3],x0);
toc;

v.playback(traj,struct('slider',true));
save('traj_lqr_with_lims.mat','traj');

end
