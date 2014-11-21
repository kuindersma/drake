function runAtlasMomentumLQR(use_mex)

if ~checkDependency('gurobi')
  warning('Must have gurobi installed to run this example');
  return;
end

path_handle = addpathTemporary({fullfile(getDrakePath,'examples','Atlas','controllers'),...
                                fullfile(getDrakePath,'examples','Atlas','frames')});

visualize = true;

if (nargin<1); use_mex = true; end

% silence some warnings
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints')
warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits')
warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits')

options.floating = true;
options.dt = 0.002;
r = Atlas('urdf/atlas_minimal_contact.urdf',options);
r = r.removeCollisionGroupsExcept({'heel','toe'});
r = compile(r);

nq = getNumPositions(r);

% set initial state to fixed point
load('data/atlas_fp.mat');
r = r.setInitialState(xstar);

x0 = xstar;
q0 = x0(1:nq);
kinsol = doKinematics(r,q0);

% build time-invariant momentum LQR controller 
H = manipulatorDynamics(r,q0,0*q0);
% Hdot = 0 since qd0=0
Hinv = inv(H);
A = zeros(nq);
B = H;
Ag = getCMM(r,kinsol);
C = Ag*Hinv;
R = 0.001*eye(nq);

options.Qy = diag([0.25 0.25 0.25 1 1 1]);
tisys = LinearSystem(A,B,[],[],C,[]);
[c,V] = tilqr(tisys,zeros(nq,1),zeros(nq,1),eps*eye(nq),R,options);

foot_support = RigidBodySupportState(r,find(~cellfun(@isempty,strfind(r.getLinkNames(),'foot'))));    

ctrl_data = QPControllerData(false,struct(...
  'acceleration_input_frame',AtlasCoordinates(r),...
  'A',A,...
  'B',B,...
  'C',C,...
  'D',zeros(6,nq),...
  'Qy',options.Qy,...
  'R',R,...
  'S',V.S,...
  's1',zeros(nq,1),...
  's2',0,...
  'x0',zeros(nq,1),...
  'u0',zeros(nq,1),...
  'y0',zeros(6,1),...
  'qtraj',q0,...
  'support_times',0,...
  'supports',foot_support,...
  'mu',1.0,...
  'ignore_terrain',false,...
  'constrained_dofs',[]));

% instantiate QP controller
options.slack_limit = 10.0;
options.w_qdd = 0.001*ones(nq,1);
options.w_grf = 0;
options.w_slack = 1;
options.debug = false;
options.use_mex = use_mex;

qp = MomentumQPController(r,{},ctrl_data,options);
clear options;

% feedback QP controller with atlas
ins(1).system = 1;
ins(1).input = 2;
ins(2).system = 1;
ins(2).input = 3;
outs(1).system = 2;
outs(1).output = 1;
sys = mimoFeedback(qp,r,[],[],ins,outs);
clear ins;

% feedback foot contact detector with QP/atlas
options.use_lcm=false;
options.contact_threshold = 0.0005;
fc = FootContactBlock(r,ctrl_data,options);
ins(1).system = 2;
ins(1).input = 1;
sys = mimoFeedback(fc,sys,[],[],ins,outs);
clear ins;

% feedback PD trajectory controller 
options.use_ik = false;
pd = IKPDBlock(r,ctrl_data,options);
ins(1).system = 1;
ins(1).input = 1;
sys = mimoFeedback(pd,sys,[],[],ins,outs);
clear ins;

qt = QTrajEvalBlock(r,ctrl_data);
sys = mimoFeedback(qt,sys,[],[],[],outs);

if visualize
  v = r.constructVisualizer;
  v.display_dt = 0.01;
  S=warning('off','Drake:DrakeSystem:UnsupportedSampleTime');
  output_select(1).system=1;
  output_select(1).output=1;
  sys = mimoCascade(sys,v,[],[],output_select);
  warning(S);
end
x0(3) = 1.3; % drop it a bit

traj = simulate(sys,[0 5],x0);
if visualize
  playback(v,traj,struct('slider',true));
end

xf = traj.eval(traj.tspan(2));

err = norm(xf(1:6)-xstar(1:6))
if err > 0.02
  error('drakeBalancing unit test failed: error is too large');
end

end
