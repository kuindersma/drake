function runLittleDogBalancing(use_mex)

if ~checkDependency('gurobi')
  warning('Must have gurobi installed to run this example');
  return;
end

path_handle = addpathTemporary(fullfile(getDrakePath(), 'examples', 'ZMP'));

visualize = true;

if (nargin<1); use_mex = true; end

options.floating = true;
options.dt = 0.001;
options.terrain = RigidBodyFlatTerrain();
w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
ld = LittleDog;
ld = ld.removeCollisionGroupsExcept({'feet'});
r = TimeSteppingRigidBodyManipulator(ld,options.dt);
% r = r.removeCollisionGroupsExcept({'feet'});
r = compile(r);
r = r.setStateFrame(LittleDogState(r));
r = r.setOutputFrame(LittleDogState(r));

nq = getNumPositions(r);

% set initial state to fixed point
load('data/littleDog_fp.mat');
xstar =[qstar;0*qstar];

x0 = xstar;
q0 = x0(1:nq);
kinsol = doKinematics(r,q0);

com = getCOM(r,kinsol);

% build TI-ZMP controller
footidx = ld.getFootIdx;
%terrain_pts = getTerrainContactPoints(r,footidx);
foot_pos = terrainContactPositions(r,kinsol);
comgoal = mean(foot_pos(1:2,:)')';
limp = LinearInvertedPendulum(com(3));
[~,V] = lqr(limp,comgoal);

foot_support = RigidBodySupportState(r,footidx);

ctrl_data = QPControllerData(false,struct(...
  'acceleration_input_frame',LittleDogCoordinates(r),...
  'D',-com(3)/9.81*eye(2),...
  'Qy',eye(2),...
  'S',V.S,...
  's1',zeros(4,1),...
  's2',0,...
  'x0',[comgoal;0;0],...
  'u0',zeros(2,1),...
  'y0',comgoal,...
  'qtraj',x0(1:nq),...
  'support_times',0,...
  'supports',foot_support,...
  'mu',1.0,...
  'ignore_terrain',false,...
  'constrained_dofs',[]));

% instantiate QP controller
options.slack_limit = 30.0;
options.w_qdd = 0.001*ones(nq,1);
options.w_grf = 0;
options.w_slack = 0.001;
options.debug = false;
options.use_mex = use_mex;

qp = QPController(r,{},ctrl_data,options);
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
options.contact_threshold = 0.002;
fc = FootContactBlock(r,ctrl_data,options);
ins(1).system = 2;
ins(1).input = 1;
sys = mimoFeedback(fc,sys,[],[],ins,outs);
clear ins;

% feedback PD trajectory controller
options.use_ik = false;
options.Kp = 160.0*ones(nq,1);
options.Kd = 0.8*2*sqrt(options.Kp);
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
x0 = xstar;
x0(3) = x0(3)+0.02; % drop it a bit

traj = simulate(sys,[0 10],x0);
if visualize
  playback(v,traj,struct('slider',true));
end

xf = traj.eval(traj.tspan(2));

err = norm(xf(1:6)-xstar(1:6))
if err > 0.02
  error('drakeBalancing unit test failed: error is too large');
end

end
