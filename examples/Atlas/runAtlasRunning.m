function runAtlasRunning(use_mex,use_angular_momentum)

if ~checkDependency('gurobi')
  warning('Must have gurobi installed to run this example');
  return;
end

if (nargin<1); use_mex = true; end
if (nargin<2); use_angular_momentum = false; end

path_handle = addpathTemporary({fullfile(getDrakePath,'examples','Atlas','controllers'),...
                                fullfile(getDrakePath,'examples','Atlas','frames')});

% silence some warnings
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints')
warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits')
warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits')

options.floating = true;
options.dt = 0.001;
options.ignore_effort_limits = true;
r = Atlas('urdf/atlas_minimal_contact.urdf',options);
r = r.removeCollisionGroupsExcept({'heel','toe'});
r = compile(r);

v = r.constructVisualizer;
v.display_dt = 0.005;

% load in running trajectory
load('data/results_2m_2mps_midstrike.mat');

xtraj = sol.xtraj_three;
ts = unique(xtraj.getBreaks);
r.setInitialState(xtraj.eval(0));
x_knots = xtraj.eval(ts);

link_indices = [findLinkInd(r,'l_foot'), findLinkInd(r,'r_foot'),...
    findLinkInd(r,'l_hand'), findLinkInd(r,'r_hand'), ...
    findLinkInd(r,'pelvis'), findLinkInd(r,'utorso')];

body_knots = cell(length(link_indices));
for j=1:length(link_indices)
  body_knots{j} = zeros(6,length(ts));
end

nq = getNumPositions(r);
A_knots = zeros(nq,nq,length(ts));
B_knots = zeros(nq,nq,length(ts));
C_knots = zeros(6,nq,length(ts));
x0_knots = zeros(nq,length(ts));
u0_knots = zeros(nq,length(ts));
y0_knots = zeros(6,length(ts));
for i=1:length(ts)
  qi = x_knots(1:nq,i);
  vi = x_knots(nq+(1:nq),i);
  kinsol = doKinematics(r,qi,true,true,vi);
  for j=1:length(link_indices)
    body_knots{j}(:,i) = forwardKin(r,kinsol,link_indices(j),[0;0;0],1);
  end
  [H,~,~,dH] = manipulatorDynamics(r,qi,vi,true);
  Hdot = matGradMult(dH(:,1:nq),vi);
  Hinv = inv(H);
  A_knots(:,:,i) = Hdot*Hinv;
  B_knots(:,:,i) = H;
  Ag = getCMM(r,kinsol,vi);
  C_knots(:,:,i) = Ag*Hinv;
  x0_knots(:,i) = H*vi;
  y0_knots(:,i) = Ag*vi;  
end

Atraj = PPTrajectory(foh(ts,A_knots));
Btraj = PPTrajectory(foh(ts,B_knots));
Ctraj = PPTrajectory(foh(ts,C_knots));

for j=1:length(link_indices)
  link_con = struct();
  link_con.link_ndx = link_indices(j);
  link_con.pt = [0;0;0];
  link_con.traj = PPTrajectory(foh(ts,body_knots{j}));
  link_con.dtraj = fnder(link_con.traj);
  link_con.ddtraj = fnder(link_con.traj,2);
  link_constraints(j) = link_con;
end

% for i=0:0.01:ts(end)
%   v.draw(i,xtraj.eval(i));
%   pause(0.01);
% end


% manually extracting the support traj for now
l_foot = r.findLinkInd('l_foot');
r_foot = r.findLinkInd('r_foot');

flight = RigidBodySupportState(r,[]);
l_foot_support = RigidBodySupportState(r,l_foot);
l_toe_support = RigidBodySupportState(r,l_foot,{{'toe'}});
r_foot_support = RigidBodySupportState(r,r_foot);
r_toe_support = RigidBodySupportState(r,r_foot,{{'toe'}});

left_phase = [flight;flight;flight;flight;l_foot_support;l_foot_support;l_foot_support; ...
  l_foot_support;l_foot_support;l_foot_support;l_toe_support; ...
  l_toe_support;l_toe_support;flight;flight;flight];
right_phase = [flight;flight;flight;flight;r_foot_support;r_foot_support;r_foot_support; ...
  r_foot_support;r_foot_support;r_foot_support;r_toe_support; ...
  r_toe_support;r_toe_support;flight;flight;flight];

supports = [left_phase; right_phase; left_phase; right_phase; left_phase; right_phase];

r = r.setInitialState(x_knots(:,1));

% build TV-LQR controller on COM dynamics
x0traj = PPTrajectory(foh(ts,x0_knots));
u0traj = PPTrajectory(foh(ts,u0_knots));
y0traj = PPTrajectory(foh(ts,y0_knots));

Qeps = 1e-8*eye(nq);
Qf = 1e-8*eye(nq);
options.Qy = diag([0.25 0.25 0.25 1 1 1]);
R = 0.001*eye(nq);
% options.tspan = ts;
% options.sqrtmethod = true;
% tv_sys = LinearSystem(Atraj,Btraj,[],[],Ctraj,[]);
% x0traj = x0traj.setOutputFrame(tv_sys.getStateFrame);
% u0traj = u0traj.setOutputFrame(tv_sys.getInputFrame);
% tic;
% [c,V] = tvlqr(tv_sys,x0traj,u0traj,Qeps,R,Qf,options);
% toc
% 
% save('c_and_V.mat','c','V');

load('c_and_V.mat');

ctrl_data = QPControllerData(true,struct(...
  'acceleration_input_frame',AtlasCoordinates(r),...
  'A',Atraj,...
  'B',Btraj,...
  'C',Ctraj,...
  'D',zeros(6,nq),...
  'Qy',options.Qy,...
  'R',R,...
  'S',V.S,...
  's1',V.s1,...
  's2',V.s2,...
  'x0',x0traj,...
  'u0',u0traj,...
  'y0',y0traj,...
  'qtraj',xtraj(1:nq),...
  'support_times',ts,...
  'supports',supports,...
  'mu',1.0,...
  'ignore_terrain',false,...
  'link_constraints',link_constraints,...
  'constrained_dofs',[]));

% instantiate QP controller
options.slack_limit = 100;
options.w_qdd = 0.1*ones(nq,1);
options.w_grf = 0;
options.w_slack = 100;
options.debug = false;
options.use_mex = use_mex;
options.contact_threshold = 0.0005;

boptions.Kp = 250*ones(6,1);
boptions.Kd = 2*sqrt(boptions.Kp);
lfoot_motion = BodyMotionControlBlock(r,'l_foot',ctrl_data,boptions);
rfoot_motion = BodyMotionControlBlock(r,'r_foot',ctrl_data,boptions);
pelvis_motion = BodyMotionControlBlock(r,'pelvis',ctrl_data,boptions);
lhand_motion = BodyMotionControlBlock(r,'l_hand',ctrl_data,boptions);
rhand_motion = BodyMotionControlBlock(r,'r_hand',ctrl_data,boptions);
% boptions.Kp(4:6) = NaN; % don't constrain orientation
% boptions.Kd(4:6) = NaN;
torso_motion = BodyMotionControlBlock(r,'utorso',ctrl_data,boptions);

motion_frames = {lfoot_motion.getOutputFrame,rfoot_motion.getOutputFrame,...
  lhand_motion.getOutputFrame,rhand_motion.getOutputFrame,...
	pelvis_motion.getOutputFrame,torso_motion.getOutputFrame};

options.body_accel_input_weights = [100 100 10 10 100 10];
qp = MomentumQPController(r,motion_frames,ctrl_data,options);

% feedback QP controller with atlas
ins(1).system = 1;
ins(1).input = 2;
ins(2).system = 1;
ins(2).input = 3;
ins(3).system = 1;
ins(3).input = 4;
ins(4).system = 1;
ins(4).input = 5;
ins(5).system = 1;
ins(5).input = 6;
ins(6).system = 1;
ins(6).input = 7;
ins(7).system = 1;
ins(7).input = 8;
ins(8).system = 1;
ins(8).input = 9;
outs(1).system = 2;
outs(1).output = 1;
sys = mimoFeedback(qp,r,[],[],ins,outs);
clear ins;

% feedback foot contact detector with QP/atlas
options.use_lcm=false;
fc = FootContactBlock(r,ctrl_data,options);
ins(1).system = 2;
ins(1).input = 1;
ins(2).system = 2;
ins(2).input = 3;
ins(3).system = 2;
ins(3).input = 4;
ins(4).system = 2;
ins(4).input = 5;
ins(5).system = 2;
ins(5).input = 6;
ins(6).system = 2;
ins(6).input = 7;
ins(7).system = 2;
ins(7).input = 8;
sys = mimoFeedback(fc,sys,[],[],ins,outs);
clear ins;  
  
% feedback PD block
options.use_ik = false;
pd = IKPDBlock(r,ctrl_data,options);
ins(1).system = 1;
ins(1).input = 1;
ins(2).system = 2;
ins(2).input = 2;
ins(3).system = 2;
ins(3).input = 3;
ins(4).system = 2;
ins(4).input = 4;
ins(5).system = 2;
ins(5).input = 5;
ins(6).system = 2;
ins(6).input = 6;
ins(7).system = 2;
ins(7).input = 7;
sys = mimoFeedback(pd,sys,[],[],ins,outs);
clear ins;

% feedback body motion control blocks
ins(1).system = 2;
ins(1).input = 1;
ins(2).system = 2;
ins(2).input = 3;
ins(3).system = 2;
ins(3).input = 4;
ins(4).system = 2;
ins(4).input = 5;
ins(5).system = 2;
ins(5).input = 6;
ins(6).system = 2;
ins(6).input = 7;
sys = mimoFeedback(lfoot_motion,sys,[],[],ins,outs);
clear ins;

ins(1).system = 2;
ins(1).input = 1;
ins(2).system = 2;
ins(2).input = 3;
ins(3).system = 2;
ins(3).input = 4;
ins(4).system = 2;
ins(4).input = 5;
ins(5).system = 2;
ins(5).input = 6;
sys = mimoFeedback(rfoot_motion,sys,[],[],ins,outs);
clear ins;

ins(1).system = 2;
ins(1).input = 1;
ins(2).system = 2;
ins(2).input = 3;
ins(3).system = 2;
ins(3).input = 4;
ins(4).system = 2;
ins(4).input = 5;
sys = mimoFeedback(lhand_motion,sys,[],[],ins,outs);
clear ins;

ins(1).system = 2;
ins(1).input = 1;
ins(2).system = 2;
ins(2).input = 3;
ins(3).system = 2;
ins(3).input = 4;
sys = mimoFeedback(rhand_motion,sys,[],[],ins,outs);
clear ins;

ins(1).system = 2;
ins(1).input = 1;
ins(2).system = 2;
ins(2).input = 3;
sys = mimoFeedback(pelvis_motion,sys,[],[],ins,outs);
clear ins;

ins(1).system = 2;
ins(1).input = 1;
sys = mimoFeedback(torso_motion,sys,[],[],ins,outs);
clear ins;

qt = QTrajEvalBlock(r,ctrl_data);
outs(1).system = 2;
outs(1).output = 1;
sys = mimoFeedback(qt,sys,[],[],[],outs);


S=warning('off','Drake:DrakeSystem:UnsupportedSampleTime');
output_select(1).system=1;
output_select(1).output=1;
sys = mimoCascade(sys,v,[],[],output_select);
warning(S);

traj = simulate(sys,[0 ts(end)],xtraj.eval(0));
playback(v,traj,struct('slider',true));

end
