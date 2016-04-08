function runAtlasSagittalDynamics
% just runs it as a passive system

options.twoD = true;
options.view = 'right';
options.floating = true;
options.use_bullet = true;
options.terrain = RigidBodyFlatTerrain();
s = 'urdf/atlas_minimal_contact.urdf';
dt = 0.001;
w = warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits');
r = TimeSteppingRigidBodyManipulator(s,dt,options);
% r = r.removeCollisionGroupsExcept({'heel','toe','back','front','knee','butt'});
% r = compile(r);
warning(w);

v = r.constructVisualizer;
v.display_dt = 0.01;

% Forward simulate dynamics with visulazation, then playback at realtime
S=warning('off','Drake:DrakeSystem:UnsupportedSampleTime');
output_select(1).system=1;
output_select(1).output=1;
sys = mimoCascade(r,v,[],[],output_select);
warning(S);
traj = simulate(sys,[0 2]);
playback(v,traj,struct('slider',true));