function runAtlasSagittalDynamics
% just runs it as a passive system

options.twoD = true;
options.view = 'right';
options.floating = true;
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;

s = 'urdf/atlas_minimal_contact.urdf';
dt = 0.005;
w = warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits');
r = TimeSteppingRigidBodyManipulator(s,dt,options);
% r = r.removeCollisionGroupsExcept({'heel','toe'});
r = compile(r);
warning(w);

v = r.constructVisualizer;
v.display_dt = 0.02;

% Run simulation, then play it back at realtime speed
tic;
xtraj = simulate(r,[0 3]);
toc;

v.playback(xtraj,struct('slider',true));
