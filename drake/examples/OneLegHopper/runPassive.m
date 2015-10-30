warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits');
options.twoD = true;
options.view = 'right';
options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.ignore_self_collisions = true;
options.use_bullet = false;
options.enable_fastqp = false;
s = 'OneLegHopper.urdf';
dt = 0.001;
r = TimeSteppingRigidBodyManipulator(s,dt,options);
r = r.setStateFrame(OneLegHopperState(r));
r = r.setOutputFrame(OneLegHopperState(r));
v=r.constructVisualizer();

traj = simulate(r,[0 2]);
playback(v,traj,struct('slider',true));
