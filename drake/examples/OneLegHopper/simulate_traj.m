function simulate_traj()

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
dt = 0.01;
r = TimeSteppingRigidBodyManipulator(s,dt,options);
r = r.setStateFrame(OneLegHopperState(r));
r = r.setOutputFrame(OneLegHopperState(r));
v = r.constructVisualizer();

load hopper_iLQR_traj;
utraj = utraj.setOutputFrame(r.getInputFrame);
sys = cascade(utraj,r);
xtraj_sim = simulate(sys,[0,xtraj.tspan(2)],xtraj.eval(0));
v.playback(xtraj_sim,struct('slider',true));

end