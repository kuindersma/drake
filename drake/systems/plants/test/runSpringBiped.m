options.enable_fastqp = false;
options.use_bullet = false;
options.floating=true;
options.terrain = RigidBodyFlatTerrain();
r = TimeSteppingRigidBodyManipulator('MassSpringBiped.urdf', 0.005,options);

v=r.constructVisualizer;

nx = r.getNumStates;
x0 = 1e-3*randn(nx,1);

x0(3)=2.0;

xtraj = simulate(r,[0 4],x0);

v.playback(xtraj,struct('slider',true));