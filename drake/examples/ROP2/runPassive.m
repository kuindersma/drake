function runPassive

options.terrain = RigidBodyFlatTerrain();
options.floating = true;
r = TimeSteppingRigidBodyManipulator('urdf/robotis_op.urdf',0.005,options);

v = r.constructVisualizer;
v.display_dt = .02;

% xtraj = simulate(r,[0 1],zeros(r.getNumStates,1));
% 
% v.playback(xtraj,struct('slider',true));
