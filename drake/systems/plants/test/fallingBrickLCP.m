function fallingBrickLCP

options.floating = true;
options.twoD = true;
options.floating = true;
options.terrain = RigidBodyFlatTerrain();
% options.ignore_self_collisions = true;
options.use_bullet = false;
s = 'FallingBrickContactPoints.urdf';
% s = 'FallingBrickBetterCollisionGeometry.urdf';
p = TimeSteppingRigidBodyManipulator(s,.01,options);
% p = p.addRobotFromURDF(s,[],[],options);
% x0 = [0;1;2;1.57;0;0;10;zeros(5,1)];
x0 = [0; 2; 1.57; 10; zeros(2,1)];

% Forward simulate dynamics with visulazation, then playback at realtime
v=p.constructVisualizer();
S=warning('off','Drake:DrakeSystem:UnsupportedSampleTime');
output_select(1).system=1;
output_select(1).output=1;
sys = mimoCascade(p,v,[],[],output_select);
warning(S);
traj = simulate(sys,[0 2],x0);
playback(v,traj,struct('slider',true));