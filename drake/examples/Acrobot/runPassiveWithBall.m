function runPassiveWithBall

% options.use_bullet = true;
options.view = 'right';
% options.ignore_self_collisions = false;
options.terrain = RigidBodyFlatTerrain();


r = TimeSteppingRigidBodyManipulator(PlanarRigidBodyManipulator('AcrobotCollision.urdf'),.001,options);
options.floating = true;
r = r.addRobotFromURDF('../../systems/plants/test/ball.urdf',zeros(3,1),zeros(3,1), options);

v = r.constructVisualizer;

%acrobot_state = [0;-pi/2;0;-1];
acrobot_pos = [0;0];
acrobot_vel = -[.5;5];
ball_pos = [1.5;-.25;0];
ball_vel = [0;-1;0];

x0 = [acrobot_pos; ball_pos; acrobot_vel; ball_vel];

% x0 = .2*randn(r.getNumStates,1);
% x0 = r.resolveConstraints(x0);
% v.drawWrapper(0,x0);

% v.draw(0,x0);

ytraj = r.simulate([0 1],x0);
% ytraj = r.simulate([0 1]);
v.playback(ytraj,struct('slider',true));

