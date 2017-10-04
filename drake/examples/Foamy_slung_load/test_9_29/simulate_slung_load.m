%r = RigidBodyManipulator('slung_load_foamy.URDF');
p = foamy_slung_load();
v = p.constructVisualizer();
% v.inspector();


traj = simulate(p,[0 1]);

t2=traj(1:9);
t2 = t2.setOutputFrame(v.getInputFrame);

playback(v,t2,struct('slider',true));