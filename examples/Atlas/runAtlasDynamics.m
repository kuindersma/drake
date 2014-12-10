function runAtlasDynamics
% Simulate the (passive) dynamics of the Atlas model 

% Load the model with a floating base
options.floating = true;
options.dt = 0.005;
options.terrain = RigidBodyFlatTerrain;
options.ignore_self_collisions = true;
r = Atlas('urdf/atlas_minimal_contact.urdf',options);
r = r.removeCollisionGroupsExcept({'heel','toe','back','front','knee','butt'});
r = compile(r);

% Initialize the viewer
v = r.constructVisualizer;
v.display_dt = 0.02;

% Compute a feasible set of initial conditions for the simulation (e.g. no
% penetration)
x0 = Point(r.getStateFrame);
x0(5) = pi/2;
x0(3) = 0.1;
x0 = resolveConstraints(r,x0);
% 
% % Forward simulate dynamics with visulazation, then playback at realtime
% S=warning('off','Drake:DrakeSystem:UnsupportedSampleTime');
% output_select(1).system=1;
% output_select(1).output=1;
% sys = mimoCascade(r,v,[],[],output_select);
% warning(S);
% traj = simulate(sys,[0 2],x0);
% playback(v,traj,struct('slider',true));

v.drawWrapper(0,x0);
tic;
xtraj = r.simulate([0 2],x0); 
toc;
v.playback(xtraj,struct('slider',true));

end
