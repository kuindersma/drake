function fallingBrickLCP(planar)

dt = 0.001;
options.floating = true;
options.terrain = RigidBodyFlatTerrain();

if planar
  options.twoD = true;
  options.view = 'right';
  x0 = [0;1.5;2;1;0;0];
  p = TimeSteppingRigidBodyManipulator('PlanarFallingBrickContactPoints.urdf',dt,options);
else
  x0 = [0;1;randn(10,1)];
  p = TimeSteppingRigidBodyManipulator('FallingBrickContactPoints.urdf',dt,options);
end

x0 = p.resolveConstraints(x0);

v = p.constructVisualizer();
v.display_dt = dt*2;
if 1 
  sys = cascade(p,v);
  sys.simulate([0 3],x0);
else
  v.drawWrapper(0,x0);
  xtraj = p.simulate([0 3],x0); 
  v.playback(xtraj,struct('slider',true));
end

