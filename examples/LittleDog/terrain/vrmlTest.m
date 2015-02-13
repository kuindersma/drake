function vrmlTest

options.floating = true;
load terrainB;
options.terrain = RigidBodyHeightMapTerrain(x,y,Z);
r = RigidBodyManipulator('../drake/systems/plants/test/FallingBrick.urdf',options);

v = r.constructVisualizer();
v.draw(0,[0;0;10;zeros(9,1)]);

opt2.viewer = 'RigidBodyWRLVisualizer';
v2 = r.constructVisualizer(opt2);
v2.draw(0,[0;0;10;zeros(9,1)]);

