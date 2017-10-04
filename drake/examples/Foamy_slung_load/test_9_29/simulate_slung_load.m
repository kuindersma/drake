%r = RigidBodyManipulator('slung_load_foamy.URDF');
p = foamy_slung_load();
v = r.constructVisualizer();
v.inspector();