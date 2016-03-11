function testUpdateWithBVForces

% silence some warnings
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints')
warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits')

options.floating = true;
options.dt = 0.002;
options.terrain = RigidBodyFlatTerrain();
options.use_bullet = false;

urdf = [getDrakePath,'/examples/Atlas/urdf/atlas_minimal_contact.urdf'];
r = TimeSteppingRigidBodyManipulator(urdf,options.dt,options);
r = r.removeCollisionGroupsExcept({'heel','toe'});
r = compile(r);

nx = r.getNumStates;
nu = r.getNumInputs;
nz = r.getNumContactPairs*4;

for i=1:100

  x = randn(nx,1);
  u = randn(nu,1);
  z = rand(nz,1);

  [xdn,df] = r.updateWithBasisForces(0,x,u,z);
  [xdn2,df2] = geval(@r.updateWithBasisForces,0,x,u,z,struct('grad_method','numerical'));

  tol=1e-4;
  valuecheck(xdn,xdn2);
  valuecheck(df,df2,tol);
end