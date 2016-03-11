function runDDP

warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits');
options.twoD = true;
options.view = 'right';
options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.ignore_self_collisions = true;
options.use_bullet = false;
options.enable_fastqp = false;
s = 'OneLegHopper.urdf';
dt = 0.005;
r = TimeSteppingRigidBodyManipulator(s,dt,options);
r = r.setStateFrame(OneLegHopperState(r));
r = r.setOutputFrame(OneLegHopperState(r));

nx = r.getNumStates;
nu = r.getNumInputs;
nc = r.getNumContactPairs;
nz = nc*2;

v = r.constructVisualizer();
options.visualizer = v;
options.enable_visualizer = true;
options.visualizer_period = 5;
options.enable_line_search = false;

xG =  zeros(nx,1); 
xG(1) = 0.0;
xG(2) = 1.0;


q0 = [0;0;.6;-1.2;.6+pi/2];
x0 = [q0;0*q0];
x0 = double(r.resolveConstraints(x0));
ddp = DifferentialDynamicProgrammingWithContacts(r,false,options);
ddp = ddp.addRunningCost(@cost);
ddp = ddp.addFinalCost(@final_cost);

T = 1.0;
N = T/dt;

xtraj = repmat(x0,1,N);
utraj = zeros(nu,N);
ztraj = zeros(nz,N);
[xtraj,utraj,ztraj] = ddp.solveTraj(xtraj,utraj,ztraj,N,dt,500);

xtraj = ddp.reconstructStateTrajectory(xtraj,0,dt);
utraj = ddp.reconstructInputTrajectory(utraj,0,dt);

if (0) % open-loop playback
  sys = cascade(utraj,r);
  xtraj = simulate(sys,[0,T],xtraj.eval(0));
end
v.playback(xtraj,struct('slider',false));



  function [g,dg,d2g] = cost(x,u,z)
    Q = diag([100*ones(nx/2,1); 1.0*ones(nx/2,1)]);
    Ru = 0.01*eye(nu);
    Rz = 0.01*eye(nz);

    g = (x-xG)'*Q*(x-xG) + u'*Ru*u + z'*Rz*z;
    if nargout > 1
      dg = [2*(x'*Q -xG'*Q), 2*u'*Ru, 2*z'*Rz];
      d2g = [2*Q, zeros(nx,nu+nz); 
             zeros(nu,nx), 2*Ru, zeros(nu,nz); 
             zeros(nz,nx+nu),2*Rz];
    end
  end


  function [g,dg,d2g] = final_cost(x)
    Q = diag([100*ones(nx/2,1); 1.0*ones(nx/2,1)]);

    g = (x-xG)'*Q*(x-xG);
    if nargout > 1
      dg = 2*(x'*Q - xG'*Q);
      d2g = 2*Q;
    end
  end

end