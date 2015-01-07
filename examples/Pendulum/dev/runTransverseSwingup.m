function ctrans=runTransverseSwingup

p = PendulumPlant();

[utraj,xtraj] = swingUpTrajectory(p);

p=setInputLimits(p,-inf,inf);

figure(1); clf; hold on;
fnplt(xtraj);

w=randn(2,1);
Nsteps = 50;

% Use orthogonal initial and final surface normals (not essential)
fs0 = p.dynamics(0,xtraj.eval(0),utraj.eval(0));
fsend = p.dynamics(utraj.tspan(end),xtraj.eval(utraj.tspan(end)),utraj.eval(utraj.tspan(end)));
init_surf_normal = fs0/norm(fs0);
final_surf_normal = fsend/norm(fsend);

disp('Designing transversal surfaces...')
transSurf = TransversalSurface.design(p,w, xtraj, utraj, Nsteps,1,init_surf_normal, final_surf_normal);

plotSurface(transSurf,xtraj,.4); drawnow;

disp('Computing transverse lqr...')
Q = 10*eye(2);
R = 1;
Vf = Q;
[ctrans,Vtraj] = transverseLQR(p,xtraj,utraj,Q,R,Vf,transSurf);

sys = feedback(p,ctrans);

Pi = getPi(transSurf,xtraj.tspan(end));
Vf = Pi*Vf*Pi';

% disp('Estimating funnel...')
% nx = p.getNumStates();
% nu = p.getNumInputs();
% xtraj_cl = MixedTrajectory({xtraj, PPTrajectory(foh(xtraj.tspan,xtraj.tspan))},{1:nx,nx+(1:nu)}); % add controller state
% psys = taylorApprox(sys,xtraj_cl,[],3,3);  %ignore var 3 (tau)
% options=struct();
% options.rho0_tau=10;
% options.max_iterations=25;
% V=sampledTransverseVerification(psys,Vf,Vtraj,Vtraj.getBreaks(),xtraj,utraj2,transSurf,options);
% figure(3); clf;
% transSurf.plotFunnel(V,xtraj);
% fnplt(xtraj);

disp('Simulating...');
v=PendulumVisualizer;

%keyboard; 

for i=1:5
%  x0=[randn(2,1);0];
  y=simulate(sys,[0 5]);%,x0);
  figure(1);
  fnplt(y);
  figure(25);
  v.playback(y);
end


end


