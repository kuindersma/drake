function robustSwingUp(N,D)

p = PendulumPlant;
v = PendulumVisualizer();

K = [-100,-30]; % PD gains

if 1
 [utraj,xtraj,z,prog] = robustSwingUpTrajectory(p,N,D);
    xs = z(prog.x_inds);
    dxs = z(prog.dx_inds);
    us = z(prog.u_inds);
    dus = z(prog.du_inds);
    hs = z(prog.h_inds);

%     xs_ = xs;
%     xs__ = xs;
%     for i=1:N-1
%       xs_(:,i+1) = xs_(:,i) + hs(i) * p.dynamics_w(0,xs_(:,i), us(i) + dus(i), lb);
%       xs__(:,i+1) = xs__(:,i) + hs(i) * p.dynamics_w(0,xs__(:,i), us(i) + K*(xs__(:,i)-xs(:,i)), lb);
%     end

    figure(12);
    plot(xs(1,:),xs(2,:),'b.-','MarkerSize',10);
    hold on;
    plot(xs(1,:)+dxs(1,:),xs(2,:)+dxs(2,:),'g.-','MarkerSize',10);
%     plot(xs_(1,:),xs_(2,:),'r.-','MarkerSize',10);
%       plot(xs__(1,:),xs__(2,:),'m.-','MarkerSize',10);
    hold off;

    keyboard
elseif 1
 
  [utraj,xtraj] = swingUpTrajectory(p,N);


else
  load nominal_20.mat

  xtraj = xtraj.setOutputFrame(p.getStateFrame);
  utraj = utraj.setOutputFrame(p.getInputFrame);
end

% v.playback(xtraj,struct('slider',true));
% keyboard;

% breaks = utraj.getBreaks();
% hs = breaks(2:end)-breaks(1:end-1);
% figure(99);
% lte = [];
% for i=1:N-1
%   lte = [lte, p.integrationCost(hs(i),xtraj.eval(breaks(i)),utraj.eval(breaks(i)))];
% end
% plot(cumsum(hs),lte,'b-');

if 0
  % run LQR on nominal plan
  Q = diag([100 10]);
  R = 0.01;
  Qf = 2*Q;
  c = tvlqr(p,xtraj,utraj,Q,R,Qf);
else
  % use fixed linear PD
  nX = p.getNumStates;
  nU = p.getNumInputs;
  % create linear feedback controller
  iframe = CoordinateFrame([p.getStateFrame.name,' - x0(t)'],nX,p.getStateFrame.prefix);
  p.getStateFrame.addTransform(AffineTransform(p.getStateFrame,iframe,eye(nX),-xtraj));
  iframe.addTransform(AffineTransform(iframe,p.getStateFrame,eye(nX),xtraj));

  oframe = CoordinateFrame([p.getInputFrame.name,' + u0(t)'],nU,p.getInputFrame.prefix);
  oframe.addTransform(AffineTransform(oframe,p.getInputFrame,eye(nU),utraj));
  p.getInputFrame.addTransform(AffineTransform(p.getInputFrame,oframe,eye(nU),-utraj));

  D = ConstantTrajectory(K);

  c = AffineSystem([],[],[],[],[],[],[],D,[]);
  c = setInputFrame(c,iframe);
  c = setOutputFrame(c,oframe);
end

count = 0;
fail_count = 0;
for i=1:1
%   p = p.setMass(1+rand*(ub(1)-lb(1))+lb(1));
%   p = p.setCOMLength(0.5+rand*(ub(2)-lb(2))+lb(2));
%   p = p.setDamping(0.1+rand*(ub(3)-lb(3))+lb(3));
  p = p.setMass(1+rand*(ub-lb)+lb);
  sys = feedback(p,c);

  if 1
    % Forward simulate dynamics with visulazation, then playback at realtime
    S=warning('off','Drake:DrakeSystem:UnsupportedSampleTime');
    output_select(1).system=1;
    output_select(1).output=1;
    sys = mimoCascade(sys,v,[],[],output_select);
    warning(S);
  end
  traj=simulate(sys,[0,utraj.tspan(2)],xtraj.eval(xtraj.tspan(1)));
  v.playback(traj,struct('slider',true));

  count = count + 1;
  if norm(double(traj.eval(traj.tspan(2)))-double(p.xG)) > 0.1
    fail_count = fail_count + 1;
  end
  fprintf('Successes: %d, Failures: %d\n',count-fail_count,fail_count);

end

keyboard
end
