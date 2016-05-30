function testDTSwingUP

load robust_200_linear_0.4.mat
% load nominal_200_linear.mat

% h = 0.01;
ts = utraj.getBreaks();
h=ts(2)-ts(1);

p = PendulumPlantDT(h);
v = PendulumVisualizer();

% K = [-100,-30]; % PD gains

xtraj = xtraj.setOutputFrame(p.getStateFrame);
utraj = utraj.setOutputFrame(p.getInputFrame);

% v.playback(xtraj,struct('slider',true));
% keyboard;

% % run LQR on nominal plan
% Q = diag([100 10]);
% R = 0.01;
% Qf = 2*Q;
% c = tvlqr(p,xtraj,utraj,Q,R,Qf);
% sys = feedback(p,c);
% traj=simulate(sys,utraj.tspan,xtraj.eval(xtraj.tspan(1)));
% v.playback(traj,struct('slider',true));

% keyboard;

nX = p.getNumStates;
nU = p.getNumInputs;
% create linear feedback controller
iframe = CoordinateFrame([p.getStateFrame.name,' - x0(t)'],nX,p.getStateFrame.prefix);
p.getStateFrame.addTransform(AffineTransform(p.getStateFrame,iframe,eye(nX),-xtraj));
iframe.addTransform(AffineTransform(iframe,p.getStateFrame,eye(nX),xtraj));

oframe = CoordinateFrame([p.getInputFrame.name,' + u0(t)'],nU,p.getInputFrame.prefix);
oframe.addTransform(AffineTransform(oframe,p.getInputFrame,eye(nU),utraj));
p.getInputFrame.addTransform(AffineTransform(p.getInputFrame,oframe,eye(nU),-utraj));

c = LinearSystem([],[],[],[],[],K);
c = setInputFrame(c,iframe);
c = setOutputFrame(c,oframe);

count = 0;
fail_count = 0;
for i=1:1
%   p = p.setMass(1+rand*(ub(1)-lb(1))+lb(1));
%   p = p.setCOMLength(0.5+rand*(ub(2)-lb(2))+lb(2));
%   p = p.setDamping(0.1+rand*(ub(3)-lb(3))+lb(3));
  p = p.setMass(1+0.5);
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
