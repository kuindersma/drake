classdef AcrobotPlant < Manipulator 
  
  properties
    % parameters from Spong95 (except inertias are now relative to the
    % joints)
    % axis)
    l1 = 1; l2 = 2;  
    m1 = 1; m2 = 1;  
    g = 9.81;
    b1=.1;  b2=.1;
%    b1=0; b2=0;
    lc1 = .5; lc2 = 1; 
    Ic1 = .083;  Ic2 = .33;
    
    xG
    uG
    w_max
  end
  
  methods
    function obj = AcrobotPlant(w_max)
      if nargin < 1
        w_max = 0;
      end
      obj = obj@Manipulator(2,1);
      obj = setInputLimits(obj,-10,10);

      obj = setInputFrame(obj,CoordinateFrame('AcrobotInput',1,'u',{'tau'}));
      obj = setStateFrame(obj,CoordinateFrame('AcrobotState',4,'x',{'theta1','theta2','theta1dot','theta2dot'}));
      obj = setOutputFrame(obj,obj.getStateFrame);
      
      obj.xG = Point(obj.getStateFrame,[pi;0;0;0]);
      obj.uG = Point(obj.getInputFrame,0);
      
      obj = setParamFrame(obj,CoordinateFrame('AcrobotParams',6,'p',...
        { 'b1','b2','lc1','lc2','Ic1','Ic2' }));
%      obj = setParamFrame(obj,CoordinateFrame('AcrobotParams',10,'p',...
%        { 'l1','l2','m1','m2','b1','b2','lc1','lc2','I1','I2' }));
%      obj = setParamFrame(obj,CoordinateFrame('AcrobotParams',1,'p',...
%        { 'lc2' }));
      obj = setParamLimits(obj,zeros(obj.getParamFrame.dim,1));
      obj.w_max = w_max;
    end
    
    function obj = addDisturbanceBound(obj,dist)
      obj.w_max = dist;
    end

    function [H,C,B] = manipulatorDynamics(obj,q,qd)
      % keep it readable:
      m1=obj.m1; m2=obj.m2; l1=obj.l1; g=obj.g; lc1=obj.lc1; lc2=obj.lc2; b1=obj.b1; b2=obj.b2;
      I1 = obj.Ic1 + obj.m1*obj.lc1^2; I2 = obj.Ic2 + obj.m2*obj.lc2^2;
      m2l1lc2 = m2*l1*lc2;  % occurs often!

      c = cos(q(1:2,:));  s = sin(q(1:2,:));  s12 = sin(q(1,:)+q(2,:));
      
      h12 = I2 + m2l1lc2*c(2);
      H = [ I1 + I2 + m2*l1^2 + 2*m2l1lc2*c(2), h12; h12, I2 ];
      
      C = [ -2*m2l1lc2*s(2)*qd(2), -m2l1lc2*s(2)*qd(2); m2l1lc2*s(2)*qd(1), 0 ];
      G = g*[ m1*lc1*s(1) + m2*(l1*s(1)+lc2*s12); m2*lc2*s12 ];
            
      % accumate total C and add a damping term:
      C = C*qd + G + [b1;b2].*qd;

      B = [0; 1];
    end
    
    function [T,U] = energy(obj,x)
      m1=obj.m1; m2=obj.m2; l1=obj.l1; g=obj.g; lc1=obj.lc1; lc2=obj.lc2; b1=obj.b1; b2=obj.b2;
      I1 = obj.Ic1 + obj.m1*obj.lc1^2; I2 = obj.Ic2 + obj.m2*obj.lc2^2;
      q = x(1:2); qd = x(3:4);
      c = cos(q(1:2,:));  s = sin(q(1:2,:));  c12 = cos(q(1,:)+q(2,:));
      
      T = .5*I1*qd(1)^2 + .5*(m2*l1^2 + I2 + 2*m2*l1*lc2*c(2))*qd(1)^2 + .5*I2*qd(2)^2 + (I2 + m2*l1*lc2*c(2))*qd(1)*qd(2);
      U = -m1*g*lc1*c(1) - m2*g*(l1*c(1)+lc2*c12);
    end
    
    % todo: also implement sodynamics here so that I can keep the
    % vectorized version?
    
    function [f,df,d2f,d3f] = dynamics(obj,t,x,u)
      
      w = rand*(2*obj.w_max) - obj.w_max;
      f = dynamics@Manipulator(obj,t,x,u+w);
      if (nargout>1)
        [df,d2f,d3f]= dynamicsGradients(obj,t,x,u,nargout-1);
      end
    end

    function [f,df] = dynamics_w(obj,t,x,u,w)
      if (nargout>1)
        [f,df] = dynamics(obj,t,x,u+w);
        nx = obj.getNumStates;
        nu = obj.getNumInputs;
        df = [df,df(:,1+nx+(1:nu))];
      else
        f = dynamics(obj,t,x,u+w);
      end
    end

    function x = getInitialState(obj)
      x = .1*randn(4,1);
    end
    
    function n = getNumPositions(obj)
      n = 2;
    end
    
    function n = getNumVelocities(obj)
      n = 2;
    end
    
    function n = getNumDisturbances(obj)
      n =1;
    end
    
    function [c,V]=balanceLQR(obj)
      Q = diag([10,10,1,1]); R = 1;
      if (nargout<2)
        c = tilqr(obj,obj.xG,obj.uG,Q,R);
      else
        if any(~isinf([obj.umin;obj.umax]))
          error('currently, you must disable input limits to estimate the ROA');
        end
        [c,V] = tilqr(obj,obj.xG,obj.uG,Q,R);
        pp = feedback(obj.taylorApprox(0,obj.xG,obj.uG,3),c);
        options.method='levelSet';
        V=regionOfAttraction(pp,V,options);
      end
    end
    
    
    
    
    
    function [utraj,xtraj]=robustSwingUpTrajectory(obj,N,M,disturbances)
      x0 = zeros(4,1); 
      xf = double(obj.xG);
      tf0 = 4;
      
      d = linspace(-disturbances,disturbances,M);
      options.integration_method = DirtranTrajectoryOptimization.MIDPOINT;
      prog = RobustDirtranTrajectoryOptimization(obj,N,M,[4 6],options);
      disp('constructor done');
      prog = prog.setDisturbances(d);
      prog = prog.addStateConstraint(ConstantConstraint(x0),1);
      prog = prog.addStateConstraint(ConstantConstraint(xf),N);
      prog = prog.addRunningCost(@running_cost);
      prog = prog.addFinalCost(@final_cost);
      prog = prog.addRobustConstraints(@robust_cost);
      
      nx=4; nu=1; nw=1;
%       for j=1:10
%         hr = randn();
%         xr = randn(nx,1);
%         xr2 = randn(nx,1);
%         ur = randn(nu,1);
%         ur2 = randn(nu,1);
%         wr = randn(nw,1);
%         gr = randn();
%         zi = abs(randn(M,1));
%         zij = abs(randn());
%       
%         [~,df1] = geval(@cost,hr,xr,ur,struct('grad_method','taylorvar'));
%         [~,df2] = cost(hr,xr,ur);
%       
%         valuecheck(df1,df2,1e-4);
%       
%         [~,df1] = geval(@robust_cost,hr,xr,ur,wr,struct('grad_method','taylorvar'));
%         [~,df2] = robust_cost(hr,xr,ur,wr);
%       
%         valuecheck(df1,df2,1e-4);
%       
%         [~,df1] = geval(@finalCost,0,xr,struct('grad_method','taylorvar'));
%         [~,df2] = finalCost(0,xr);
%       
%         valuecheck(df1,df2,1e-4);
%         
%         
%         [~,df1] = geval(@prog.forward_constraint_fun,hr,xr,xr2,ur,wr,struct('grad_method','taylorvar'));
%         [~,df2] = prog.forward_constraint_fun(hr,xr,xr2,ur,wr);
%       
%         valuecheck(df1,df2,1e-4);
%         
% 
%         [~,df1] = geval(@prog.backward_constraint_fun,hr,xr,xr2,ur,wr,struct('grad_method','taylorvar'));
%         [~,df2] = prog.backward_constraint_fun(hr,xr,xr2,ur,wr);
%       
%         valuecheck(df1,df2,1e-4);
%         
%         
%         tmp1 = @(gamma,x0,x1,u0,u1,w1) prog.robust_constraint_fun(@cost,@robust_cost,1,gamma,x0,x1,u0,u1,w1);
% 
%         [~,df1] = geval(tmp1,gr,xr,xr2,ur,ur2,wr,struct('grad_method','taylorvar'));
%         [~,df2] = tmp1(gr,xr,xr2,ur,ur2,wr);
%       
%         valuecheck(df1,df2,1e-4);
%         
%         tmp2 = @(gamma,x0,x1,u0,u1,w1,z) prog.complementarity_fun(@cost,@robust_cost,1,gamma,x0,x1,u0,u1,w1,z);
% 
%         [~,df1] = geval(tmp2,gr,xr,xr2,ur,ur2,wr,zij,struct('grad_method','taylorvar'));
%         [~,df2] = tmp2(gr,xr,xr2,ur,ur2,wr,zij);
%       
%         valuecheck(df1,df2,1e-4);
% 
%         [~,df1] = geval(@prog.running_gamma_cost,gr,struct('grad_method','taylorvar'));
%         [~,df2] = prog.running_gamma_cost(gr);
%       
%         valuecheck(df1,df2,1e-4);
%         
%         [~,df1] = geval(@prog.z_sum_constr,zi,struct('grad_method','taylorvar'));
%         [~,df2] = prog.z_sum_constr(zi);
%       
%         valuecheck(df1,df2,1e-4);
% 
%         [~,df1] = geval(@prog.w_equality,wr,zi,struct('grad_method','taylorvar'));
%         [~,df2] = prog.w_equality(wr,zi);
%       
%         valuecheck(df1,df2,1e-4);
%         
%         
%       end
% 
%       keyboard



      
      % add a display function to draw the trajectory on every iteration
      function displayStateTrajectory(t,x,u)
        plot(x(1,:),x(3,:),'b.-','MarkerSize',10);
        drawnow;
      end
      prog = addTrajectoryDisplayFunction(prog,@displayStateTrajectory);
   

     
      traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));
      
      for attempts=1:1
        disp('Running solve');
        tic
        [xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
        toc
        if info==1, break; end
      end
      keyboard
    
      
      function [g,dg] = running_cost(dt,x,u)
        R = 10;
        g = sum((R*u).*u,1);
        dg = [zeros(1,1+nx),2*u'*R];
        return;
        
        xd = repmat([pi;0;0;0],1,size(x,2));
        xerr = x-xd;
        xerr(1,:) = mod(xerr(1,:)+pi,2*pi)-pi;
        
        Q = diag([10,10,1,1]);
        R = 100;
        g = sum((Q*xerr).*xerr + (R*u).*u,1);
        
        if (nargout>1)
          dgddt = 0;
          dgdx = 2*xerr'*Q;
          dgdu = 2*u'*R;
          dg = [dgddt,dgdx,dgdu];
        end
      end
      
      function [g,dg] = robust_cost(dx,du,w)
        W = 0*eye(length(w));
        Qw = 500*eye(4);
        Rw = 0.1;
        g = dx'*Qw*dx + du'*Rw*du + w'*W*w;
        dg = [2*dx'*Qw 2*du'*Rw, 2*w'*W];
      end
      
      function [h,dh] = final_cost(t,x)
        h = t;
        dh = [1,zeros(1,nx)];
      end
     
    end
    
    
    function [utraj,xtraj]=swingUpTrajectory(obj,N)
      x0 = zeros(4,1); 
      xf = double(obj.xG);
      tf0 = 4;
      
      prog = DirtranTrajectoryOptimization(obj,N,[2 6]);
      prog = prog.addStateConstraint(ConstantConstraint(x0),1);
      prog = prog.addStateConstraint(ConstantConstraint(xf),N);
      prog = prog.addRunningCost(@cost);
      prog = prog.addFinalCost(@finalCost);
      
      traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));
      
      for attempts=1:10
        tic
        [xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
        toc
        if info==1, break; end
      end
      
      function [g,dg] = cost(dt,x,u)
        R = 1;
        g = sum((R*u).*u,1);
        dg = [zeros(1,1+size(x,1)),2*u'*R];
        return;
        
        xd = repmat([pi;0;0;0],1,size(x,2));
        xerr = x-xd;
        xerr(1,:) = mod(xerr(1,:)+pi,2*pi)-pi;
        
        Q = diag([10,10,1,1]);
        R = 100;
        g = sum((Q*xerr).*xerr + (R*u).*u,1);
        
        if (nargout>1)
          dgddt = 0;
          dgdx = 2*xerr'*Q;
          dgdu = 2*u'*R;
          dg = [dgddt,dgdx,dgdu];
        end
      end
      
      function [h,dh] = finalCost(t,x)
        h = t;
        dh = [1,zeros(1,size(x,1))];
        return;
        
        xd = repmat([pi;0;0;0],1,size(x,2));
        xerr = x-xd;
        xerr(1,:) = mod(xerr(1,:)+pi,2*pi)-pi;
        
        Qf = 100*diag([10,10,1,1]);
        h = sum((Qf*xerr).*xerr,1);
        
        if (nargout>1)
          dh = [0, 2*xerr'*Qf];
        end
      end      

    end

    
  end
  
end
