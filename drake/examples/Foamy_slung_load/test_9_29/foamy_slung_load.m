classdef foamy_slung_load < DrakeSystem

  properties

    plane_rbm
    plane_link_id=2
  end
  
  methods
  
    function obj = foamy_slung_load(options)
   
      obj = obj@DrakeSystem(17,0,4,17,0,1);

      obj = obj.setStateFrame(CoordinateFrame('FoamyState',17,'x',{'x1','x2','x3','q0','q1','q2','q3','theta','phi','v1','v2','v3','w1','w2','w3','theta_dot','phi_dot'}));
      obj = obj.setInputFrame(CoordinateFrame('FoamyInput',4,'u',{'thr','ail','elev','rud'}));
      obj = obj.setOutputFrame(obj.getStateFrame);       

      options.floating = 'quat';
      urdf = 'slung_load_foamy.URDF';

      warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits');
      warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
      warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
      warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits');
      rbm = RigidBodyManipulator(urdf,options);
%       rbm = rbm.setStateFrame(obj.getStateFrame);
%       rbm = rbm.setOutputFrame(obj.getStateFrame);  

      obj.plane_rbm = rbm;

    end
    
    function nw = getNumDisturbances(obj)
      switch obj.disturbance_type
        case 1
          nw = 3;
        case 2
          nw = obj.getNumContStates();
        case 3
          nw = obj.getNumInputs();
        otherwise
          error('Unknown disturbance type');
      end
    end
    
    function [f,df] = dynamics(obj,t,x,u)
        quat = x(4:7); %Quaternion rotation from body to lab frame
        float_v = x(10:12); %Lab-frame velocity vector
        float_w = x(13:15); %Body-frame angular velocity

        v_plant_to_manip = [4 5 6 1 2 3 7 8];
        q_manip = x(1:9);
        v_manip = x(9+v_plant_to_manip);
        [H,C,~] = obj.plane_rbm.manipulatorDynamics(q_manip,v_manip);
        [~,F,T] = foamy_dynamics(t,x,u);
        kinsol = doKinematics(obj.plane_rbm, q_manip, [], struct('compute_gradients', true));
  
        body = obj.plane_rbm.findLinkId('airplane_body');
        base = 1;
        J = geometricJacobian(obj.plane_rbm, kinsol, base, body , base);
        
%         F_gen = J'*[F;T]; 
        F_gen = J'*[0;0;0;0;0;10]; % TODO: haven't had a chance to debug this 
        qdd = H\(-C+[F_gen;0;0]); %[wx,wy,wz,vx,vy,vz,theta_dot,phi_dot]'
        f = [float_v;.5*[-quat(2:4)'*float_w; quat(1)*float_w + cross(quat(2:4),float_w)];v_manip(7:8);qdd(4:6);qdd(1:3);qdd(7:8)];
        df = 0;
    end
    
    function y = output(obj,t,x,u)
      y=x;
    end
        
    function [xtraj,utraj] = runDircol(obj,display)
            % initial conditions:
            disp('slung_initial')
            [x0, u0] = findTrim(obj,6) %find trim conditions for level flight at 6 m/s
            x0(1) = -5;
            x0(3) = 1.5;

            % final conditions:
            xf = x0;
            xf(1) = 5; %translated in x

            tf0 = (xf(1)-x0(1))/6; % initial guess at duration 

            N = 10;
            prog = DircolTrajectoryOptimization(obj,N,[0 tf0]);
            prog = addStateConstraint(prog,ConstantConstraint(x0),1);
            prog = addStateConstraint(prog,ConstantConstraint(xf),N);
            prog = addStateConstraint(prog,QuadraticConstraint(.5,.5,eye(4),zeros(4,1)),1:N,4:7);
            prog = addInputConstraint(prog,BoundingBoxConstraint([0; 0; 0; 0], [255; 255; 255; 255]),1:N);
            prog = addRunningCost(prog,@cost);
            %prog = addFinalCost(prog,@(t,x) finalCost(t,x,xf));

            %--- snopt options ---%
            %prog = setSolver(prog,'snopt');
            %prog = prog.setSolverOptions('snopt','majoroptimalitytolerance',1e-5);
            %prog = prog.setSolverOptions('snopt','majorfeasibilitytolerance',1e-5);

            t_init = linspace(0,tf0,N);

            %Set initial guess for controls to be trim conditions
            traj_init.u = setOutputFrame(PPTrajectory(foh(t_init,kron(ones(1,N),u0))),getInputFrame(obj));
            %Simulate with u0 input to generate initial guess
            [t_guess, x_guess] = ode45(@(t,x) obj.dynamics(t,x,u0),t_init,x0);
            traj_init.x = setOutputFrame(PPTrajectory(foh(t_guess,x_guess')),getStateFrame(obj));

            tic
            [xtraj,utraj,~,~,info]=solveTraj(prog,t_init,traj_init);
            toc
            
            if nargin == 2 && display
                v = FoamyVisualizer(obj);
                v.playback(xtraj);
            end

            function [g,dg] = cost(dt,x,u)
                R = eye(4);
                g = 0.5*(u-u0)'*R*(u-u0);
                dg = [0, zeros(1,17), (u-u0)'*R];
            end

            function [h,dh] = finalCost(t,x,xf)
                hx = .5*(x(1:3)-xf(1:3))'*(x(1:3)-xf(1:3));
                %hv = .5*(x(8:10)-xf(8:10))'*(x(8:10)-xf(8:10));
                %hw = .5*(x(11:13)-xf(11:13))'*(x(11:13)-xf(11:13));

                %This is to handle the quaternion double cover issue
                hq1 = 1 - xf(4:7)'*x(4:7);
                hq2 = 1 + xf(4:7)'*x(4:7);
                if hq1 <= hq2
                    h = hx + hq1;
                    dh = [0,x(1:3)'-xf(1:3)',-xf(4:7)',zeros(1,6)];
                else
                    h = hx + hq2;
                    dh = [0,x(1:3)'-xf(1:3)',xf(4:7)',zeros(1,6)];
                end
            end
            
        end
        
        function [xtrim, utrim] = findTrim(obj,v,varargin)
            %Trim the plane for level forward flight at speed v

            v0 = [v 0 0]';
            w0 = [0 0 0]';

            u0 = [127; 127; 127; 127; 0; 0; 0];

            function r = trimResidual(x)
                if isempty(varargin)
                    q = [x(5) 1 x(6) x(7)]';
                else
                    q = [1 x(5) x(6) x(7)]';
                end
                q = q/sqrt(q'*q);
                y = [0;0;0;q;0;0;v0;w0;0;0];
                xdot = obj.dynamics(0,y,x(1:4));
                r = xdot(4:16);
            end

            options = optimoptions('fsolve','Algorithm','levenberg-marquardt','FiniteDifferenceType','central','Display','off','TolFun',1e-10,'TolX',1e-10);
            y = fsolve(@trimResidual,u0,options);

            if isempty(varargin)
                qtrim = [y(5) 1 y(6) y(7)]';
            else
                qtrim = [1 y(5) y(6) y(7)]';
            end
            qtrim = qtrim/sqrt(qtrim'*qtrim);
            xtrim = [0;0;0;qtrim;0;0;v0;w0;0;0];
            utrim = y(1:4);

        end

        function x = getInitialState(obj)
            x = zeros(17,1);
            x(4) = .707;
            x(7) = .707;%Quaternion - plane oriented right-side up
        end
        
        function q = getZeroConfiguration(obj)
            q = zeros(9,1);
            q(5) = 1; %Quaternion - plane oriented right-side up
        end
        
        function G = omega_to_qdot(obj,q)
            G = [-q(2) -q(3) -q(4);...
                q(1) -q(4) q(3);...
                q(4) q(1) -q(2);...
                -q(3) q(2) q(1)];
            
        end
        
        function v = constructVisualizer(obj)
          
          v = obj.plane_rbm.constructVisualizer;
%           v = v.setInputFrame(obj.getStateFrame);
        end
  end
    
  
end

