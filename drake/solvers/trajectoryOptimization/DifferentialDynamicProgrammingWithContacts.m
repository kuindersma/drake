classdef DifferentialDynamicProgrammingWithContacts < DifferentialDynamicProgramming
  %DIFFERENTIALDYNAMICPROGRAMMINGWITHCONTACTS A class for performing DDP or iLQR
  %   including contact forces and constraints
  %
  % This class assumes that there are a fixed number (N) time steps, and
  % that the trajectory is discreteized into timesteps h (N-1), state x
  % (N), control input u (N), and contact force coefficients f (N)
   
  properties
    gurobi_options
  end
  methods
    function obj = DifferentialDynamicProgrammingWithContacts(plant,compute_second_derivatives,options)
      % @param plant the robot
      % @param N the number of time samples
      % @param options 

      if nargin < 3
        options = struct();
      end

      obj = obj@DifferentialDynamicProgramming(plant,compute_second_derivatives,options);

      
      obj.gurobi_options.outputflag = 0; % not verbose
      obj.gurobi_options.method = 1; % -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier
      obj.gurobi_options.presolve = 0;
      % obj.gurobi_options.prepasses = 1;

      if obj.gurobi_options.method == 2
        obj.gurobi_options.bariterlimit = 100; % iteration limit
        obj.gurobi_options.barhomogeneous = 0; % 0 off, 1 on
%         obj.gurobi_options.barconvtol = 5e-5;
      end

     end

    function [xtraj,utraj,ztraj] = solveTraj(obj,xtraj0,utraj0,ztraj0,N,h,iter)
      xtraj=xtraj0;
      utraj=utraj0;
      ztraj=ztraj0;

      p=obj.plant;
      nx=size(xtraj,1);
      nu=size(utraj,1);
      nz=size(ztraj,1);

      phi_traj = zeros(nz,N);

      if p.twoD % using hard-coded d for now
        d=2;
      else
        d=4;
      end
      
      prev_cost=inf;
    	for j=1:iter
        [V,Vx,Vxx] = feval(obj.final_cost_function,xtraj(:,N));
        Ku = cell(1,N-1);
        ku = cell(1,N-1);
        Kz = cell(1,N-1);
        kz = cell(1,N-1);
        Qx = cell(1,N-1);
        Qu = cell(1,N-1);
        Qz = cell(1,N-1);
        Qxx = cell(1,N-1);
        Qux = cell(1,N-1);
        Qzx = cell(1,N-1);
        Quu = cell(1,N-1);
        Qzu = cell(1,N-1);
        Qzz = cell(1,N-1);
          
        % below two terms for computing expect cost reduction in line search
        sum_kT_Quu_k = 0;
        sum_kT_Qu = 0;
        for i=N-1:-1:1
          x = xtraj(:,i); 
          u = utraj(:,i);
          z = ztraj(:,i);

          x_ind = 1:nx;
          u_ind = nx+(1:nu);
          z_ind = nx+nu+(1:nz);
          
%           [~,df] = geval(@p.updateWithBasisForces,i*p.timestep,x,u,z,struct('grad_method','numerical'));

          [~,df,phi] = p.updateWithBasisForces(i*p.timestep,x,u,z);
          Fx = full(df(:,1+x_ind)); % +1 to skip time argument
          Fu = full(df(:,1+u_ind));
          Fz = full(df(:,1+z_ind));
                    
          phi_traj(:,i) = repmat(max(eps,phi),d,1); 
          
%           [~,dg,d2g] = geval(obj.running_cost_function,x,u,z,struct('grad_method','taylorvar'));
          [~,dg,d2g] = feval(obj.running_cost_function,x,u,z);
          
          d2g =reshape(d2g,nx+nu+nz,[]);
          gx = dg(x_ind);
          gu = dg(u_ind);
          gz = dg(z_ind);
          Gxx = d2g(x_ind,x_ind);
          Gux = d2g(u_ind,x_ind);
          Guu = d2g(u_ind,u_ind);
          Gzx = d2g(z_ind,x_ind);
          Gzu = d2g(z_ind,u_ind);
          Gzz = d2g(z_ind,z_ind);
          
          mu = 0.1;
          Qx{i} = gx + Vx*Fx;
          Qu{i} = gu + Vx*Fu;
          Qz{i} = gz + Vx*Fz;
          Qxx{i} = Gxx + Fx'*(Vxx + mu*eye(nx))*Fx;
          Qux{i} = Gux + Fu'*(Vxx + mu*eye(nx))*Fx;
          Qzx{i} = Gzx + Fz'*(Vxx + mu*eye(nx))*Fx;
          Quu{i} = Guu + Fu'*(Vxx + mu*eye(nx))*Fu;
          Qzu{i} = Gzu + Fz'*(Vxx + mu*eye(nx))*Fu;
          Qzz{i} = Gzz + Fz'*(Vxx + mu*eye(nx))*Fz;
          
          zub = -1./(1-exp(phi_traj(:,i)));
          if any(isinf(zub))
            keyboard
          end
%           Ain = [-eye(nu+nz); eye(nu+nz)]; 
%           bin = [-(p.umin-utraj(:,i)); -ztraj(:,i); (p.umax-utraj(:,i)); zub-ztraj(:,i)];
          
          Q=[Quu{i},Qzu{i}'; Qzu{i},Qzz{i}];
%           [result,info_fqp] = fastQPmex({Q},[Qu{i},Qz{i}]',Ain,bin,[],[],[]);
          if 1%info_fqp<0
            % then call gurobi
%             disp('QPController: failed over to gurobi');
            model.Q = sparse(Q);
            model.A = sparse(zeros(1,nu+nz));
            model.rhs = 0;
            model.sense = '=';
            model.lb = [p.umin-utraj(:,i); -ztraj(:,i)];
            model.ub = [p.umax-utraj(:,i); zub-ztraj(:,i)];
            model.obj = [Qu{i},Qz{i}]';
            if obj.gurobi_options.method==2
              % see drake/algorithms/QuadraticProgram.m solveWGUROBI
              model.Q = .5*model.Q;
            end

            result = gurobi(model,obj.gurobi_options);

            result = result.x;
          end
          

%           a = result(1:nu);
%           b = result(nu+(1:nz));

          ku{i} = result(1:nu);
          kz{i} = result(nu+(1:nz));
%           active_ind = sum(reshape(abs(Ain*result - bin)<1e-6,nu+nz,2),2) > 0;
%           Quu_free = Quu{i};
%           Quu_free(active_ind(1:nu),:) = 0;
%           Qzz_free = Qzz{i};
%           Qzz_free(active_ind(nu+(1:nz)),:) = 0;
%           Ku{i} = -pinv(Quu_free)*Qux{i};
%           Kz{i} = -pinv(Qzz_free)*Qzx{i};
        
%           sum_kT_Quu_k = sum_kT_Quu_k + ku{i}'*Quu*ku{i} + kz{i}'*Qzz*kz{i};
%           sum_kT_Qu = sum_kT_Qu + ku{i}'*Qu + kz{i}'*Qz;

          Quu_inv = inv(Quu{i});
          Qzz_inv = inv(Qzz{i});

          A = (eye(nu) - Quu_inv*Qzu{i}'*Qzz_inv*Qzu{i})\Quu_inv;
%           ku{i} = -A*(Qu{i}' - Qzu{i}'*Qzz_inv*Qz{i}');
          Ku{i} = -A*(Qux{i} - Qzu{i}'*Qzz_inv*Qzx{i});
           
%           kz{i} = -Qzz_inv*(Qz{i}' + Qzu{i}*ku{i});
          Kz{i} = -Qzz_inv*(Qzx{i} + Qzu{i}*Ku{i});

%           i
%           try
%             valuecheck(a,ku{i},1e-5);
%             valuecheck(b,kz{i},1e-5);
%           catch
%             keyboard
%           end
          
          V = V - 0.5*(ku{i}'*Quu{i}*ku{i} + kz{i}'*Qzz{i}*kz{i});
          Vx = Qx{i} - (Ku{i}'*Quu{i}*ku{i} + Kz{i}'*Qzz{i}*kz{i})';
          Vxx = Qxx{i} - Ku{i}'*Quu{i}*Ku{i} - Kz{i}'*Qzz{i}*Kz{i}; 

        end

        xtraj_new = xtraj;
        utraj_new = utraj;
        ztraj_new = ztraj;

        alpha = 0.005;
        line_search_incomplete = true;
        while (line_search_incomplete)
          cost = 0;
          for i=1:N-1
            dx = xtraj_new(:,i)-xtraj(:,i);
            du = alpha*ku{i} + Ku{i}*dx;
            dz = alpha*kz{i} + Kz{i}*dx;
            utraj_new(:,i) = utraj(:,i) + du;
            utraj_new(:,i) = min(p.umax,max(p.umin,utraj_new(:,i)));
            
            ztraj_new(:,i) = ztraj(:,i) + dz;
            ztraj_new(:,i) = max(0,ztraj_new(:,i));
            
%             norm_u = norm(utraj_new(:,i))
%             norm_z = norm(ztraj_new(:,i))

%             norm(dx)

%             du2 = ku{i} + Ku{i}*dx;
%             dz2 = kz{i} + Kz{i}*dx;
%             if ~valuecheck(Qu{i}' + Qux{i}*dx + Qzu{i}'*dz2 + Quu{i}*du2,0,1e-5)
%               keyboard
%             end
%             if ~valuecheck(Qz{i}' + Qzx{i}*dx + Qzu{i}*du2 + Qzz{i}*dz2,0,1e-5)
%               keyboard
%             end


            try
%               xtraj_new(:,i+1) = p.updateWithBasisForces(i*p.timestep,xtraj_new(:,i),utraj_new(:,i),ztraj_new(:,i));
              xtraj_new(:,i+1) = p.update(i*p.timestep,xtraj_new(:,i),utraj_new(:,i));
            catch
              blah=1;
            end
            cost = cost + obj.running_cost_function(xtraj_new(:,i),utraj_new(:,i),ztraj_new(:,i));
            
          end
          cost = cost + obj.final_cost_function(xtraj_new(:,N));

          if obj.options.enable_line_search
            expected_cost_reduction = 0.5*alpha^2 * sum_kT_Quu_k + alpha*sum_kT_Qu;
            z=-(prev_cost-cost)/expected_cost_reduction;
            if z > 0.5 || alpha < 1e-5
              line_search_incomplete = false;
            else
              alpha = alpha*.5;
            end
          else
            line_search_incomplete = false;
          end
          
        end
        cost
        prev_cost=cost;
    
%         norm(xtraj-xtraj_new)
        
        xtraj=xtraj_new;
        utraj=utraj_new;
        ztraj=ztraj_new;
        
        if obj.options.enable_visualizer && ~mod(j,obj.options.visualizer_period)
          xtraj_pp = obj.reconstructStateTrajectory(xtraj,0,h);
          obj.options.visualizer.playback(xtraj_pp);
        end
      end

    end
  end  
end
