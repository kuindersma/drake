classdef LinearComplementarityProgram < NonlinearProgram
%
% Solve the standard linear complementarity problem:
%     z >= 0, Mz + q >= 0, z'*(Mz + q) = 0
%
% Required input:
%     M(n,n)  - matrix
%     q(n)    - vector
% 
% Output:
%     z(n)    - solution
%
% Optional input:
%     l(n)    - lower bounds                       default: zero
%     u(n)    - upper bounds                       default: infinity
%     z(n)    - starting point                     default: zero
%
% The optional lower and upper bounds are used to define a linear mixed 
% complementarity problem (box constrained variational inequality).
%       l <= z <= u
%       where l_i < z_i < u_i  => (Mz + q)_i = 0
%             l_i = z          => (Mz + q)_i >= 0
%             u_i = z          => (Mz + q)_i <= 0
% 
  


end