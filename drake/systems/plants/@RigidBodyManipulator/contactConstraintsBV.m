function [phi,normal,B,xA,xB,idxA,idxB,mu,JB,dJB] = contactConstraintsBV(obj,varargin)
% function [phi,normal,B,xA,xB,idxA,idxB,mu,JB,dJB] = contactConstraintsBV(obj,kinsol,allow_multiple_contacts,active_collision_options)
% a variation on the contactConstraints function that uses a different
% representation of the friction cone--basis vectors are the extreme rays
% of the cone (no basis for the normal force)
%
% @retval phi (m x 1) Vector of gap function values (typically contact distance), for m possible contacts
% @retval normal (3 x m) Contact normal vector in world coordinates, points from B to A
% @retval B {2k}(3 x m) friction polyhedron basis vectors
% @retval xA (3 x m) The closest point on body A to contact with body B, relative to body A origin and in body A frame
% @retval xB (3 x m) The closest point on body B to contact with body A, relative to body B origin and in body B frame
% @retval idxA (m x 1) The index of body A. 0 is the special case for the environment/terrain
% @retval idxB (m x 1) The index of body B. 0 is the special case for the environment/terrain
% @retval mu (m x 1) Coefficients of friction
% @retval JB {2k}{m x n} parameterization of the polyhedral approximation of the
%    friction cone in joint coordinates
% @retval dJB {2k}{nm x n} dJBdq
% @param See contactConstraints for argument list description

if nargout > 9,
  [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = contactConstraints(obj,varargin{:});
elseif nargout > 8,
  [phi,normal,d,xA,xB,idxA,idxB,mu,n,D] = contactConstraints(obj,varargin{:});
else
  [phi,normal,d,xA,xB,idxA,idxB,mu] = contactConstraints(obj,varargin{:});
end

if nargout > 1,
  nk = length(d);  
  B = cell(1,2*nk);
  JB = cell(1,2*nk);
  dJB = cell(1,2*nk);
  
  muI = sparse(diag(mu));
  norm_mat = sparse(diag(1./sqrt(1 + mu.^2)));
  if nargout > 9
    nv = obj.getNumVelocities;
    muI2 = sparse(diag(repmat(mu,nv,1)));
    norm_mat2 = sparse(diag(repmat(diag(norm_mat),nv,1)));
  end
  for k=1:nk,
    B{k} = (normal + d{k}*muI)*norm_mat;
    B{nk+k} = (normal - d{k}*muI)*norm_mat;
    if nargout > 8
      JB{k} = norm_mat*(n + muI*D{k});
      JB{nk+k} = norm_mat*(n - muI*D{k});
    end
    if nargout > 9
      dJB{k} = norm_mat2*(dn + muI2*dD{k});
      dJB{nk+k} = norm_mat2*(dn - muI2*dD{k});
    end
  end
end
end
