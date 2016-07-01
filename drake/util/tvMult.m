function X = tvMult(T,v,ind)
% [X] = tvMult(T,x,ind)
% Returns the product of a rank-3 tensor T and a vector v with contraction
% over any index 1, 2, or 3. This corresponds to T_{ijk} v_i, T_{ijk} v_j, 
% or T_{ijk} v_k in index notation. Dimensions must be compatible.
% @param T a tensor (an n x m x p array)
% @param v a vector
% @param ind an index to sum over (1, 2, or 3)
% @retval X the resulting matrix

if ind == 1
    X = zeros(size(T,2),size(T,3));
    for j = 1:size(T,2)
        for k = 1:size(T,3)
            X(j,k) = T(:,j,k)'*v;
        end
    end
elseif ind == 2
    X = zeros(size(T,1),size(T,3));
    for j = 1:size(T,1)
        for k = 1:size(T,3)
            X(j,k) = T(j,:,k)*v;
        end
    end
elseif ind == 3
    X = zeros(size(T,1),size(T,2));
    for j = 1:size(T,1)
        for k = 1:size(T,2)
            X(j,k) = squeeze(T(j,k,:))'*v;
        end
    end
else
    error('Invalid index');
end
end