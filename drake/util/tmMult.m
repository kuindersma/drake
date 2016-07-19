function X = tmMult(T,M,ind_t,ind_m)
% [X] = tmMult(T,x,ind)
% Returns the product of a rank-3 tensor T and a matrix M with contraction
% over any pair of indeces. Dimensions must be compatible.
% @param T a tensor (an n x m x p array)
% @param M a matrix (an q x r array)
% @param ind_t an index to sum over (1, 2, or 3)
% @param ind_m and index to sum over (1 or 2)
% @retval X the resulting rank-3 tensor

if ind_t == 1
    if ind_m == 1
        X = zeros(size(T,2),size(T,3),size(M,2));
        for j = 1:size(T,2)
            for k = 1:size(T,3)
                X(j,k,:) = T(:,j,k)'*M;
            end
        end
    elseif ind_m == 2
        X = zeros(size(T,2),size(T,3),size(M,1));
        for j = 1:size(T,2)
            for k = 1:size(T,3)
                X(j,k,:) = T(:,j,k)'*M';
            end
        end
    else
        error('Invalid index');
    end
elseif ind_t == 2
    if ind_m == 1
        X = zeros(size(T,1),size(T,3),size(M,2));
        for j = 1:size(T,1)
            for k = 1:size(T,3)
                X(j,k,:) = T(j,:,k)*M;
            end
        end
    elseif ind_m == 2
        X = zeros(size(T,1),size(T,3),size(M,1));
        for j = 1:size(T,1)
            for k = 1:size(T,3)
                X(j,k,:) = T(j,:,k)*M';
            end
        end
    else
        error('Invalid index');
    end
    
elseif ind_t == 3
    if ind_m == 1
        X = zeros(size(T,1),size(T,2),size(M,2));
        for j = 1:size(T,1)
            for k = 1:size(T,2)
                X(j,k,:) = squeeze(T(j,k,:))'*M;
            end
        end
    elseif ind_m ==2
        X = zeros(size(T,1),size(T,2),size(M,1));
        for j = 1:size(T,1)
            for k = 1:size(T,2)
                X(j,k,:) = squeeze(T(j,k,:))'*M';
            end
        end
    else
        error('Invalid index');
    end
else
    error('Invalid index');
end
end