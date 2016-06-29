 function y = tvMult(T,x,ind)
    if ind == 1
        y = zeros(size(T,2),size(T,3));
        for j = 1:size(T,2)
            for k = 1:size(T,3)
                y(j,k) = T(:,j,k)'*x;
            end
        end
    elseif ind == 2
        y = zeros(size(T,1),size(T,3));
        for j = 1:size(T,1)
            for k = 1:size(T,3)
                y(j,k) = T(j,:,k)*x;
            end
        end
    else %ind == 3
        y = zeros(size(T,1),size(T,2));
        for j = 1:size(T,1)
            for k = 1:size(T,2)
                y(j,k) = squeeze(T(j,k,:))'*x;
            end
        end
    end
end