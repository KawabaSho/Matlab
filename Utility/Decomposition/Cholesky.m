function [flag, S] = Cholesky(A)
    % A = S*S' , flag 1 : A is not positive definite...
    flag = 0;
    mn = size(A);
    for i = 1 : mn(1)
        flag = flag + (A(i,i)<=0);
    end
    if flag
        flag = 1;
        S = [];
        return
    end
    S  = zeros(mn);
    S(1,1) = sqrt(A(1,1));
    for i = 2:mn(1)
        % (i,j)
        for j = 1 : i - 1
            sig1 = 0;
            for k = 1 : j - 1
                sig1 = sig1 + S(i,k)*S(j,k);
            end
            S(i,j) = (A(i,j) - sig1)/S(j,j);
        end
        % (i,i)
        sig2 = 0;
        for k = 1 : i - 1
            sig2 = sig2 + S(i,k)*S(i,k);
        end
        F = A(i,i) - sig2;
        if F<0
            flag = 1;
            S = [];
            return
        end
        S(i,i) = sqrt(F);
    end
end