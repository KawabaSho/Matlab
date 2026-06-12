% 2024/05/04 Kawarabayshi
% Ref : https://cattech-lab.com/science-tools/lecture-mini-qr-decomposition/
%     : https://people.inf.ethz.ch/arbenz/ewp/Lnotes/2010/chapter3.pdf
%     : https://en.wikipedia.org/wiki/Hessenberg_matrix
% ex)
A = [1,5,4; 2,4,-7; 2,7,14;];
A = magic(5);
% A = [4,1;1,4];
% tic
[U,T] = Schur_decomposition2(A)
% toc
% tic
% [Q_,R_] = qr(A)
% toc
% Q-Q_
% R-R_
% A = Q*R;


function [U,T] = Schur_decomposition2(A_)
    % m = n, A = U*T*U'
    m_ = size(A_,1);
    I  = eye(m_);
    A = A_;
    U = I;
    for i = 1 : 1000000
    % while true
        [Q,R] = Householder_transformation(A,m_,I);

        A = R*Q;
        U = U*Q;
        % U*U'
        % U*A*U'
        % pause;
        % U = U*Q;
        % max(abs(U0-U))
        % if max(abs(U0-U)) <(1e-15); break;end
        % 1
        % U = U0;
    end
    T=A;
    % Q*R
    % 
    % [Q_,R_] = qr(A_);
    % Q
    % Q_
    % R
    % R_
    % % toc
    % Q+Q_
    % R+R_

    % Q_ = -Q;
    % R_ = -R;
    % U_ = U;
    % T_ = A;
end
function B = Hess(A_,m_,I)
    % m = n, A = Q*R (Q : Orthogonal matrix, R : upper(right) triangular matrix)
    B = A_;
    for i = 2 : m_ - 1
        m = m_ + 1 - i;
        u = B(i:end,i-1);
        norm_ai = sqrt(u'*u);
        u = sign(u(1))*u;
        u(1) = u(1) + norm_ai;
        H = I(1:m,1:m) - (2/(u'*u))*(u*u');
        H_ = I;
        H_(i:end,i:end) = H;
        B = H_*B*H_;
    end
end
function [Q,R] = Householder_transformation(A_,m_,I)
    % m>=n, A = Q*R (Q : Orthogonal matrix, R : upper(right) triangular matrix)
    Q = I;
    R = A_;
    for i = 1 : m_ - 1
        m = m_ + 1 - i;
        u = R(i:end,i);
        norm_ai = sqrt(u'*u);
        u(1) = u(1) - norm_ai;
        % u(1) = sign(u(1))*u(1) + norm_ai;
        H = I(1:m,1:m) - 2/(u'*u)*(u*u');
        H_ = I;
        H_(i:end,i:end) = H;
        R = H_*R;
        Q = Q*H_;
    end
end
