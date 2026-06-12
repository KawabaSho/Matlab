% 2024/05/04 Kawarabayshi
% Ref : https://cattech-lab.com/science-tools/lecture-mini-qr-decomposition/
%     : https://people.inf.ethz.ch/arbenz/ewp/Lnotes/2010/chapter3.pdf
% ex)
% A = [1,5,4; 2,4,-7; 2,7,14; 1,1,1];
% [Q,R] = Householder_transformation(A)
% A = Q*R;

function [Q,R] = Householder_transformation(A_)
    % m>=n, A = Q*R (Q : Orthogonal matrix, R : upper(right) triangular matrix)
    mn = size(A_);
    if mn(1) < mn(2); Q = [];R = []; return; else; m_ = mn(1); end
    Q = eye(m_);
    R = A_;
    for i = 1 : m_ - 1
        m = m_ + 1 - i;
        ai = R(i:end,i);
        norm_ai = sqrt(ai'*ai);
        u = ai - [norm_ai; zeros(m-1,1)];
        H = eye(m) - 2*(u*u')./(u'*u);
        mh = i-1;
        H_ = [eye(mh), zeros(mh,m); zeros(m,mh), H];
        R = H_*R;
        Q = Q*H_;
    end
end