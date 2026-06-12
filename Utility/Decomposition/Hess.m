% Compute Hessenberg matrix
% 
% 2024/05/04 Kawarabayshi
% Ref : https://cattech-lab.com/science-tools/lecture-mini-qr-decomposition/
%     : https://people.inf.ethz.ch/arbenz/ewp/Lnotes/2010/chapter3.pdf
%     : https://en.wikipedia.org/wiki/Hessenberg_matrix
% ex)
% A = [1,5,4; 2,4,-7; 2,7,14;];
% B = Hess2(A);
% V = hess(A); % matlab function
% B - V % error

function B = Hess(A_)
    % m=n, A = Q*R (Q : Orthogonal matrix, R : upper(right) triangular matrix)
    m_ = size(A_,1);
    I  = eye(m_);
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
function bool = sign(x)
    if  x >= 0
        bool = 1;
    else 
        bool = -1;
    end
end
