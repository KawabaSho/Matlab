function s_ = B_spline(PP, dt, n)
    m = size(PP,2) + n + 1;
    tt_knot = open_uniform_knot_vector(m, n);
    
    t_    = tt_knot(n + 1) : dt : tt_knot(m - n);
    t_num = size(t_,2);
    s_    = zeros(size(PP,1),t_num);
    
    for iter = 1 : t_num-1
        t = t_(iter);
        % B-spline basis function
        b = zeros(n + 1, m - 1);
        for j = 1 : m - 1
            b(1,j) = (tt_knot(j) <= t) && (t < tt_knot(j+1));
        end
        for k = 1 : n
            for j = 1 : m - k - 1
                t_b1 = tt_knot(j + k) - tt_knot(j);
                t_b2 = tt_knot(j+k+1) - tt_knot(j+1);
                if t_b1; t_b1 = 1/t_b1; end
                if t_b2; t_b2 = 1/t_b2; end
                
                b(k+1,j) = (t-tt_knot(j))*t_b1*b(k,j) + ...
                       (tt_knot(j+k+1) - t)*t_b2*b(k,j+1);
            end
        end
        bn = b(end,:);
        ss = 0;
        for i = 1 : m - n - 1
            ss = ss + PP(:,i)*bn(i);
        end
        s_(:,iter) = ss;
    end
    s_(:,end) = PP(:,end);
end

function tt = open_uniform_knot_vector(m, n)
    tt = zeros(1,m);
    tt(end - n : end) = 1;
    midle_m = m - 2*(n+1);
    for i = 1 : midle_m
        tt(n + 1 + i) = 1/(midle_m + 1)*i;
    end
end