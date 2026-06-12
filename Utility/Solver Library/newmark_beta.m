function x0_curr = newmark_beta(x0,u0,dt, f, dfdx, gamma, beta, IterationMax, Tolarance, Index_x, Index_v, Index_a)
    % x0 = [x1;x2;v1;v2;a1;a2], f(x) = 0
    x0_curr = zeros(size(x0));

    x_pre  = x0(Index_x);
    v_pre  = x0(Index_v);
    a_pre  = x0(Index_a);

    x0_curr(Index_x) = x_pre + v_pre*dt + a_pre*(0.5 - beta)*dt^2+a_pre*beta*dt^2;
    x0_curr(Index_v) = v_pre + a_pre*(1 - gamma)*dt + a_pre*gamma*dt;
    x0_curr(Index_a) = a_pre;%zeros(size(a_pre));
    eps = f(x0_curr,u0);

    % newton raphson method
    for i = 1 : IterationMax
        dx  = - dfdx(x0_curr)\eps;
        err = norm(dx,'inf');
        if err < Tolarance; break; end
        dx = min(1,Tolarance/err)*dx;
        x0_curr(Index_x) = x0_curr(Index_x) + dx;
        x0_curr(Index_v) = x0_curr(Index_v) + dx/gamma/dt;
        x0_curr(Index_a) = x0_curr(Index_a) + dx/beta/dt^2;
        eps = f(x0_curr,u0);
    end
    if i == IterationMax; fprintf("Newton Raphson Method : Iteration has reached maximum.\n"); end
end
