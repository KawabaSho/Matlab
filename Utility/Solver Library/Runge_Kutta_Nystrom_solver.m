function x0_curr = Runge_Kutta_Nystrom_solver(x0,u0,dt,f,Index_x,Index_v)
    x0_curr = zeros(size(x0));
    x1 = x0_curr;
    x2 = x0_curr;
    x3 = x0_curr;

    x0_x = x0(Index_x);
    x0_v = x0(Index_v);
    dx_v = dt*x0_v;
    dt2  = dt*dt;

    dfdx1 = f(x0,u0);
    a1    = dfdx1(Index_v);
    dx_a1 = dt2*a1;
    k1    = x0_x + 0.4*dx_v + 0.08*dx_a1;
    
    x1(Index_x) = k1;
    x1(Index_v) = x0_v;
    dfdx2 = f(x1,u0);
    a2 = dfdx2(Index_v);
    dx_a2 = dt2*a2;
    k2 = x0_x + 2/3*dx_v + 2/9*dx_a1;

    x2(Index_x) = k2;
    x2(Index_v) = x0_v;
    dfdx3 = f(x2,u0);
    a3 = dfdx3(Index_v);
    k3 = x0_x + 0.8*dx_v + 0.16*(dx_a1+dx_a2);

    x3(Index_x) = k3;
    x3(Index_v) = x0_v;
    dfdx4 = f(x3,u0);
    a4 = dfdx4(Index_v);

    x0_curr(Index_x) = x0_x + dx_v + dt2/192*(23*a1 + 75*a2 - 27*a3 + 25*a4);
    x0_curr(Index_v) = x0_v + dt/192*(23*a1 + 125*a2 - 81*a3 + 125*a4);

end