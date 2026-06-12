function x0_curr = Runge_Kutta_Nystrom_solver(x0,u0,dt,f,Index_x,Index_v)
    x0_curr = zeros(numel(Index_x)+numel(Index_v),1);
    x1 = x0_curr;
    x2 = x0_curr;
    x3 = x0_curr;

    x0_x = x0(Index_x);
    x0_v = x0(Index_v);

    dfdx1 = f(x0,u0);
    a1 = dfdx1(Index_v);
    k1 = x0_x + 0.4*dt*x0_v + 0.08*dt*dt*a1;
    
    x1(Index_x) = k1;
    x1(Index_v) = x0_v;
    dfdx2 = f(x1,u0);
    a2 = dfdx2(Index_v);
    k2 = x0_x + 2/3*dt*x0_v + 2/9*dt*dt*a1;

    x2(Index_x) = k2;
    x2(Index_v) = x0_v;
    dfdx3 = f(x2,u0);
    a3 = dfdx3(Index_v);
    k3 = x0_x + 0.8*dt*x0_v + 0.16*dt*dt*(a1+a2);

    x3(Index_x) = k3;
    x3(Index_v) = x0_v;
    dfdx4 = f(x3,u0);
    a4 = dfdx4(Index_v);

    x0_curr(Index_x) = x0_x + dt*(x0_v + (dt/192)*(23*a1 + 75*a2 - 27*a3 + 25*a4));
    x0_curr(Index_v) = x0_v + dt/192*(23*a1 + 125*a2 - 81*a3 + 125*a4);

end