function  u0_opt = CheckRun(xx)
    goal = [9;9];
    rr = goal - xx(1:2);
    r = sqrt(rr'*rr);
    nn = rr/r;
    if r < 0.4
        v = r;
    else
        v = 0.4;
    end
    u0_opt = [v; -0.5*(nn(1)*sin(xx(3)) - nn(2)*cos(xx(3)))];
end