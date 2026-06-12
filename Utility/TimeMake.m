function tt = TimeMake(t0, dt, time_num)
    % tt = [t0, t0 + dt, ..., t0 + dt*(time_num - 1)]
    tt = zeros(1,time_num);
    for i = 1 : time_num - 1
        tt(i + 1) = dt*i;
    end
    tt = tt + t0;
end


