function flag = Timer(tt_, t0, Ts)
    dt = rem(tt_ - t0, Ts);
    % x = fix((tt_ - t0)/Ts);
    % et = tt_ - t0
    % dt = et - Ts*fix(et/Ts);
    % disp(dt)
    if dt<1e-8 || Ts-1e-8<dt
        flag = 1;
    else
        flag = 0;
    end
end