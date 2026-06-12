function u = Saturation(u,umax,umin)
    Max = numel(u);
    for i = 1 : Max
       if u(i) > umax; u(i) = umax;end
    end
    for i = 1 : Max
       if u(i) < umin; u(i) = umin;end 
    end
end