function dfdx = DDE_4thOrderCentral(func, h)
% Discreate Difference Equation
% Richardson extrapolation
% func = [f(i-2), f(i-1), f(i+1), f(i+2)]
    dfdx = [1, (-func(4) + 8*func(3) - 8*func(2) + func(1))./h/12];
end