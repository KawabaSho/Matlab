function dfdx = CDE_4thOrderCentral(func, x0, h)
% Continuous Difference Equation
% Richardson extrapolation
    dfdx = (-func(:, x0 + 2*h) + 8*func(:, x + h) - 8*func(:, x0 - h) + func(:, x - 2*h))./h/12;
end