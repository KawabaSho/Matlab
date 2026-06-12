function dfdx = DDE_2thOrderCentral(func, h)
% Discreate Difference Equation
% func = [f(i-1), f(i+1)]
    dfdx = (func(:, 2) - func(:, 1))./h*0.5;
end