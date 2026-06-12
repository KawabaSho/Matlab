function dfdx = DDE_2thOrder(func, h)
% Discreate Difference Equation
% func = [f(i-2), f(i-1), f(i)]

    dfdx = (3*func(:, 3) - 4*func(:, 2) + func(:, 1))./h*0.5;
end