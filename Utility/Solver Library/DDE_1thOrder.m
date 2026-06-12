function dfdx = DDE_1thOrder(func, h)
% Discreate Difference Equation
% func = [f(i-1), f(i)]
% dfdx = [dx dy]./dx

    dfdx = [1, (func(2) - func(1))./h];
end