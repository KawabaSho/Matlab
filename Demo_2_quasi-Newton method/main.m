clear,close all
%% compute u given func(u) = 0
func = @(u) (1-u)^2;

u0 = 0;

Opt = quasi_Newton_method(func);

% fprintf("%20.19f\n",fmincon(func,u0,[],[]))
fprintf("%20.19f\n",Opt.solve(u0))