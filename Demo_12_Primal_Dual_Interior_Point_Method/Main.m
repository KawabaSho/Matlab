%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Demo_12_Primal_Dual_Interior_Point_Method
%
%   Finding optimal solution of constrained optimization problems using the
%   primal-dual interior point method.
%
%
%   Copyright © 2026 Kawarabayashi
%   Released under the MIT license
%   https://opensource.org/licenses/mit-license.php
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%{
    Reference: ISBN978-4-274-23385-2
%}
clear,close all

flag_solver = 1; % 1:solver_PDIP, 2:fmincon(interior-point)

%% Problem
f    = @(x) 0.1*(x-1).*(x-3)+sin(x); % Objective function
f_x  = @(x) 0.2*(x-2)+cos(x);        % Gradient of objective function
x0 = 0.5;                              % Initial guess
% subject to  
x_f = [-3,8];                        % a1 <= x <= a2. % sample: x_f = [-10,-5; 2, 6;]; 
function c = g(x,a1,a2)              % format of constraint
    c = [  x - a2;  
           a1 -  x;];  
end  
function c = g_x(x,a1,a2)            % format of constraint gradient
    sizemax = numel(a2);
    sizemin = numel(a1);
    c = [eye(sizemax);-eye(sizemin)];
end

%% solve
tic
switch flag_solver
    case 1
        op = solver_PDIP(f,f_x,x0,"con",@(x)g(x,x_f(:,1),x_f(:,2)),"con_x",@(x)g_x(x,x_f(:,1),x_f(:,2))); % Initializing
        % op = solver_PDIP(@f,@f_jacobian,x0);                                                            % non-constrainted optimization
        op.Run();
        x_opt = op.getSolution;
    case 2
        x_opt = fmincon(f,x0,[],[],[],[],[],[],@(x)G(x,x_f(:,1),x_f(:,2)));
    otherwise
end
toc

fprintf("Optimal solution (x, f(x)) = (%8f, %8f)\n", x_opt,f(x_opt))

%% Rendering
dx = 0.01;
x  = -6:dx:10;
plot(x,f(x),"Color","b","DisplayName","Objective function")
hold on
plot(x0,f(x0),".","Color","g","DisplayName","$x_{0}$(Initial guess)","MarkerSize",18)
plot(x_opt,f(x_opt),".","Color","r","DisplayName","$x_{\rm{opt}}$","MarkerSize",18)
for i = 1 : size(x_f,1)
    if i == 1
        xregion(x_f(i,1),x_f(i,2),"FaceColor",[0.73 0.83 0.95],"DisplayName","Feasible area")
    else
        set(get(get(xregion(x_f(i,1),x_f(i,2),"FaceColor",[0.73 0.83 0.95]),"Annotation"),'LegendInformation'),'IconDisplayStyle','off')
    end
end

hold off
% grid on
legend("FontSize",12,"Interpreter","latex")
xlabel("$x$","Interpreter","latex","FontSize",18)
ylabel("$f(x)$","Interpreter","latex","FontSize",18)


% constraints for fmincon
function [con, ceq] = G(x,a1,a2)
    con = g(x,a1,a2);
    ceq = [];
end





