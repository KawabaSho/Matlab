%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Primal_Dual_Interior_Point_Method
%
%   Copyright © 2026 Kawarabayashi
%   Released under the MIT license
%   https://opensource.org/licenses/mit-license.php
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
classdef solver_PDIP < handle
    % Primal-dual interior point method class
    % Input: function handle f(x), df/dx and initial guess x0
    % option: equality and inequality constraints 
    % (If you consider these constraints, please also input the gradient)
    %{
        example 1:
        op = solver_PDIP(@f,@f_x,x0,"con",@g,"con_x",@g_x);
        op.Run();
        x_opt = op.getSolution;
        
        example 2:
        op = solver_PDIP(@f,@f_x,x0);
        op.Run();
        x_opt = op.getSolution;
        min_f = op.getCost;

        Reference: ISBN978-4-274-23385-2
    %}
    properties
        % Problem
        f        % cost
        f_x      % Jacobian
        f_xx     % Hessian
        f_name
        ceq      % equal
        ceq_x    % jacobian
        ceq_xx
        ceq_dim
        con      % inequal
        con_x
        con_xx
        con_dim
        x0       % initial guess
        x_dim

        % solution
        Ff   = []
        Xf   = []
        Yf   = []
        Zf   = []
    end
    properties (Hidden)
        F_        % cost
        F_X
        H_        % equal [h, g]
        H_X
        H_dim
        X0       % initial guess [state; state slack; g slack]
        p_s      % state pointer of X
        p_ss     % slack for X
        p_gs     % slack for Inequality
        X_OPT    % optimal solution
        X_dim
        
        pX_dim
        pY_dim
        pZ_dim

        L_
        L_X
        L_XX

        J_k   % coefficient matrix      J_k(X,Y,Z)
        r_0   % KKT conditions          r_0(X,Y,Z)
        r     % centerd KKT conditions  r_0(X,Y,Z,nu)

        P_     % log_barrier_function
        P_L    % L1_barrier_penalty_function
        dP_    % dP = P_L - P

        % flag
        UnconstrainedFlag
        constraintsID
        BFGSisActivated
    end
    properties (Constant,Hidden)
        e_init     = 1.0  % Initial guess X0 in PDIP initializing
        nu_init    = 0.1 % Barrier parameter

        M_c        = 1    % Mc > 0
        tau        = 0.25

        xi         = 0.01
        gamma      = 0.99

        beta       = 0.5

        M_L        = 1
        M_U        = 2

        nu_tolerance = 1e-8

        outerloopMaxiteration = 10
        innerloopMaxiteration = 10

        armijo_maxiteration   = 20
    end
    properties
        SystemReport        = 0
        SystemReport_string = ""
    end

    methods
        function obj = solver_PDIP(f, f_x, x0,option)
            %solver_PDIP
            %   f   is a cost function.
            %   f_x is a gradient of the cost function (1,n)vector.
            %   h   is a eqaul.    (n,1)vector ceq
            %   g   is an inequal. (n,1)vector con
            arguments
                f 
                f_x
                x0
                option.ceq   = []
                option.ceq_x = []
                option.con   = []
                option.con_x = []
                option.f_xx  = []
            end
            if isempty(option.f_xx)
                obj.BFGSisActivated = 1;
            else
                % unsupported
                obj.BFGSisActivated = 0;
                obj.f_xx = option.f_xx;
                obj.SystemReport        = 0x01;
                obj.SystemReport_string = obj.SystemReport_string + "f_xx is not supported. Please use modified BFGS\n";
            end
            if isempty(option.ceq)
                obj.ceq_dim = 0;
            else
                obj.ceq = option.ceq;
                obj.ceq_x = option.ceq_x;
                obj.ceq_dim = length(option.ceq(x0));
            end            
            if isempty(option.con)
                obj.con_dim = 0;
            else
                obj.con = option.con;
                obj.con_x = option.con_x;
                obj.con_dim = length(option.con(x0));
            end
            xdim = length(x0);
            obj.x_dim = xdim;
            obj.f = f;
            obj.f_x = f_x; % gradient
            obj.f_name = func2str(f);
            obj.p_s  = 1           : xdim;      % decition variables
            obj.p_ss = xdim+1      : xdim+xdim; % slack variables for decition variables

            % Re-definition in PDIP
            obj.H_dim = obj.ceq_dim + obj.con_dim;
            obj.X_dim = xdim + xdim + obj.con_dim;
            obj.F_ = @(X) f(X(obj.p_s)-X(obj.p_ss));
            obj.X0 = [x0.*(x0>0) + obj.e_init.*(x0<=0);
                      obj.e_init.*(x0<=0) - x0.*(x0<=0);
                      zeros(obj.con_dim,1)]; % The operation of setting the initial guess to a positive value within optimization
            
            % Conditional Branching Based on Constraint Type
            obj.UnconstrainedFlag = 0;
            dxdX = [eye(obj.x_dim), -eye(obj.x_dim)]; % Decision variables are [main problem variables, their slack variables, slack variables of inequality constraints].
            obj.F_X = @(X) [f_x(X(obj.p_s)-X(obj.p_ss))*dxdX, zeros(1,obj.con_dim)];
            if obj.ceq_dim
                if obj.con_dim % ceq con
                    obj.p_gs = xdim+xdim+1 : xdim+xdim+obj.con_dim;
                    obj.H_   = @(X) [obj.ceq(X(obj.p_s)-X(obj.p_ss)); obj.con(X(obj.p_s)-X(obj.p_ss))+X(obj.p_gs)];
                    obj.H_X  = @(X) [obj.ceq_x(X(obj.p_s)-X(obj.p_ss))*dxdX, zeros(obj.ceq_dim,obj.con_dim);
                                     obj.con_x(X(obj.p_s)-X(obj.p_ss))*dxdX, eye(obj.con_dim)];
                    obj.X0 = [x0.*(x0>0) + obj.e_init.*(x0<=0);
                              obj.e_init.*(x0<=0) - x0.*(x0<=0);
                              -obj.con(x0)];
                else % ceq
                    obj.H_   = @(X) obj.ceq(X(obj.p_s)-X(obj.p_ss));
                    obj.H_X  = @(X) obj.ceq_x(X(obj.p_s)-X(obj.p_ss))*dxdX;
                end
            else
                if obj.con_dim % con
                    obj.p_gs  = xdim+xdim+1 : xdim+xdim+obj.con_dim;
                    obj.H_    = @(X) obj.con(X(obj.p_s)-X(obj.p_ss))+X(obj.p_gs);
                    obj.H_X   = @(X) [obj.con_x(X(obj.p_s)-X(obj.p_ss))*dxdX, eye(obj.con_dim)];
                    obj.X0 = [x0.*(x0>0) + obj.e_init.*(x0<=0);
                              obj.e_init.*(x0<=0) - x0.*(x0<=0);
                              -obj.con(x0)];
                else % Only X >= 0
                    obj.UnconstrainedFlag = 1;
                    obj.H_    = @(X)[];
                    obj.H_X   = @(X)[];
                end
            end
            if obj.UnconstrainedFlag
                obj.L_  = @(X,Y,Z,F,H)            F - Z'*X;
                obj.L_X = @(X,Y,Z,F_X,H_X)        F_X - Z';
                obj.J_k = @(X,Z,Lxx,Hx)           [Lxx,-eye(obj.X_dim);diag(Z),diag(X)];
                obj.r_0 = @(X,Z,Lx)               [Lx'; X.*Z];
                obj.r   = @(X,Z,Lx,H,nu)          [Lx'; X.*Z - nu];
                obj.P_   = @(X,F,H,nu,rho)        F - nu*sum(log(X),"all");
                obj.P_L = @(X,dX,F,Fx,H,Hx,nu,rho)F - nu*sum(log(X),"all") + Fx*dX - nu*sum(dX./X,"all");
                obj.dP_  = @(X,dX,Fx,H,Hx,nu,rho) Fx*dX - nu*sum(dX./X,"all");
            else
                obj.L_  = @(X,Y,Z,F,H)            F + Y'*H - Z'*X;
                obj.L_X = @(X,Y,Z,F_X,H_X)        F_X + Y'*H_X - Z';
                obj.J_k = @(X,Z,Lxx,Hx)           [[Lxx,Hx',-eye(obj.X_dim)];...
                                                  [Hx,zeros(obj.H_dim,obj.H_dim + obj.X_dim)];...
                                                  [diag(Z),zeros(obj.X_dim,obj.H_dim),diag(X)]];
                obj.r_0 = @(X,Z,Lx,H)             [Lx'; H; X.*Z];
                obj.r   = @(X,Z,Lx,H,nu)          [Lx'; H; X.*Z - nu];
                obj.P_   = @(X,F,H,nu,rho)        F - nu*sum(log(X),"all") + rho*sum(abs(H),"all");
                obj.P_L = @(X,dX,F,Fx,H,Hx,nu,rho)F + Fx*dX - nu*sum(log(X) + dX./X,"all") + rho*sum(abs(H + Hx*dX),"all");
                obj.dP_  = @(X,dX,Fx,H,Hx,nu,rho) Fx*dX - nu*sum(dX./X,"all") + rho*sum(abs(H + Hx*dX),"all") - rho*sum(abs(H),"all");
            end
            obj.pX_dim = 1 : obj.X_dim;
            obj.pY_dim = (obj.X_dim + 1) : (obj.X_dim + obj.H_dim);
            obj.pZ_dim = (obj.X_dim + obj.H_dim + 1) : (obj.X_dim + obj.H_dim + obj.X_dim);

            if obj.SystemReport
                obj.SystemReport_string = obj.SystemReport_string + "Initialization failed.\n";
                error(obj.SystemReport_string);
            end
        end
        % function out = Lagrangian_log_barrier_function(obj, x, y, v)
        %     out = obj.F_(x) - v*sum(log(x),"all") + sum(y'*obj.H_(x));
        % end
        % function out = log_barrier_function(obj, x, v)
        %     out = obj.F_(x) - v*sum(log(x),"all");
        % end
        % function out = L1_barrier_penalty_function(obj, x, v, ro)
        %     out = obj.F_(x) - v*sum(log(x),"all") + ro*sum(abs(obj.H_(x)),"all");
        % end
        % function r = centered_KKT(x,y,z,nu)
        %     r = [obj.L_x(x,y,z);
        %          obj.H(x);
        %          x.*z - nu;]; % x>0, z>0
        % end
        function xopt = getSolution(obj)
            xopt = obj.Xf(obj.p_s)-obj.Xf(obj.p_ss);
        end
        function xopt = getCost(obj)
            xopt = obj.Ff;
        end
        function xopt = getValue(obj,X)
            xopt = X(obj.p_s)-X(obj.p_ss);
        end
        
        function Run(obj)
            nu     = obj.nu_init;
            Xk     = obj.X0;
            Yk     = zeros(obj.H_dim,1);
            Zk     = zeros(obj.X_dim,1)+.1;
            xdim   = obj.X_dim;
            p_xdim = obj.pX_dim;
            p_ydim = obj.pY_dim;
            p_zdim = obj.pZ_dim;
            
            % Initialize
            BE       = eye(xdim);
            B        = BE;
            F        = obj.F_(Xk);
            H        = obj.H_(Xk);
            Fx       = obj.F_X(Xk);
            Hx       = obj.H_X(Xk);
            Lx       = obj.L_X(Xk,Yk,Zk,Fx,Hx);
            rk       = obj.r(Xk,Zk,Lx,H,nu);
            dW       = -inv(obj.J_k(Xk,Zk,BE,Hx))*rk; % []
            dX       = dW(p_xdim);
            dY       = dW(p_ydim);
            dZ       = dW(p_zdim);
            nx       = -Xk./dX;
            rho      = max(abs(Yk+dY)) + 1;
            % Armijo condition
            alpha_xk = min([1, obj.gamma*min(nx(dX<0))]);
            P        = obj.P_(Xk,F,H,nu,rho);
            xi_dP    = obj.xi*obj.dP_(Xk,dX,Fx,H,Hx,nu,rho);
            
            % report_flag = [0x00,0x00,0x00];   %%%% debug mode %%%%
            for iter_outer = 1 : obj.outerloopMaxiteration
                for iter_inner = 1 : obj.innerloopMaxiteration
                    % primal problem
                    % Armijo condition
                    for iter_almijo = 1 : obj.armijo_maxiteration
                        Xnew  = Xk + alpha_xk*dX;
                        Pcurr = obj.P_(Xnew,obj.F_(Xnew),obj.H_(Xnew),nu,rho);
                        if Pcurr <= P + alpha_xk*xi_dP
                            break;
                        end
                        alpha_xk = alpha_xk * obj.beta;
                    end
                    % if iter_almijo == obj.armijo_maxiteration; report_flag(1) = 0x01; end    %%%% debug mode %%%%
                    
                    % Dual problem
                    XZ                      = (Xnew.*Zk); % (n,1)vector
                    minNM                   = nu/obj.M_L; % (1,1)vector
                    maxNM                   = nu*obj.M_L; % (1,1)vector
                    C_lk                    = XZ;
                    C_uk                    = XZ;
                    C_lk(XZ>minNM)          = minNM;
                    C_uk(XZ<maxNM)          = maxNM;
                    flag_positive           = dZ>0;
                    minAlpha_ng             = (C_lk./Xnew - Zk)./dZ; % dZ<=0
                    maxAlpha_ps             = (C_uk./Xnew - Zk)./dZ; % dZ>0. (dZ=0 -> Inf)
                    Alpha_yk                = minAlpha_ng;           % (n,1)vector
                    Alpha_yk(flag_positive) = maxAlpha_ps(flag_positive);
                    % alpha_zk                = min([1;Alpha_yk]);          % % (dZ=0 -> alpha_zk=1)

                    % new solution
                    Znew    = Zk + min([1;Alpha_yk])*dZ; % Znew = Zk + alpha_zk*dZ;
                    Ynew    = Yk + dY;                 % Ynew = Yk + alpha_yk*dY; (alpha_yk = 1)
                    F0new   = obj.F_(Xnew);
                    H0new   = obj.H_(Xnew);
                    Fx0new  = obj.F_X(Xnew);
                    Hx0new  = obj.H_X(Xnew);
                    Lx0new  = obj.L_X(Xnew,Ynew,Znew,Fx0new,Hx0new);
                    % BFGS updating
                    s  = Xnew - Xk;
                    y  = Lx0new - obj.L_X(Xk,Ynew,Znew,Fx,Hx);
                    B  = obj.modifiedBFGS(B,s,y');

                    % Updating
                    Xk  = Xnew;
                    Yk  = Ynew;
                    Zk  = Znew;

                    F  = F0new;
                    H  = H0new;
                    Fx = Fx0new;
                    Hx = Hx0new;
                    Lx = Lx0new;
                    rk = obj.r(Xk,Zk,Lx,H,nu);
                    dW = -inv(obj.J_k(Xk,Zk,B,Hx))*rk;

                    dX  = dW(p_xdim);
                    dY  = dW(p_ydim);
                    dZ  = dW(p_zdim);

                    nx  = -Xk./dX;
                    rho = max(abs(Yk+dY)) + 1;
                    alpha_xk = min([1, obj.gamma*min(nx(dX<0))]);
                    P    = obj.P_(Xk,F,H,nu,rho);
                    xi_dP   = obj.xi*obj.dP_(Xk,dX,Fx,H,Hx,nu,rho);
                    % fprintf("Iteration (%03d,%03d)," + ...
                    %         " (x, f(x)) = (%8f, %8f)\n", iter_outer, iter_inner, getValue(obj,Xk),F)    %%%% debug mode %%%%
                    % inner convergence condition
                    if norm(rk) <= obj.M_c*nu
                        break;
                    end
                end
                if obj.nu_tolerance > nu
                    break;
                end
                nu = nu*obj.tau;
                rk = obj.r(Xk,Zk,Lx,H,nu);
                dW       = -inv(obj.J_k(Xk,Zk,BE,Hx))*rk; % []
                dX       = dW(p_xdim);
                dY       = dW(p_ydim);
                dZ       = dW(p_zdim);
                nx       = -Xk./dX;
                rho      = max(abs(Yk+dY)) + 1;
                % Armijo condition
                alpha_xk = min([1, obj.gamma*min(nx(dX<0))]);
                P        = obj.P_(Xk,F,H,nu,rho);
                xi_dP    = obj.xi*obj.dP_(Xk,dX,Fx,H,Hx,nu,rho);
                B  = BE; % initial BFGS matrix 
                
            end
            obj.Ff   = F;
            obj.Xf   = Xk;
            obj.Yf   = Yk;
            obj.Zf   = Zk;
            % SystemCheck_Run(obj,report_flag);    %%%% debug mode %%%%
        end
        function SystemCheck_Run(obj,report_flag)
            if sum(report_flag,"all")
                if report_flag(1)
                    obj.SystemReport_string = obj.SystemReport_string ...
                        + "Iteration has reached armijo_maxiteration\n";
                end
                if report_flag(2)
                    obj.SystemReport_string = obj.SystemReport_string ...
                        + "";
                end
                if report_flag(3)
                    obj.SystemReport_string = obj.SystemReport_string ...
                        + "";
                end
                fprintf(obj.SystemReport_string)
            end
        end
    end
    methods (Static)
        function B_out = modifiedBFGS(B,s,y,omega)
            % s = x(k+1)-x(k)                                            : vector(n,1)
            % y = L_x(x(k+1),y(k+1),z(x+1))' - L_x(x(k),y(k+1),z(x+1))'  : vector(n,1)
            arguments
                B
                s
                y
                omega = 0.2;
            end
            Bs  = B*s;
            sBs = s'*Bs;
            sy  = s'*y;
            if sy < omega*sBs
                phi = (1-omega)*sBs/(sBs-sy);
            else
                phi = 1.0;
            end
            z = phi*y + (1 - phi)*Bs;
            B_out = B - Bs*(s'*B)./sBs + z*z'./s'*z;
        end
    end
end