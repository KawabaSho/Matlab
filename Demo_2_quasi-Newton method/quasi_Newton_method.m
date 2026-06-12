classdef quasi_Newton_method < handle
    %% quasi-Newton method for dot_x = f(x)
    %
    % Calculate x_opt given f(x_opt) = 0
    %
    properties
        gradFunc

        Tolerance
        MaxIteration
        alpha_d
        messsage

        core_run
    end

    methods
        function obj = quasi_Newton_method(f,option)
            % Calculate x_opt given f(x_opt) = 0
            % option.alpha : x_{k+1} = x_{k} + alpha*d_{k}
            arguments
                f                                              = [];
                option.Tolerance (1,1) double {mustBePositive} = 1e-6
                option.MaxIteration (1,1) {mustBeInteger}      = 1.5*1e+7
                option.alpha_d                                 = 1
                option.message                                 = true
            end
            obj.gradFunc = f;
            obj.Tolerance    = option.Tolerance;
            obj.MaxIteration = option.MaxIteration;
            obj.alpha_d        = option.alpha_d;
            obj.messsage     = option.message;
            obj.core_run = @(x)obj.solver(x,option.MaxIteration,option.Tolerance,f,option.message);

        end
        function xopt = solve(obj,x)
            xopt = obj.core_run(x);
        end
    end
    methods (Static)
        function xopt = solver(x,MaxIteration,Tolerance,gradf,flag)
            x0  = x;
            x1  = x;
            dim = numel(x);
            II  = eye(dim);
            H   = II; % inv(B0)

            
            for Iter = 1 : MaxIteration
                gradfx0 = gradf(x0);
                if max(abs(gradfx0),[],"all")<Tolerance;break;end

                d0 = -H*gradfx0; % d0 = s0

                x1 = x0 + d0;
                s0 = x1 - x0;
                y0 = gradf(x1) - gradfx0;
                
                sy = s0'*y0;
                H = (II - s0*y0'/sy)*H*(II - y0*s0'/sy) + s0*s0'/sy;
            end
            if flag
                if Iter == MaxIteration
                    fprintf(' quasi_Newton_method:solve : Iter is MaxIteration.\n')
                end
            end
            xopt = x1;
        end
        % function a_opt = Line_Search(~,f,x,d,gradf,gradfx)
        %         alpha = 0 : 0.001 : 1;  
        %         num = numel(alpha); 
        %         c1 = 1e-4;
        %         c2 = 0.9;
        %         for i = 2 : num
        %             a = alpha(i);
        %             if (f(x+a*d) <= f(x) + c1+a*gradfx'*d)
        %                 if (gradf(x+a*d)'*d >= c2*gradfx'*d)
        %                     a_opt = alpha(i);
        %                     break;
        %                 end
        %             end
        %         end
        %         if i == num;a_opt = 1;end
        % end
    end
end