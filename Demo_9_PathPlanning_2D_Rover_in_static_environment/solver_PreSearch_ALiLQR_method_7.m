classdef solver_PreSearch_ALiLQR_method_7 < handle
    properties
        NN
        tt
        Q_k
        q_k
        R_k
        r_k
        H_k
        Q_N
        q_N
        Lxx
        Luu
        Lxu
        Lux
        XX
        UU
        size_input
        size_state
        x_tf
        u_tf
        hf_CalcAB
        hf_DiscreteDynamics

        LLambda_t
        LLambda_tf
        con_u_dim
        ceq_u_dim
        con_x_dim
        ceq_x_dim
        con_x_tf_dim
        ceq_x_tf_dim

        C_tf        
        C_t         
        Cx_tf
        Cx_t
        Cu_t

        constraints_t_dim
        constraints_tf_dim
        con_dim
        con_tf_dim
        ceq_dim
        ceq_tf_dim

        con_pointer
        con_tf_pointer
        con_x_pointer   
        con_x_tf_pointer

        E_t    
        E_tf   
        E_t0   
        E_tf0  

        update_LM_t 
        update_LM_tf

        penalty_multipliers_t
        penalty_multipliers_tf
        penalty_gain_phi
        constraints_tolerance
        Constraint_MaxIteration
        J_Log
        Jc_Log 
        LLambda_t2
        LLambda_tf2
        con_u_dim2
        ceq_u_dim2
        con_x_dim2
        ceq_x_dim2
        con_x_tf_dim2
        ceq_x_tf_dim2

        C_tf2   
        C_t2         
        Cx_tf2
        Cx_t2
        Cu_t2
        constraints_t_dim2
        constraints_tf_dim2
        con_dim2
        con_tf_dim2
        ceq_dim2
        ceq_tf_dim2

        con_pointer2
        con_tf_pointer2
        con_x_pointer2 
        con_x_tf_pointer2

        E_t2   
        E_tf2  
        E_t02  
        E_tf02 

        update_LM_t2 
        update_LM_tf2 

        penalty_multipliers_t2
        penalty_multipliers_tf2
        penalty_gain_phi2
        constraints_tolerance2
        Constraint_MaxIteration2
        J_Log2
        Tolerance
        MaxIteration
        UserInterface
        UserInterface_end

        
    end
    methods
        function obj = solver_PreSearch_ALiLQR_method_7(Q,q,R,r,H,Q_N,q_N,Ts,Tf,CalcAB,DiscreteDynamics,Tolerance,MaxIteration)
            arguments
                Q 
                q 
                R 
                r 
                H
                Q_N
                q_N
                Ts 
                Tf 
                CalcAB 
                DiscreteDynamics
                Tolerance    = 0.1
                MaxIteration = 50
            end

            obj.tt = Ts : Ts : Tf;
            obj.NN = size(obj.tt,2);
            obj.Q_k = Q;
            obj.R_k = R;

            obj.size_input = size(R,2);
            obj.size_state = size(Q,2);
            
            if isempty(q)
                q = zeros(obj.size_state,1);
            end
            if isempty(r)
                r = zeros(obj.size_input,1);
            end
            if isempty(H)
                H = zeros(obj.size_state,obj.size_input);
            end
            if isempty(Q_N)
                Q_N = zeros(obj.size_state,obj.size_state);
            end
            if isempty(q_N)
                q_N = zeros(obj.size_state,1);
            end

            obj.q_k = q;
            obj.r_k = r;
            obj.H_k = H;
            obj.Q_N = Q_N;
            obj.q_N = q_N;

            obj.Lxx = Q;
            obj.Luu = R;
            obj.Lxu = H;
            obj.Lux = H';

            obj.UU = zeros(obj.size_input, obj.NN);
            obj.XX = zeros(obj.size_state, obj.NN+1);

            obj.hf_CalcAB           = CalcAB;
            obj.hf_DiscreteDynamics = DiscreteDynamics;

            obj.Tolerance    = Tolerance;
            obj.MaxIteration = MaxIteration;
        end
        function setconstraint(obj,con_u,ceq_u,GradCon_u,GradCeq_u,...
                                con_x,ceq_x,GradCon_x,GradCeq_x,...
                                con_x_tf,ceq_x_tf,GradCon_x_tf,GradCeq_x_tf,...
                                option)
            arguments
                obj 
                con_u 
                ceq_u 
                GradCon_u 
                GradCeq_u 
                con_x
                ceq_x
                GradCon_x 
                GradCeq_x 
                con_x_tf 
                ceq_x_tf 
                GradCon_x_tf 
                GradCeq_x_tf 
                option.Tolerance = 1e-6
                option.PenaltyMultipliers_t    = 20
                option.PenaltyMultipliers_tf   = 20
                option.Constraint_MaxIteration = 10;
                option.penalty_gain_phi        = 2.2;
            end
            obj.constraints_tolerance = option.Tolerance;
            obj.penalty_multipliers_t = option.PenaltyMultipliers_t;
            obj.penalty_multipliers_tf = option.PenaltyMultipliers_tf;

            if isempty(con_u)
                con_u = @(x,u)[];
                GradCon_u = @(x,u)[];
                con_u_dim_ = 0;
            else
                con_u_dim_ = size(con_u(obj.XX,obj.UU),1);
            end
            if isempty(ceq_u)
                ceq_u = @(x,u)[];
                GradCeq_u = @(x,u)[];
                ceq_u_dim_ = 0;
            else
                ceq_u_dim_ = size(ceq_u(obj.XX,obj.UU),1);
            end
            if isempty(con_x)
                con_x = @(x,u)[];
                GradCon_x = @(x,u)[];
                con_x_dim_ = 0;
            else
                con_x_dim_ = size(con_x(obj.XX,obj.UU),1);
            end
            if isempty(ceq_x)
                ceq_x = @(x,u)[];
                GradCeq_x = @(x,u)[];
                ceq_x_dim_ = 0;
            else
                ceq_x_dim_ = size(ceq_x(obj.XX,obj.UU),1);
            end
            if isempty(con_x_tf)
                con_x_tf = @(x)[];
                GradCon_x_tf = @(x)[];
                con_x_tf_dim_ = 0;
            else
                con_x_tf_dim_ = size(con_x_tf(obj.XX(:,end)),1);
            end
            if isempty(ceq_x_tf)
                ceq_x_tf = @(x)[];
                GradCeq_x_tf = @(x)[];
                ceq_x_tf_dim_ = 0;
            else
                ceq_x_tf_dim_ = size(ceq_x_tf(obj.XX(:,end)),1);
            end
            obj.C_tf = @(x_tf)  [ceq_x_tf(x_tf); con_x_tf(x_tf)];
            obj.C_t  = @(XX,UU) [ceq_x(XX,UU);ceq_u(XX,UU);
                                 con_x(XX,UU);con_u(XX,UU);];

            obj.Cx_tf = @(x_tf)  [GradCeq_x_tf(x_tf); GradCon_x_tf(x_tf)];
            obj.Cx_t  = @(XX,UU)...
                [GradCeq_x(XX,UU);zeros(ceq_u_dim_,obj.size_state,obj.NN);
                 GradCon_x(XX,UU);zeros(con_u_dim_,obj.size_state,obj.NN)];
            
            obj.Cu_t  = @(XX,UU)...
                [zeros(ceq_x_dim_,obj.size_input,obj.NN); GradCeq_u(XX,UU);
                 zeros(con_x_dim_,obj.size_input,obj.NN); GradCon_u(XX,UU)];

            obj.con_u_dim = con_u_dim_;
            obj.ceq_u_dim = ceq_u_dim_;
            obj.con_x_dim = con_x_dim_;
            obj.ceq_x_dim = ceq_x_dim_;
            obj.con_x_tf_dim = con_x_tf_dim_;
            obj.ceq_x_tf_dim = ceq_x_tf_dim_;

            obj.constraints_t_dim  = con_u_dim_ + ceq_u_dim_ + con_x_dim_ + ceq_x_dim_;
            obj.constraints_tf_dim = con_x_tf_dim_ + ceq_x_tf_dim_;
            obj.con_dim    = con_u_dim_ + con_x_dim_;
            obj.con_tf_dim = con_x_tf_dim_;
            obj.ceq_dim    = ceq_u_dim_ + ceq_x_dim_;
            obj.ceq_tf_dim = ceq_x_tf_dim_;
            obj.E_t   = ones(obj.constraints_t_dim,obj.NN).*option.PenaltyMultipliers_t;
            obj.E_tf  = ones(obj.constraints_tf_dim,1).*option.PenaltyMultipliers_tf;
            obj.E_t0  = zeros(obj.constraints_t_dim,obj.NN);
            obj.E_tf0 = zeros(obj.constraints_tf_dim,1);
            if obj.con_dim
                obj.update_LM_t = @(lambda_t,Ic)...
                obj.static_update_Lagrange_multipliers(lambda_t,Ic,...
                @(lambda)obj.static_update_Lagrange_multipliers_con(lambda,...
                    (obj.ceq_dim+1):(obj.ceq_dim+obj.con_dim)) );
            elseif obj.ceq_dim
                obj.update_LM_t = @(lambda_t,Ic)...
                obj.static_update_Lagrange_multipliers(lambda_t,Ic,...
                @(x)x);
            else
                obj.update_LM_t  = @(lambda,I) lambda;
            end
            if obj.con_tf_dim
                obj.update_LM_tf = @(lambda_tf,Ic)...
                obj.static_update_Lagrange_multipliers(lambda_tf,Ic,...
                @(lambda)obj.static_update_Lagrange_multipliers_con(lambda,...
                    (obj.ceq_tf_dim+1):(obj.ceq_tf_dim+obj.con_tf_dim)) );
            elseif obj.ceq_tf_dim
                obj.update_LM_tf = @(lambda_tf,Ic)...
                obj.static_update_Lagrange_multipliers(lambda_tf,Ic,...
                @(x)x);
            else
                obj.update_LM_tf  = @(lambda,I) lambda;
            end
            if obj.constraints_t_dim == 0
                res0   = zeros(1,obj.NN);
                resx00 = zeros(1,obj.size_state,obj.NN);
                resu00 = zeros(1,obj.size_input,obj.NN);
                obj.C_t     = @(XX,UU) res0;
                obj.Cx_t    = @(XX,UU) resx00;
                obj.Cu_t    = @(XX,UU) resu00;
                obj.constraints_t_dim = 1;
                obj.con_dim = 1;
                obj.ceq_dim = 0;
                obj.E_t   = zeros(1,obj.NN);
                obj.E_t0  = zeros(1,obj.NN);
            end
            if obj.constraints_tf_dim == 0
                res0   = 0;
                resx00 = zeros(1,obj.size_state);
                obj.C_tf  = @(XX) res0;
                obj.Cx_tf = @(XX) resx00;
                obj.constraints_tf_dim = 1;
                obj.con_tf_dim = 1;
                obj.ceq_tf_dim = 0;
                obj.E_tf   = 0;
                obj.E_tf0  = 0;
            end
            if (obj.constraints_tf_dim + obj.constraints_t_dim) > 0
                obj.Constraint_MaxIteration = option.Constraint_MaxIteration;
                obj.penalty_gain_phi = option.penalty_gain_phi;
            else
                obj.Constraint_MaxIteration = 1;
                obj.penalty_gain_phi = 0;
            end
            obj.LLambda_t    = NaN(size(obj.E_t0));
            obj.J_Log        = NaN(obj.Constraint_MaxIteration,obj.MaxIteration+1);
            obj.con_pointer = (obj.ceq_dim+1):(obj.ceq_dim+obj.con_dim);
            obj.con_tf_pointer = (obj.ceq_tf_dim+1):(obj.ceq_tf_dim+obj.con_tf_dim);
            obj.con_x_pointer    = (obj.ceq_dim+1):(obj.ceq_dim+con_x_dim_);
            obj.con_x_tf_pointer = (obj.ceq_tf_dim+1):(obj.ceq_tf_dim+obj.con_tf_dim);
        end
        function setconstraint2(obj,con_u,ceq_u,GradCon_u,GradCeq_u,...
                                con_x,ceq_x,GradCon_x,GradCeq_x,...
                                con_x_tf,ceq_x_tf,GradCon_x_tf,GradCeq_x_tf,...
                                option)
            arguments
                obj 
                con_u 
                ceq_u 
                GradCon_u 
                GradCeq_u 
                con_x
                ceq_x
                GradCon_x 
                GradCeq_x 
                con_x_tf 
                ceq_x_tf 
                GradCon_x_tf 
                GradCeq_x_tf 
                option.Tolerance = 1e-6
                option.PenaltyMultipliers_t    = 20
                option.PenaltyMultipliers_tf   = 20
                option.Constraint_MaxIteration = 10;
                option.penalty_gain_phi        = 2.2;
                % option.Constraint_MaxIteration = 15;
                % option.penalty_gain_phi        = 1.75;
            end
            obj.constraints_tolerance2 = option.Tolerance;
            obj.penalty_multipliers_t2 = option.PenaltyMultipliers_t;
            obj.penalty_multipliers_tf2 = option.PenaltyMultipliers_tf;

            if isempty(con_u)
                con_u = @(x,u)[];
                GradCon_u = @(x,u)[];
                con_u_dim_ = 0;
            else
                con_u_dim_ = size(con_u(obj.XX,obj.UU),1);
            end
            if isempty(ceq_u)
                ceq_u = @(x,u)[];
                GradCeq_u = @(x,u)[];
                ceq_u_dim_ = 0;
            else
                ceq_u_dim_ = size(ceq_u(obj.XX,obj.UU),1);
            end
            if isempty(con_x)
                con_x = @(x,u)[];
                GradCon_x = @(x,u)[];
                con_x_dim_ = 0;
            else
                con_x_dim_ = size(con_x(obj.XX,obj.UU),1);
            end
            if isempty(ceq_x)
                ceq_x = @(x,u)[];
                GradCeq_x = @(x,u)[];
                ceq_x_dim_ = 0;
            else
                ceq_x_dim_ = size(ceq_x(obj.XX,obj.UU),1);
            end
            if isempty(con_x_tf)
                con_x_tf = @(x)[];
                GradCon_x_tf = @(x)[];
                con_x_tf_dim_ = 0;
            else
                con_x_tf_dim_ = size(con_x_tf(obj.XX(:,end)),1);
            end
            if isempty(ceq_x_tf)
                ceq_x_tf = @(x)[];
                GradCeq_x_tf = @(x)[];
                ceq_x_tf_dim_ = 0;
            else
                ceq_x_tf_dim_ = size(ceq_x(obj.XX(:,end)),1);
            end
            obj.C_tf2 = @(x_tf)  [ceq_x_tf(x_tf); con_x_tf(x_tf)];
            obj.C_t2  = @(XX,UU) [ceq_x(XX,UU);ceq_u(XX,UU);
                                 con_x(XX,UU);con_u(XX,UU);];
            obj.Cx_tf2 = @(x_tf)  [GradCeq_x_tf(x_tf); GradCon_x_tf(x_tf)];
            obj.Cx_t2  = @(XX,UU)...
                [GradCeq_x(XX,UU);zeros(ceq_u_dim_,obj.size_state,obj.NN);
                 GradCon_x(XX,UU);zeros(con_u_dim_,obj.size_state,obj.NN)];
            obj.Cu_t2  = @(XX,UU)...
                [zeros(ceq_x_dim_,obj.size_input,obj.NN); GradCeq_u(XX,UU);
                 zeros(con_x_dim_,obj.size_input,obj.NN); GradCon_u(XX,UU)];

            obj.con_u_dim2 = con_u_dim_;
            obj.ceq_u_dim2 = ceq_u_dim_;
            obj.con_x_dim2 = con_x_dim_;
            obj.ceq_x_dim2 = ceq_x_dim_;
            obj.con_x_tf_dim2 = con_x_tf_dim_;
            obj.ceq_x_tf_dim2 = ceq_x_tf_dim_;
            obj.constraints_t_dim2  = con_u_dim_ + ceq_u_dim_ + con_x_dim_ + ceq_x_dim_;
            obj.constraints_tf_dim2 = con_x_tf_dim_ + ceq_x_tf_dim_;
            obj.con_dim2    = con_u_dim_ + con_x_dim_;
            obj.con_tf_dim2 = con_x_tf_dim_;
            obj.ceq_dim2    = ceq_u_dim_ + ceq_x_dim_;
            obj.ceq_tf_dim2 = ceq_x_tf_dim_;
            obj.E_t2   = ones(obj.constraints_t_dim2,obj.NN).*option.PenaltyMultipliers_t;
            obj.E_tf2  = ones(obj.constraints_tf_dim2,1).*option.PenaltyMultipliers_tf;
            obj.E_t02  = zeros(obj.constraints_t_dim2,obj.NN);
            obj.E_tf02 = zeros(obj.constraints_tf_dim2,1);
            if obj.con_dim2
                obj.update_LM_t2 = @(lambda_t,Ic)...
                obj.static_update_Lagrange_multipliers(lambda_t,Ic,...
                @(lambda)obj.static_update_Lagrange_multipliers_con(lambda,...
                    (obj.ceq_dim2+1):(obj.ceq_dim2+obj.con_dim2)) );
            elseif obj.ceq_dim2
                obj.update_LM_t2 = @(lambda_t,Ic)...
                obj.static_update_Lagrange_multipliers(lambda_t,Ic,...
                @(x)x);
            else
                obj.update_LM_t2  = @(lambda,I) lambda;
            end
            if obj.con_tf_dim2
                obj.update_LM_tf2 = @(lambda_tf,Ic)...
                obj.static_update_Lagrange_multipliers(lambda_tf,Ic,...
                @(lambda)obj.static_update_Lagrange_multipliers_con(lambda,...
                    (obj.ceq_tf_dim2+1):(obj.ceq_tf_dim2+obj.con_tf_dim2)) );
            elseif obj.ceq_tf_dim2
                obj.update_LM_tf2 = @(lambda_tf,Ic)...
                obj.static_update_Lagrange_multipliers(lambda_tf,Ic,...
                @(x)x);
            else
                obj.update_LM_tf2  = @(lambda,I) lambda;
            end
            if obj.constraints_t_dim2 == 0
                res0   = zeros(1,obj.NN);
                resx00 = zeros(1,obj.size_state,obj.NN);
                resu00 = zeros(1,obj.size_input,obj.NN);
                obj.C_t2     = @(XX,UU) res0;
                obj.Cx_t2    = @(XX,UU) resx00;
                obj.Cu_t2    = @(XX,UU) resu00;
                obj.constraints_t_dim2 = 1;
                obj.con_dim2 = 1;
                obj.ceq_dim2 = 0;
                obj.E_t2   = zeros(1,obj.NN);
                obj.E_t02  = zeros(1,obj.NN);
            end
            if obj.constraints_tf_dim2 == 0
                res0   = 0;
                resx00 = zeros(1,obj.size_state);
                obj.C_tf2  = @(XX) res0;
                obj.Cx_tf2 = @(XX) resx00;
                obj.constraints_tf_dim2 = 1;
                obj.con_tf_dim2 = 1;
                obj.ceq_tf_dim2 = 0;
                obj.E_tf2   = 0;
                obj.E_tf02  = 0;
            end
            if (obj.constraints_tf_dim2 + obj.constraints_t_dim2) > 0
                obj.Constraint_MaxIteration2 = option.Constraint_MaxIteration;
                obj.penalty_gain_phi2 = option.penalty_gain_phi;
            else
                obj.Constraint_MaxIteration2 = 1;
                obj.penalty_gain_phi2 = 0;
            end
            obj.LLambda_t2    = NaN(size(obj.E_t02));
            obj.J_Log2        = NaN(obj.Constraint_MaxIteration2,obj.MaxIteration+1);
            obj.con_pointer2 = (obj.ceq_dim2+1):(obj.ceq_dim2+obj.con_dim2);
            obj.con_tf_pointer2 = (obj.ceq_tf_dim2+1):(obj.ceq_tf_dim2+obj.con_tf_dim2);
            obj.con_x_pointer2    = (obj.ceq_dim2+1):(obj.ceq_dim2+con_x_dim_);
            obj.con_x_tf_pointer2 = (obj.ceq_tf_dim2+1):(obj.ceq_tf_dim2+obj.con_tf_dim2);
        end
        function setUserInterface(obj,func)
            if isempty(func)
                obj.UserInterface = [];
            else
                obj.UserInterface = func;
            end
        end
        function setUserInterface_end(obj,func)
            if isempty(func)
                obj.UserInterface_end = [];
            else
                obj.UserInterface_end = func;
            end
        end
        function updateTerminalState(obj,x_tf)
            obj.x_tf = x_tf;
        end
        function updateTerminalValues(obj,x_tf,u_tf)
            obj.x_tf = x_tf;
            obj.u_tf = u_tf;
        end
        function set_InitialGuess(obj,uu)
            obj.UU = uu;
        end
        function set_InitialState(obj,xx)
            obj.XX = xx;
        end
        function set_InitialGuess_RepMat(obj,u)
            obj.UU = repmat(u,1,obj.NN);
        end

        function Run(obj,x0)
            lambda_t  = obj.E_t0;
            lambda_tf = obj.E_tf0;
            I_t0      = obj.E_t;
            I_tf0     = obj.E_tf;
            phi       = obj.penalty_gain_phi;
            lambda_t2  = obj.E_t02;
            lambda_tf2 = obj.E_tf02;
            I_t02      = obj.E_t2;
            I_tf02     = obj.E_tf2;
            phi2       = obj.penalty_gain_phi2;
            Jlog = NaN(obj.Constraint_MaxIteration,obj.NN+1);
            Jclog = NaN(obj.Constraint_MaxIteration,obj.NN+1);
            uu = obj.UU;
            xx = zeros(obj.size_state,obj.NN+1);
            xx(:,1) = x0;
            for i = 1 : obj.NN
                xx(:,i+1) = obj.hf_DiscreteDynamics(xx(:,i),uu(:,i));
            end
            [c_t, c_tf, cx_t, cx_tf, cu_t] = obj.constraints_update(xx,uu);
            [I_t, I_tf]  = obj.activation_constraints(c_t,c_tf,lambda_t,lambda_tf,I_t0,I_tf0);
            [c_t2, c_tf2, cx_t2, cx_tf2, cu_t2] = obj.constraints_update2(xx,uu);
            [I_t2, I_tf2]  = obj.activation_constraints2(c_t2,c_tf2,lambda_t2,lambda_tf2,I_t02,I_tf02);
            [J_cost,J_constraint,Ic_t,Ic_tf,Ic_t2,Ic_tf2] = obj.Cost(...
                xx,uu,  I_t,I_tf,c_t,c_tf,lambda_t,lambda_tf,I_t2,I_tf2,c_t2,c_tf2,lambda_t2,lambda_tf2);
            Iter_con = 1;
            J_pre     = J_cost + J_constraint;
            Jlog(1,1) = J_cost;
            Jclog(1,1) = J_constraint;
            while true
                Iter  = 1;
                while true 
                    [KK,dd,dV] = obj.Backward(xx,uu,...
                        I_t,I_tf,Ic_t,Ic_tf,cx_t,cx_tf,cu_t,lambda_t,lambda_tf,...
                        I_t2,I_tf2,Ic_t2,Ic_tf2,cx_t2,cx_tf2,cu_t2,lambda_t2,lambda_tf2);
                    [xx,uu,J_cost,J_constraint, I_t,I_tf, Ic_t,Ic_tf,cx_t,cx_tf,cu_t,c_t,...
                        I_t2,I_tf2, Ic_t2,Ic_tf2,cx_t2,cx_tf2,cu_t2,c_t2] ...
                        = obj.Forward(xx,uu,KK,dd,dV,J_pre, I_t0,I_tf0,lambda_t,lambda_tf,...
                        I_t02,I_tf02,lambda_t2,lambda_tf2);
                    J_curr = J_cost + J_constraint;
                    if ~isreal(J_curr)
                        fprintf("J_curr error\n");
                    end
                    Jlog(Iter_con, Iter + 1) = J_cost;
                    Jclog(Iter_con, Iter + 1) = J_constraint;
                    if (abs(J_curr - J_pre) < obj.Tolerance) || Iter == obj.MaxIteration
                        break;
                    end
                    J_pre = J_curr;
                    Iter = Iter + 1;
                end
                if Iter_con == obj.Constraint_MaxIteration
                    break;
                end
                lambda_t  = obj.update_LM_t(lambda_t,Ic_t);
                lambda_tf = obj.update_LM_tf(lambda_tf,Ic_tf);
                lambda_t2  = obj.update_LM_t2(lambda_t2,Ic_t2);
                lambda_tf2 = obj.update_LM_tf2(lambda_tf2,Ic_tf2);
                I_t0      = I_t0 .*phi;
                I_tf0     = I_tf0.*phi;
                I_t02      = I_t02 .*phi2;
                I_tf02     = I_tf02.*phi2;

                Iter_con = Iter_con + 1;
            end
            obj.LLambda_t   = lambda_t;
            obj.LLambda_tf  = lambda_tf;
            obj.LLambda_t2  = lambda_t2;
            obj.LLambda_tf2 = lambda_tf2;
            obj.XX = xx;
            obj.UU = uu;
            obj.J_Log  = Jlog;
            obj.Jc_Log = Jclog;
        end
        function [KK,dd,dV] = Backward(obj,xx,uu, ...
                    I_t,I_tf,Ic_t,Ic_tf,cx_t,cx_tf,cu_t,lambda_t,lambda_tf,...
                    I_t2,I_tf2,Ic_t2,Ic_tf2,cx_t2,cx_tf2,cu_t2,lambda_t2,lambda_tf2)
            Num = obj.NN;
            Qk = obj.Q_k;
            Rk = obj.R_k;
            xtf = obj.x_tf;
            utf = obj.u_tf;
            Vxx0 = obj.Q_N + cx_tf'*(I_tf.*cx_tf) + cx_tf2'*(I_tf2.*cx_tf2);
            Vx0  = obj.Q_N*(xx(:,end) - xtf) + cx_tf'*(Ic_tf + lambda_tf) + cx_tf2'*(Ic_tf2 + lambda_tf2);
            Vxx  = Vxx0;
            Vx   = Vx0;
            KK  = zeros(obj.size_input,obj.size_state,Num);
            dd = zeros(obj.size_input,Num);
            dV = [0,0];
            FirstOrder = 0;
            SecondOrder = 0;
            rho = 0;
            Eye_input = eye(obj.size_input);
            while true
                for i = Num : -1 : 1
                    x_k = xx(:,i);
                    u_k = uu(:,i);
                    cx     = cx_t(:,:,i);
                    cu     = cu_t(:,:,i);
                    Ic     = Ic_t(:,i);
                    I      = I_t(:,i);
                    lambda = lambda_t(:,i);
                    L_Ic   = lambda + Ic;
                    cx2     = cx_t2(:,:,i);
                    cu2     = cu_t2(:,:,i);
                    Ic2     = Ic_t2(:,i);
                    I2      = I_t2(:,i);
                    lambda2 = lambda_t2(:,i);
                    L_Ic2   = lambda2 + Ic2;

                    [A,B]  = obj.hf_CalcAB(x_k,u_k);
                    Lx = Qk*(x_k - xtf);
                    Lu = Rk*(u_k - utf);
                    Qxx = obj.Lxx + A'*Vxx*A + cx'*(I.*cx) + cx2'*(I2.*cx2);
                    Quu = obj.Luu + B'*Vxx*B + cu'*(I.*cu) + cu2'*(I2.*cu2);
                    Qux = obj.Lux + B'*Vxx*A + cu'*(I.*cx) + cu2'*(I2.*cx2);
                    Qx  = Lx + A'*Vx + cx'*L_Ic + cx2'*L_Ic2;
                    Qu  = Lu + B'*Vx + cu'*L_Ic + cu2'*L_Ic2;
                    [flag,id] = Cholesky(Quu);
                    if flag
                        Vxx = Vxx0;
                        Vx  = Vx0;
                        FirstOrder  = 0;
                        SecondOrder = 0;
                        rho = rho + abs(Quu(id,id))*10;
                        fprintf("iLQR::Backward::Cholesky, Quu is not positive definite.\n")
                        break;
                    end
                    inv_Quu_reg = inv(Quu + rho*Eye_input);
                    K = -inv_Quu_reg*Qux;
                    d = -inv_Quu_reg*Qu;
                    % update
                    Vxx = Qxx + K'*Quu*K + K'*Qux + Qux'*K;
                    Vx  = Qx  + K'*Quu*d + K'*Qu + Qux'*d;
                    FirstOrder = FirstOrder + d'*Qu;
                    SecondOrder = SecondOrder + d'*Quu*d;
                    KK(:,:,i) = K;
                    dd(:,i)   = d;
                end
                if i==1;break;end
            end
            dV(1) = FirstOrder;
            dV(2) = 0.5*SecondOrder;
        end
        function [xx,uu,J_cost,J_constraint, I_t,I_tf, Ic_t,Ic_tf,cx_t,cx_tf,cu_t,c_t,...
                I_t2,I_tf2, Ic_t2,Ic_tf2,cx_t2,cx_tf2,cu_t2,c_t2]...
                = Forward(obj,xx,uu,KK,dd,dV,J_pre, I_t0,I_tf0,lambda_t,lambda_tf,...
                I_t02,I_tf02,lambda_t2,lambda_tf2)
            dV1 = dV(1);
            dV2 = dV(2);
            Alpha = 1;
            Num = obj.NN;
            xx_ = zeros(obj.size_state,Num + 1);
            uu_ = zeros(obj.size_input,Num);
            xx_(:,1) = xx(:,1);
            while true
                for i = 1 : Num
                    xx_0  = xx_(:,i);
                    u_bar = uu(:,i) + KK(:,:,i)*(xx_0 - xx(:,i)) ...
                           + dd(:,i)*Alpha;
                    xx_(:,i+1) = obj.hf_DiscreteDynamics(xx_0,u_bar);
                    uu_(:,i)  = u_bar;
                end
                [c_t, c_tf, cx_t, cx_tf, cu_t] = obj.constraints_update(xx_,uu_);
                [I_t, I_tf]  = obj.activation_constraints(c_t,c_tf,lambda_t,lambda_tf,I_t0,I_tf0);
                [c_t2, c_tf2, cx_t2, cx_tf2, cu_t2] = obj.constraints_update2(xx_,uu_);
                [I_t2, I_tf2]  = obj.activation_constraints2(c_t2,c_tf2,lambda_t2,lambda_tf2,I_t02,I_tf02);
                [J_cost,J_constraint,Ic_t,Ic_tf,Ic_t2,Ic_tf2]  = obj.Cost(xx_,uu_, I_t,I_tf,c_t,c_tf,lambda_t,lambda_tf,...
                    I_t2,I_tf2,c_t2,c_tf2,lambda_t2,lambda_tf2);
                J_curr = J_cost + J_constraint;
                z = (J_curr - J_pre)/(Alpha*dV1 + Alpha*Alpha*dV2);
                if (1e-4<=z)&&(z<=10)
                    xx = xx_;
                    uu = uu_;
                    break;
                end
                if Alpha<(1e-4)
                    xx = xx_;
                    uu = uu_;
                    break;
                end
                Alpha = 0.5*Alpha;
            end
        end
        function [J,J_c, Ic_t,Ic_tf, Ic_t2,Ic_tf2] = Cost(obj,xx,uu, I_t,I_tf,c_t,c_tf,lambda_t,lambda_tf...
                , I_t2,I_tf2,c_t2,c_tf2,lambda_t2,lambda_tf2) 

            dxx = xx - obj.x_tf;
            duu = uu - obj.u_tf;
            Q = obj.Q_k;
            q = obj.q_k;
            R = obj.R_k;
            r = obj.r_k;
            H = obj.H_k;
            Ic_t   = c_t.*I_t;
            Ic_tf  = c_tf.*I_tf;
            Ic_t2   = c_t2.*I_t2;
            Ic_tf2  = c_tf2.*I_tf2;
            cost_c1 = c_tf'*Ic_tf + c_tf2'*Ic_tf2;
            cost_c2 = lambda_tf'*c_tf + lambda_tf2'*c_tf2;
            cost1 = dxx(:,end)'*obj.Q_N*dxx(:,end);
            cost2 = obj.q_N'*dxx(:,end);
            for i = 1 : obj.NN
                dx = dxx(:,i);
                du = duu(:,i);
                cost_c1 = cost_c1 + c_t(:,i)'*Ic_t(:,i) + c_t2(:,i)'*Ic_t2(:,i);
                cost_c2 = cost_c2 + lambda_t(:,i)'*c_t(:,i)+lambda_t2(:,i)'*c_t2(:,i);
                cost1 = cost1 + dx'*Q*dx + du'*R*du;
                cost2 = cost2 + dx'*H*du + q'*dx + r'*du;
            end
            J_c = cost_c1*0.5 + cost_c2;
            J   = cost1*0.5 + cost2;

        end
        function [I_t, I_tf] = activation_constraints(obj,cc_t,cc_tf,lambda_t,lambda_tf,I_t,I_tf)
            id_con = obj.con_pointer;
            Icon = ~( (cc_t(id_con,:)<=0).*(lambda_t(id_con,:)==0));
            I_t(id_con,:) = I_t(id_con,:).*Icon;
            id_con_tf = obj.con_tf_pointer;
            Icon_tf = ~((cc_tf(id_con_tf,:)<=0).*(lambda_tf(id_con_tf,:)==0));
            I_tf(id_con_tf,:) = I_tf(id_con_tf,:).*Icon_tf;
        end
        function [I_t, I_tf] = activation_constraints2(obj,cc_t,cc_tf,lambda_t,lambda_tf,I_t,I_tf)
            id_con = obj.con_pointer2;
            Icon = ~( (cc_t(id_con,:)<=0).*(lambda_t(id_con,:)==0));
            I_t(id_con,:) = I_t(id_con,:).*Icon;
            id_con_tf = obj.con_tf_pointer2;
            Icon_tf = ~((cc_tf(id_con_tf,:)<=0).*(lambda_tf(id_con_tf,:)==0));
            I_tf(id_con_tf,:) = I_tf(id_con_tf,:).*Icon_tf;
        end
        function [c_t, c_tf, cx_t, cx_tf, cu_t] = constraints_update(obj,xx,uu)
            c_t   = obj.C_t  (xx,uu);    
            c_tf  = obj.C_tf (xx(:,end));
            cx_t  = obj.Cx_t (xx,uu);    
            cx_tf = obj.Cx_tf(xx(:,end));
            cu_t  = obj.Cu_t (xx,uu);    
        end
        function [c_t, c_tf, cx_t, cx_tf, cu_t] = constraints_update2(obj,xx,uu)
            c_t   = obj.C_t2  (xx,uu);     
            c_tf  = obj.C_tf2 (xx(:,end)); 
            cx_t  = obj.Cx_t2 (xx,uu);     
            cx_tf = obj.Cx_tf2(xx(:,end)); 
            cu_t  = obj.Cu_t2 (xx,uu);     
        end
    end
    methods (Static)
        function lambda_t...
                = static_update_Lagrange_multipliers(lambda_t,Ic,Lag_mult_con)
            lambda_t  = lambda_t  + Ic;
            lambda_t  = Lag_mult_con(lambda_t);
        end
        function lambda = static_update_Lagrange_multipliers_con(lambda,pointer)
            lambda(pointer,:) = max(lambda(pointer,:),0);
        end
    end
end