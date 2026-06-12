classdef class_mppi < handle
    % MPPIを使た最適化です．
    % sample数（K）が10000以上はParallel Computing Toolbox で計算した方が速いです．
    % 

    properties
        alpha
        K
        T
        Cov
        Cov_
        L
        lambda
        index1
        index2

        clamp
        phi_tf
        phi_t
        Cost
        Func

        dim_u
        dim_x
        tra2u_opt;
        tra2x_opt;
        tra2data_u;
        tra2time_u;
        tra2data_x;
        tra2time_x;

        U0
        U_opt
        X_opt
        Cost_opt

        UU_sample
        XX_sample
        Cost_sample
    end
    properties (Constant)
        gamma = 10;
    end

    methods
        function obj = class_mppi(Cost_tf,Cost_t,T,K,U0,Cov,lambda,SystemDynamics,dim_states,options)
            % Cost(x(t))
            % T is the number of timesteps.
            % K is the number of samples.
            % U0 is the initial guess: [u(0),u(1),...,u(T-1)].
            % Cov is the covariance matrix.
            % Lambda is temperature.
            % SystemDynamics is x(t+1) = SystemDynamics(x(t),u(t))
            % x(t) = [x_1, x_2, ..., x_dim_states]
            % dim_states: x(t) in R^(dim_states)
            % X_opt = [x(0),x(1),...,x(T)]
            arguments
                Cost_tf function_handle
                Cost_t function_handle
                T 
                K 
                U0 
                Cov 
                lambda
                SystemDynamics
                dim_states
                options.umax = [];
                options.umin = [];
            end

            al = 1 - obj.gamma/lambda;
            % if ~((0<=al)&&(al<=1))
            %     error("class_mppi.m クラス生成に失敗しました. \n0<=α<=1を満たしていません．\n")
            % end
            obj.alpha  = al;
            obj.K      = K;
            obj.T      = T;
            obj.Cov    = Cov;
            obj.Cov_   = inv(Cov);
            obj.lambda = lambda;
            
            obj.phi_tf = Cost_tf;
            obj.phi_t  = Cost_t;
            obj.Cost   = @(x)cost(x,Cost_tf,Cost_t,T);
            obj.Func   = SystemDynamics;
    
            dimu  = size(U0,1);
            dimx  = dim_states;
            obj.dim_u  = dimu;
            obj.dim_x  = dimx;

            obj.U_opt  = reshape(U0,[1,dimu*T]);
            obj.U0     = reshape(U0,[1,dimu*T]);
            % [u_1(0),u_2(0),...,u_dimu(0), u_1(1),u_2(1),...,u_dimu(1),...]
            obj.X_opt  = zeros(1,dimx*(T+1));
            obj.Cost_opt = Inf;
    
            % (K,dimu,T) -> reshape -> (K,dim*T)
            tra2u_opt = [dimu,T];
            tra2x_opt = [dimx,T+1];
            tra2data_u = [K,dimu,T];
            tra2time_u = [K,dimu*T];
            tra2data_x = [K,dimx,T+1];
            tra2time_x = [K,dimx*(T+1)];
            obj.tra2u_opt  = tra2u_opt;
            obj.tra2x_opt  = tra2x_opt;
            obj.tra2data_u = tra2data_u;
            obj.tra2time_u = tra2time_u;
            obj.tra2data_x = tra2data_x;
            obj.tra2time_x = tra2time_x;

            obj.UU_sample   = NaN(tra2time_u);
            obj.XX_sample   = NaN(tra2time_x);
            obj.Cost_sample = NaN(K,1);

            % for sampling
            obj.L = chol(Cov);
            
            % clamping
            umin = options.umin; % umin = [u1,u2,u3,...]
            umax = options.umax;
            if isempty(umax)
                if isempty(umin)
                    obj.clamp = @(x)x;
                else
                    obj.clamp = @(x)obj.clamp_min(x,umin,T,dimu);
                end
            else
                if isempty(umin)
                    obj.clamp = @(x)clamp_max(x,umax,T,dimu);
                else
                    obj.clamp = @(x)obj.clamp_minmax(x,umin,umax,T,dimu);
                end
            end

            % index
            % al
            g1         = fix((1-al)*K);
            obj.index1 = 1 : g1;
            obj.index2 = (g1+1) : K;

            % Parallel Computing setting
            % parpool
        end

        function  u0_opt = run(obj,x0)
            K_      = obj.K;
            T_      = obj.T;
            lambda_ = obj.lambda;
            dim_u_  = obj.dim_u;
            dim_x_  = obj.dim_x;
            size_uT = dim_u_*T_;
            L_      = obj.L;
            U_      = obj.U0; % [1,dimuT]
            Uk      = reshape(U_,[dim_u_,T_])';
            % U = [u_1(t0),u_2(t0),...,u_dimu(t0); u_1(t1),u_2(t1),...,u_dimu(t1);...]
            % g1      = obj.index1; % 1-k番目のサンプルに対する入力の生成
            % g2      = obj.index2;
            SystemDynamics = @obj.Func;
            Clamp   = @obj.clamp;
            q_tf    = obj.phi_tf;
            q_t     = obj.phi_t;
            % gam_    = obj.gamma;
            invcov    = obj.Cov_;
            

            % obj.UU_sample;
            % XX   = obj.XX_sample;

            % state-cost S
            Sk = NaN(K_,1);

            % sampling
            ee = (NaN(K_,size_uT)); 
            % ee = [ [e(k1,u0(t0)); e(k2,u0(t0));...],[e(k1,u0(t1)); e(k2,u0(t1));...],...]
            for ie = 1 : T_
                ee(:,((ie-1)*dim_u_+1):(ie*dim_u_)) = (randn(K_,dim_u_)*L_); % L = chol(Cov)
            end
            % V = zeros(K_,size_uT);
            V = U_ + ee;
            
            % V(g1,:) = Uk + ee(g1,:);
            % V(g2,:) = ee(g2,:);

            
            
            instUU = zeros(K_,size_uT);
            instXX = zeros(K_,dim_x_*(T_+1));
            for k = 1 : K_ % parfor
                Vk = reshape(V(k,:),[dim_u_,T_])'; % Vk = [v_1(t0),v_2(t0); v_1(t1),v_2(t1); ...]
                Vk = Clamp(Vk);
                Ek = reshape(ee(k,:),[dim_u_,T_])';
                xx = NaN(T_+1,dim_x_);
                S = 0;
                xx(1,:) = x0'; % 行ベクトル
                for t = 1 : T_
                    vt = Vk(t,:);
                    ut = Uk(t,:);
                    et = Ek(t,:);
                    xx_t = xx(t,:);
                    S = S + q_t(xx_t) + lambda_*ut*invcov*et';

                    xx_t = SystemDynamics(xx_t,vt); % 行ベクトル入出力
                    xx(t+1,:) = xx_t;
                end
                S = S + q_tf(xx_t); % 行ベクトル入出力
                Sk(k) = S;
                instUU(k,:) = reshape(Vk',[1,dim_u_*T_]);
                instXX(k,:) = reshape(xx',[1,dim_x_*(T_+1)]);
            end
            % weight
            beta  = min(Sk);
            f     = exp((beta-Sk)/lambda_);
            eta   = sum(f,"all");
            ww    = (f/eta);
            % solution
            UU_opt = sum(ww.*instUU); % [1,uT]
            u0_opt = UU_opt(1:dim_u_)';


            obj.X_opt(1:dim_x_) = x0';
            for t2 = 1 : T_
                indu = ((t2-1)*dim_u_+1):(t2*dim_u_);
                indx = ((t2-1)*dim_x_+1):(t2*dim_x_);
                indx_= ((t2)*dim_x_+1):((t2+1)*dim_x_);
                xx0 = obj.X_opt(indx);
                obj.X_opt(indx_) = SystemDynamics(xx0,UU_opt(indu));

            end
            % update
            obj.U0    = UU_opt;%[UU_opt((dim_u_+1):end),UU_opt((end-dim_u_+1):end)];
            obj.U_opt = UU_opt;
            obj.UU_sample = instUU;
            obj.XX_sample = instXX;

            % obj.UU_sample = V;
            % u0_opt = CheckRun(x0);
        end

        function out = getData(obj,name)
            arguments
                obj
                name (1,1) string ...
                    {mustBeMember(name,["XData","UData","Cost","X_opt","U_opt","Cost_opt"])}
            end

            if strcmpi(name,"XData")
                out = obj.XX_sample;
            elseif strcmpi(name,"UData")
                out = obj.UU_sample;
            elseif strcmpi(name,"Cost")
                out = obj.Cost_sample;
                
            elseif strcmpi(name,"X_opt")
                out = reshape(obj.X_opt,obj.tra2x_opt);
            elseif strcmpi(name,"U_opt")
                out = reshape(obj.U_opt,obj.tra2u_opt);
            elseif strcmpi(name,"Cost_opt")
                out = obj.Cost_opt;
            end
            % if strcmpi(name,"XData")
            %     out = reshape(obj.XX_sample,obj.tra2data_x);
            % elseif strcmpi(name,"UData")
            %     out = reshape(obj.UU_sample,obj.tra2data_u);
            % elseif strcmpi(name,"Cost")
            %     out = obj.Cost_sample;
            % 
            % elseif strcmpi(name,"X_opt")
            %     out = reshape(obj.X_opt,obj.tra2x_opt);
            % elseif strcmpi(name,"U_opt")
            %     out = reshape(obj.U_opt,obj.tra2u_opt);
            % elseif strcmpi(name,"Cost_opt")
            %     out = obj.Cost_opt;
            % end
            % if strcmpi(name,"XData")
            %     out = reshape([obj.X_opt,obj.XX_sample],obj.tra2data_x);
            % elseif strcmpi(name,"UData")
            %     out = reshape([obj.U_opt,obj.UU_sample],obj.tra2data_u);
            % elseif strcmpi(name,"Cost")
            %     out = [obj.Cost_opt;obj.Cost_sample];
            % end
        end
    end
    methods (Static)
        function S = cost(x,Cost_tf,Cost_t,T)
            S = Cost_tf(x(end,:));
            for t = 1 : T
                S = S + Cost_t(x(t,:));
            end
        end
        function x = clamp_min(x,xmin,T,dimx)
            % x(T,u)
            % bx = repmat(xmin,[1,T]);
            for t = 1 : T
            for i = 1 : dimx
                if x(t,i)<xmin(i)
                    x(t,i)=xmin(i);
                end
            end
            end
        end
        function x = clamp_max(x,xmax,T,dimx)
            % x(T,u)
            for t = 1 : T
            for i = 1 : dimx
                if x(t,i)>xmax(i)
                    x(t,i)=xmax(i);
                end
            end
            end
        end
        function x = clamp_minmax(x,xmin,xmax,T,dimx)
            for t = 1 : T
            for i = 1 : dimx
                if x(t,i)>xmax(i)
                    x(t,i)=xmax(i);
                end
                if x(t,i)<xmin(i)
                    x(t,i)=xmin(i);
                end
            end
            end
        end
    end
end

