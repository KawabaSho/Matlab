%{
    Dec.9th,2023
    Koboensyu V  Swing-by EDVEGA(Electric Delta-V Earth Gravity Assist)

    Reference : 「高校数学・物理で学ぶスイングバイ軌道の作り方
                                    〜スイングバイの軌道計算をしてみよう〜」
    https://edu.jaxa.jp/contents/english/homework/swingby2.html

    Author : Kawarabayashi
%}

clear, close all
%  Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fig   = FigureManager("Name","Sample","Position",[1100,500, 510, 400]);
Phys  = PhysicsManeger;

%  Params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% _e : earth
% _s : sun
% _c : chaser, satellite
AU      = 0.149597870697400004e+9;
sec2day = 24*60*60;
r_e     = 6378.1363/AU;

mu_e   = 398600.4418/AU^3 * sec2day^2;     % [AU^3/day^2]
mu_s   = 0.295912208285591149e-3;          % [AU^3/day^2]

T_e    = 365.2422;                         % [day] Period
w_e    = 1.9910637973e-7 * sec2day;        % [rad/day]
v_e    = w_e*1;                            % 地球公転速度

%  initial condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% earth
x_e  = [1, 0, 0]';                 % 地球の初期位置
x0_e = [x_e;...
        [0, sqrt(mu_s), 0]'];
% chaser,satellite
vinf = 5/AU*sec2day;               % 地球公転速度に対する相対速度V∞, AU/day
vv_s = [-sqrt(4*v_e*v_e - vinf*vinf)*vinf*0.5/v_e, (2*v_e*v_e - vinf*vinf)*0.5/v_e,0]';
x0_c = [1,0,0,vv_s']';

% index
p_xc = 1:3;
p_vc = 4:6;
p_xe = 7:9;
p_ve = 10:12;

%  Setting Actors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = [x0_c; x0_e];
dynamics = @(x,u)Dynamics_2Body(x,u,p_xc,p_vc,p_xe,p_ve,mu_e,mu_s);
Phys.addActor(Physics(x0, ...
    @(x0,u0,dt)Runge_Kutta_Nystrom_solver(x0,u0,dt,dynamics,[p_xc,p_xe],[p_vc,p_ve])))


%  Setting Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fig.addList(Plot3D(x0_c(1),x0_c(2),x0_c(3),'Color','g'));
hold on
Fig.addList(Plot3D(x0_e(1),x0_e(2),x0_e(3),'Color','b'));
Fig.addList(Plot3D(x0_c(1),x0_c(2),x0_c(3),'o','Color','g'));
Fig.addList(Plot3D(x0_e(1),x0_e(2),x0_e(3),'o','Color','b'));
Fig.addList(Plot3D(0,0,0,'o','Color','r'));
Fig.setLegend(["Probe","Earth"])
Fig.setLabel("X [AU]", "Y [AU]")

grid on
axis equal
lightangle(-45,70)
view(0,90)
xlim([-2 2])
ylim([-2 2])
zlim([-2 2])

%  Run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt  = 0.01; %[day]
tt_ = 0;
flag_swingby = true;
for i = 1 : 20000000
    % initial constrains %%%%%%%%%
    if tt_ > 0.1
        dynamics = @(x,u)Dynamics_3Body(x,u,p_xc,p_vc,p_xe,p_ve,mu_e,mu_s);
        Phys.PhysicsHandles{1}.Integrator = ...
        @(x0,u0,dt)Runge_Kutta_Nystrom_solver(x0,u0,dt,dynamics,[p_xc,p_xe],[p_vc,p_ve]);
    end

    % sekibun
    Phys.Integration(0, dt);

    % Swing-by %%%%%%%%%%%%%%%%%%%%
    if (tt_ > 0.20*T_e)&&(flag_swingby)
        x_curr = Phys.getValue(1);
        x_c = x_curr(p_xc);
        v_c = x_curr(p_vc);

        % Swing-by position after 250 days %%%%%%%%%%%%
        % tf = 100*250;
        % ratio = 100;
        % optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dV = [0;0;0];
        % f = @(dV)performanceIndex(dV,x_curr,v_c,dt,tf,Phys,p_xc,p_vc,p_xe,p_ve,r_e,ratio);
        % options = optimoptions('fmincon','Algorithm','sqp',...
        %                        'MaxFunctionEvaluations',90000,...
        %                        'MaxIterations',20000,...
        %                        'Display','none',...
        %                        'UseParallel',true);
        % dV = fmincon(f,dV,[],[],[0,0,1],0,[],[],[],options);

        % Swing-by sample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ratio = 100;tf = 100*250;
        % dV = [1.829646517778439e-04, 9.948943339143992e-04,0]'; % normal swing-by after 250 days r_e*100
        % ratio = 100;tf = 100*50;
        % dV = [0.004806343762960, 0.004827698782542,0]'; % short1 swing-by after 50 days r_e*100
        % ratio = 100;tf = 100*100;
        dV = [0.002303762634479, 0.003785758390392,0]'; % short2 swing-by after 100 days r_e*100
        % ratio = 0.01;tf = 100*200;
        % dV = [5.926304214337331e-04, 0.001945349294923,0]'; % crazy swing-by after 200 days r_e*0.01

        x_curr(p_vc) = v_c + dV;
        Phys.PhysicsHandles{1}.values = x_curr;
        
        Fig.addList(Plot3D(x_c(1),x_c(2),x_c(3),'*','Color','k'));
        Fig.setLegend(["Probe","Earth"])
        flag_swingby = false;
    end

    % draw stock
    if ~mod(i,100)
        x_curr = Phys.getValue(1);
    
        % satellite
        x_c = x_curr(p_xc);
        % v_c = x_curr(p_vc);
    
        % earth
        x_e = x_curr(p_xe);
        % v_e = x_curr(p_ve);
    
        Fig.addData(1,x_c);
        Fig.addData(2,x_e);
        Fig.setData(3,x_c);
        Fig.setData(4,x_e);
    end
    % draw
    if ~mod(i,200)
        pause(0);
    end
    tt_ = i*dt;
end


function dxdt = Dynamics_2Body(x0,u,p_xc,p_vc,p_xe,p_ve,mu_e,mu_s)
    % satellite
    x_c = x0(p_xc);
    v_c = x0(p_vc);
    % earth
    x_e = x0(p_xe);
    v_e = x0(p_ve);
    
    mur_se3 = mu_s/(x_e(1)^2+x_e(2)^2+x_e(3)^2)^1.5;
    mur_sc3 = mu_s/(x_c(1)^2+x_c(2)^2+x_c(3)^2)^1.5;

    dxdt = [v_c;
            - mur_sc3*x_c;
            v_e;
            - mur_se3*x_e];
end
function dxdt = Dynamics_3Body(x0,u,p_xc,p_vc,p_xe,p_ve,mu_e,mu_s)
    % satellite
    x_c = x0(p_xc);
    v_c = x0(p_vc);
    % earth
    x_e = x0(p_xe);
    v_e = x0(p_ve);
    
    mur_se3 = mu_s/(x_e(1)^2+x_e(2)^2+x_e(3)^2)^1.5;
    mur_sc3 = mu_s/(x_c(1)^2+x_c(2)^2+x_c(3)^2)^1.5;
    e = x_c - x_e;
    mur_ec3 = mu_e/(e(1)^2+e(2)^2+e(3)^2)^1.5;

    dxdt = [v_c;
            - mur_sc3*x_c - mur_ec3*e;
            v_e;
            - mur_se3*x_e];
end
function J = performanceIndex(dV,x_curr,v_c,dt,tf,Phys,p_xc,p_vc,p_xe,p_ve,r_e,ratio)
        % Swing-by position after tf steps
        x_curr(p_vc) = v_c + dV;
        Phys.PhysicsHandles{1}.values = x_curr;
        for i_m = 1 : tf
            Phys.Integration(0, dt);
        end
        xhat = Phys.getValue(1);
        xhat_e = xhat(p_xe);
        xhat_c = xhat(p_xc);
        J = norm(xhat_c - xhat_e*(1 + r_e*ratio));
end







