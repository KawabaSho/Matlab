clear,close all
%{
    差動二輪モデルの障害物回避
    モデル：
    x = [position_x,position_y,angle], u = [velocity, angular_rate]
    dx/dt = [velocity*cos(angle);
             velocity*sin(angle);
             angular_rate]

    03/21/2025 Kawarabayashi
%}
%% ERROR FLAG
%#ok<*SAGROW>

%% Configuration
NN_ = 50000;     % iteration
dt_ = 0.01;      % sec

plotperiod = 10; %
% t0  = 0;
tt_ = 0;

x_dim = 3;
u_dim = 2;

%% Log
log_tt = NaN(1,NN_+1);
log_xx = NaN(x_dim,NN_+1);
log_uu = NaN(u_dim,NN_);

%% Initialize
xx_ = [-9;-9;pi/4];  % [x,y,theta]
% xx_ = [0;0;0];

%% save
log_tt(1,1) = tt_;
log_xx(:,1) = xx_;

%% Controller initialize
% MPPI
T      = 10; % Number of timesteps
K      = 50; % Number of samples
UU     = zeros(u_dim,T);
Cov    = diag([0.24,0.24]);
lambda = 2;
x_goal = [9;9];

random = 70;
N_obs  = 8 + random;
xx_obs = ...
    [[0;0],[-1.4;0],[0;-1.4],[-1;-2],[4;2],...
    [-5.72981858196456;-1.25181947797960],...
    [ 6.65176840302507; 4.67531727295287],...
    [ 7.35187882228645; 2.49185118650539],...
    (2*rand(2,random)-1)*8];
r_obs  = [1.1,0.5,0.5,0.2,0.8,0.6,...
          0.7,0.6,0.5,...
          rand(1,random)];
r_veh  = 0.1;
Cost = @(x)MPPI_Cost(x,x_goal',N_obs,xx_obs',r_obs,r_veh); % Cost(x(t))
mppi = class_mppi(Cost,Cost,T,K,UU,Cov,lambda,...
    @(x,u)Func4mppi(x,u,dt_*100),length(xx_),...
    "umax",[0.5,pi/2/2/2],"umin",[0.5,-pi/2/2/2]);

%% Render initialize
fig1 = figure;
theta = linspace(0, 2*pi, 100);

% goal
graph_g = fill(...
        0.1*cos(theta)+x_goal(1,1),...
        0.1*sin(theta)+x_goal(2,1),...
        [0.75,0,0], 'EdgeColor', [0.75,0,0], 'LineWidth', 0.1);
hold on

% sample
for ip = 1 : K
    graph_mppi{ip} = ...
        plot(NaN(1,T+1),NaN(1,T+1),'Color',[[0.496,0.702,1]*0.8,0.5]); 
end
graph_mppi{K+1} = ...
        plot(NaN(1,T+1),NaN(1,T+1),'Color',[0,1,0]);

% obstacles
for i_ro = 1 : N_obs
    gra{i_ro} = fill(...
        r_obs(i_ro)*cos(theta)+xx_obs(1,i_ro),...
        r_obs(i_ro)*sin(theta)+xx_obs(2,i_ro),...
        [0.75,0.75,0.75], 'EdgeColor', 'k', 'LineWidth', 0.1);
end

% vehicle
graph1 = plot(log_xx(1,:),log_xx(2,:),'k');
graph2 = quiver(log_xx(1,1),log_xx(2,1),cos(log_xx(3,1)),sin(log_xx(3,1)),0.5,"filled",'r','LineWidth',0.5);
graph3 = quiver(log_xx(1,1),log_xx(2,1),-sin(log_xx(3,1)),cos(log_xx(3,1)),0.5,"filled",'b','LineWidth',0.5);
theta = linspace(0, 2*pi, 20);
rc = r_veh*cos(theta);
rs = r_veh*sin(theta);
graph4 = fill(rc+log_xx(1,1),rs+log_xx(2,1),...
    [0.496,0.702,1], 'EdgeColor', [0.496,0.702,1], 'LineWidth', 0.01);



hold off
axis equal
xlim([-10 10])
ylim([-10 10])

%% Render loop
i = 1;
% try
for i_render = 1 : NN_/plotperiod
    for i_time = 1 : plotperiod

        xx = xx_;
        uu = mppi.run(xx);
        % uu = [0.1 + i*0.00001; 0.1];
        
        k1 = dt_*Func(xx,uu);
        k2 = dt_*Func(xx + 0.5*k1,uu);
        k3 = dt_*Func(xx + 0.5*k2,uu);
        k4 = dt_*Func(xx + k3,uu);
        xx_ = xx + (k1 + k4 + 2*(k2 + k3))/6;
        tt_ = tt_ + dt_;

        log_tt(1,i+1) = tt_;
        log_xx(:,i+1) = xx_;
        log_uu(:,i) = uu;
        i = i + 1;
        CollisionCheck(xx,N_obs,xx_obs,r_obs,r_veh)
    end

    ix = xx_(1);
    iy = xx_(2);
    graph1.XData(i_render) = ix;
    graph1.YData(i_render) = iy;
    c = cos(xx_(3));
    s = sin(xx_(3));
    graph2.XData    = ix;
    graph2.YData    = iy;
    graph2.UData    = c;
    graph2.VData    = s;
    graph3.XData    = ix;
    graph3.YData    = iy;
    graph3.UData    = -s;
    graph3.VData    = c;
    graph4.XData    = rc+ix;
    graph4.YData    = rs+iy;

    % sample 
    XX = mppi.getData("XData");
    for ip = 1 : K
        xy = reshape(XX(ip,:),[x_dim,T+1]);
        graph_mppi{ip}.XData = xy(1,:);
        graph_mppi{ip}.YData = xy(2,:);
    end
    
    xy2 = mppi.getData("X_opt");
    graph_mppi{K+1}.XData = xy2(1,:);
    graph_mppi{K+1}.YData = xy2(2,:);
    
    pause(0);
end
% catch e
%     fprintf(2,e.message + "\n")
%     return;
% end

function dxdt = Func(x,u,~)
    % x = [x,y,theta], u = [v, omega]
    v = u(1);
    w = u(2);
    theta = x(3);
    c = cos(theta);
    s = sin(theta);
    dxdt = [v*c; v*s; w];
end
function x_ = Func4mppi(x,u,dt)
    % x = [x,y,theta], u = [v, omega]
    v = u(1);
    w = u(2);
    theta = x(3);
    c = cos(theta);
    s = sin(theta);
    x_ = x + [v*c, v*s, w]*dt;
end
%%
% function [A,B] = ContinuousLinearFunc(x,u,t)
%     v = u(1);
%     s = sin(x(3));
%     c = cos(x(3));
%     A = [0,0,-v*s;0,0,v*c;0,0,0;];
%     B = [c,0;s,0;0,1];
% end
% 
% function [A,B] = DiscreteLinearFunc(x,u,t,dt)
%     v = u(1);
%     s = sin(x(3))*dt;
%     c = cos(x(3))*dt;
%     A = eye(3) + [0,0,-v*s;0,0,v*c;0,0,0;];
%     B = [c,0;s,0;0,dt];
% end
function CollisionCheck(xx,N_obs,xx_obs,r_obs,r_veh)
    x  = xx(1:2);
    for i = 1 : N_obs
        dxo = x - xx_obs(:,i);
        ro  = r_obs(i) + r_veh;
        if (ro*ro - dxo'*dxo) > 0
            ErrorProcess(i);
        end
    end
end
function ErrorProcess(i)
    mes = "障害物 "+num2str(i)+" さんに衝突しました．制御器に問題があるようです．";
    dtime = datetime("now","Format","MM/dd/uuuu HH:mm:ss");
    
    % read
    fileID = fopen('Log_error.txt','r');
    if fileID == -1
        fileID = fopen('Log_error.txt','w');
        fprintf(fileID,"Error History:\n");
        count = 1;
    else
        while ~feof(fileID)
            line = fgetl(fileID);
        end
        count = str2double(line(1:4)) + 1;
    end
    fclose(fileID);
    
    % write
    fileID = fopen('Log_error.txt','a');
    fprintf(fileID,num2str(count,'%04d') + " " + string(dtime)  + " " + "(0x41C) Control failure\n");
    fclose(fileID);
    
    error("\n\t" + mes + "\n\t通算衝突回数は %d です．",count);
end

%% MPPI 
function c = MPPI_Cost(xx,x_goal,N_obs,xx_obs,r_obs,r_veh)
    % x = [x(t),y(t),theta(t)], x_goal = [x,y]
    x  = xx(1:2);
    dx = x - x_goal;
    c  = sqrt(dx*dx');
    
    % obstacles
    for i = 1 : N_obs
        dxo = x - xx_obs(i,:);
        ro  = r_obs(i) + r_veh;
        if (ro*ro - dxo*dxo') > 0
            c = c + 1000;
        end
    end
end




















