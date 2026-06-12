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
tt_ = 0;

% 状態量と制御入力の次元
x_dim = 3;
u_dim = 2;

% 障害物配置と衝突半径の設定
random = 0;
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

% 目標位置
x_goal = [9;9];

%% Log
log_tt = NaN(1,NN_+1);
log_xx = NaN(x_dim,NN_+1);
log_uu = NaN(u_dim,NN_);

%% Initialize
xx_ = [-9;-9;pi/4*0];  % [x,y,theta]
% xx_ = [0;0;0];

%% save
log_tt(1,1) = tt_;
log_xx(:,1) = xx_;

%% Controller initialize
controller = @contorller_P;

%% Render initialize
fig1 = figure;
theta = linspace(0, 2*pi, 100);

% goal
graph_g = fill(...
        0.1*cos(theta)+x_goal(1,1),...
        0.1*sin(theta)+x_goal(2,1),...
        [0.75,0,0], 'EdgeColor', [0.75,0,0], 'LineWidth', 0.1);
hold on

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
try
for i_render = 1 : NN_/plotperiod
    for i_time = 1 : plotperiod

        xx = xx_;
        uu = controller(xx,x_goal);
        
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
    pause(0);
end
catch e
    fprintf(2,e.message)
    return;
end


%% functions
function dxdt = Func(x,u,~)
    % x = [x,y,theta], u = [v, omega]
    v = u(1);
    w = u(2);
    theta = x(3);
    c = cos(theta);
    s = sin(theta);
    dxdt = [v*c; v*s; w];
end
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
    mes = "\t障害物 "+num2str(i)+" さんに衝突しました．制御器に問題があるようです．\n\t";
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
    
    error(mes + "通算衝突回数は %d です．\n",count);
end
%% user functions
function u0_opt = contorller_P(xx,goal)
    rr = goal - xx(1:2);
    r = sqrt(rr'*rr);
    nn = rr/r;
    if r < 0.4
        v = r;
    else
        v = 0.4;
    end
    u0_opt = [v; -0.5*(nn(1)*sin(xx(3)) - nn(2)*cos(xx(3)))];
end
















