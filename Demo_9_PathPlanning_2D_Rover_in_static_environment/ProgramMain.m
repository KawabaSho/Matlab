%% Collision avoidance with 2D Rover 
%{
    ローバーの障害物回避を行います．
    最適化ソルバは，SQP,AL-iLQRが使えます．flagで選択してください．
    N = 1000，dt = 0.01は実行に時間がかかります（1分間くらい．SQPでは30分くらい）．
    N =  100, dt = 0.1は高速ですが，軌道が離散的になりすぎて制約を満たしにくいです．

    state variables = [x,y,theta, dot_x,dot_y,dot_theta]
    input variables = [left_acceleration,　right_acceleration]
    2025/10/31 Kawarabayashi
%}

%%
clear, close all
%%
flag = 1; % solver setting 1: al-ilqr 0: sqp
%% pridiction parms
NN = 1000;
dt = 0.01;
tf = 10.0;
% alilqr NN = size(dt:dt:tf), size(uu) -> NN

%% actor params
x0 = [0;0;pi/180*20; 0;0;0];        % initial states
x0_1 = [0.6;-100.1];  % obstacle 1　ここの第3成分を0.1にすると，実行不可能解が得られます．
x0_2 = [0.45;0.3]; % obstacle 2
r_0 = 0.075; % radius
r_1 = 0.1;   % obstacle 1 radius 
r_2 = r_1;   % obstacle 2 radius

%% rober dynamics
% x = [x,y,theta, dot_x,dot_y,dot_theta], u = [left_acceleration,　right_acceleration]
function xpp = Dynamics(x,u,dt)
    L = 0.1;
    uu = [cos(x(3)),sin(x(3))]*(u(1)+u(2));
    dxdt = [x(4);
            x(5);
            x(6);
            -10*x(4) + uu(1);
            -10*x(5) + uu(2);
            -0.001*x(6) + L*(-u(1)+u(2))];
    xpp = x + dt*dxdt;
end
function [A,B] = linearization(x,u,dt)
    A = eye(6) +...
        [0,0,0,1,0,0;
         0,0,0,0,1,0;
         0,0,0,0,0,1;
         0,0,-(u(1)+u(2))*sin(x(3)),-10,0,0;
         0,0,(u(1)+u(2))*cos(x(3)),0,-10,0;
         0,0,0,0,0,-0.001;]*dt;
    L = 0.1;
    % uu = [(u(1)+u(2))*cos(x(3)),(u(1)+u(2))*sin(x(3))];
    B = [0,0;
         0,0;
         0,0;
         cos(x(3)),cos(x(3));
         sin(x(3)),sin(x(3));
         -L,L]*dt;

end
function c = con_position(xx,~,NN,r0,r1,x1,x2)
    % con_position(xx,uu,NN)
    c = zeros(2,NN);
    xx_veh = xx(1:2,1:end-1);
    r2 = (r0+r1);
    x = [x1,x2];
    for i = 1 : 2
        dx = (x(1:2,i)- xx_veh)';
        for k = 1 : NN
            dxi = dx(k,:);
            c(i,k) = 1*(r2-sqrt(dxi*dxi'));
        end
    end
end
function c = con_position_grad(xx,~,NN,x1,x2)
    % con_position_grad(xx,uu,NN)
    xx_veh = xx(1:2,1:end-1);
    c = zeros(2,6,NN);
    x = [x1,x2];
    for i = 1 : 2
        dx = (x(1:2,i)- xx_veh)';
        for k = 1 : NN
            dxi = dx(k,:);
            r = sqrt(dxi*dxi');
            c(i,1:2,k) = 1*dxi./r;
        end
    end
end
function c = con_input(~,uu)
    % con_input(xx,uu,umax,umin)
    umax = [10;10];
    umin = [0;0];
    c = [uu-umax;umin-uu];
end
function cu = con_input_grad(~,~,NN)
    % con_input_grad(xx,uu,NN)
    cu_i = [  1,0;
              0,1;
              -1,0;
              0,-1;];
    cu = repmat(cu_i,1,1,NN);
end

%% optimization problem
function J = criterion(uu,x0,xg,ug,NN,dt,Q,q,R,r)
    xx = prediction(uu,x0,NN,dt);
    dxx = xx - xg;
    duu = uu - ug;
    cost1 = dxx(:,end)'*Q*dxx(:,end)/dt; % Qには予めdtが掛けられているため終端コストをdtで割る必要がある
    cost2 = q'*dxx(:,end);
    for i = 1 : NN
        dx = dxx(:,i);
        du = duu(:,i);
        cost1 = cost1 + dx'*Q*dx + du'*R*du;
        cost2 = cost2 + q'*dx + r'*du;
    end
    J   = cost1*0.5 + cost2;
end
function [con, ceq] = con(uu,x0,NN,dt, r0,r1,x1,x2)
    xx = prediction(uu,x0,NN,dt);
    con = [con_input([],uu);con_position(xx,[],NN,r0,r1,x1,x2)];
    ceq = [];
end
function xxo = prediction(uu,x0,NN,dt)
    xx = zeros(6,NN+1);
    xx(:,1) = x0;
    for i = 1 : NN
        xx(:,i+1) = Dynamics(xx(:,i),uu(:,i),dt);
    end
    xxo = xx;
end
xg = [1;0.4;0; 0;0;0];   % goal
ug = [0;0];              % goal input
Q  = [5,0,0,0,0,0;
      0,5,0,0,0,0;
      0,0,1e-8,0,0,0;
      0,0,0,1,0,0;
      0,0,0,0,1,0;
      0,0,0,0,0,1e-8;];  % *dt は含まれる
q  = [0;0;0;0;0;0];
R  = 0.1*[1,0;0,1];      % *dt は含まれる
r  = [0;0];
%% initial guess
u0 = zeros(2,NN)+[1;1];

if flag
    %% AL-iLQR
    S = Q/dt; % 終端の重み　fiminconとal-ilqrで同じ評価関数にするための補正
    tic
    AL = solver_PreSearch_ALiLQR_method_7(Q,[],R,[],[],S,[],dt,tf,...
                                            @(x,u)linearization(x,u,dt),...
                                            @(x,u)Dynamics(x,u,dt),1e-18,100);
    AL.updateTerminalValues(xg,ug)
    AL.set_InitialGuess_RepMat(u0(:,1))
    uu0 = AL.UU;
    xx0 = zeros(AL.size_state,AL.NN+1);
    for ip = 1 : AL.NN
    xx0(:,ip+1) = AL.hf_DiscreteDynamics(xx0(:,ip),uu0(:,ip));
    end
    AL.set_InitialState(xx0)
    AL.setUserInterface([]);
    AL.setUserInterface_end([])
    AL.setconstraint(...
    @(xx,uu)con_input(xx,uu),...    % con_u
    [],...                                                     % ceq_u
    @(xx,uu)con_input_grad(xx,uu,AL.NN),...    % Gradient con_u
    [],...                                                     % Gradient ceq_u
    @(xx,uu)con_position(xx,uu,AL.NN,r_0,r_1,x0_1,x0_2),...    % con_x
    [],...                                                     % ceq_x
    @(xx,uu)con_position_grad(xx,uu,AL.NN,x0_1,x0_2),...       % Gradient con_x
    [],...                                                     % Gradient ceq_x
    [],...                                                     % con_x(tf)
    [],...%@(x)ceq_state_tf(x),...                             % ceq_x(tf)
    [],...                                                     % Gradient con_x(tf)
    [],"Tolerance",1e-18,"Constraint_MaxIteration",20)%@ceq_state_grad_tf)                                    % Gradient ceq_x(tf)
    AL.setconstraint2(...
    [],... % con_u
    [],...                                                     % ceq_u
    [],...                    % Gradient con_u
    [],...                                                     % Gradient ceq_u
    [],...                      % con_x
    [],...                         % ceq_x
    [],...                 % Gradient con_x
    [],...                    % Gradient ceq_x
    [],...                                                     % con_x(tf)
    [],...                                                     % ceq_x(tf)
    [],...                                                     % Gradient con_x(tf)
    [])                                                        % Gradient ceq_x(tf)
    AL.Run(x0);                    % AL-iLQR 実行
    toc
    xx_opt = AL.XX;
    uu_opt = AL.UU;
    cc = con(uu_opt,x0,NN,dt, r_0,r_1,x0_1,x0_2);
    fprintf("AL-iLQRの評価関数は　%10f\n",criterion(uu_opt,x0,xg,ug,NN,dt,Q,q,R,r))
    fprintf("制約の最大違反は　%10e\n",max(cc.*(cc>0),[],"all"))
else
    %% fmincon
    % opts = optimoptions('fmincon','Algorithm','interior-point');
    opts = optimoptions('fmincon','Algorithm','sqp');%,'PlotFcn','optimplotconstrviolation');
    tic
    uu_opt = fmincon(@(u)criterion(u,x0,xg,ug,NN,dt,Q,q,R,r),u0,[],[],[],[],[],[],...
        @(u)con(u,x0,NN,dt, r_0,r_1,x0_1,x0_2),opts);
    toc
    xx_opt = prediction(uu_opt,x0,NN,dt);
    fprintf("Fminconの評価関数は　%10f\n",criterion(uu_opt,x0,xg,ug,NN,dt,Q,q,R,r))
    cc = con(uu_opt,x0,NN,dt, r_0,r_1,x0_1,x0_2);
    fprintf("制約の最大違反は　%10e\n",max(cc.*(cc>0),[],"all"))
    % 4.52693e-11
end
%% rendering
cirD = 0 : pi/180 : 2*pi;
nxr = cos(cirD);
nyr = sin(cirD);
Cmap = Colormap_Turbo;

figure(1)
plot(x0(1)+r_0*nxr,x0(2)+r_0*nyr,'k')
hold on
plot(x0_1(1)+r_1*nxr,x0_1(2)+r_1*nyr,'k')
plot(x0_2(1)+r_2*nxr,x0_2(2)+r_2*nyr,'k')
plot(xg(1),xg(2),'ob')
cd = scatter(xx_opt(1,:),xx_opt(2,:),6,Cmap.CData(NN+1),"filled");
hold off
axis equal
ylim([-0.2 0.6])
xlim([-0.2 1.3])

figure(2);
plot((dt:dt:tf)-dt,uu_opt(1,:),'DisplayName','Left input')
hold on
plot((dt:dt:tf)-dt,uu_opt(2,:),'DisplayName','Right input')
hold off
legend
xlabel("Time")

figure(3);
plot(0:dt:tf,xx_opt(1,:),'DisplayName','X')
% ylabel("[-]","Color",'k')
hold on
plot(0:dt:tf,xx_opt(2,:),'DisplayName','Y')
% ylabel("","Color",'k')
hold off
legend
xlabel("Time")



