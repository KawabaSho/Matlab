clear, close all

addpath include\

% Initialization
Fig   = FigureManager("Sample");
Phys  = PhysicsManeger;

%  Setting Actors
dynamics = @Dynamics;
ActNum   = 100;
for i = 1 : ActNum
    x0       = [2*rand(9,1)-1; 1;0;0;0]; % initial condition
    Phys.addActor(RigidBody(x0, dynamics, [], []))
end

% Setting Plot
Fig.addList(Plot3D(Phys.getValue(1)));
hold on
for i = 2 : ActNum
    Fig.addList(Plot3D(Phys.getValue(i)));
end
Fig.setGraphic;

% Run
for i = 1 : 999
    Phys.Integration(rand(4,1), 0.01);
    for k = 1 : Fig.ListNum
        Fig.addData(k,Phys.getValue(k));
    end

    if ~mod(i,3)
        pause(0);
    end
end


function dxdt = Dynamics(x0,u)
    F_b = [0; 0; u(1) + u(2) + u(3) + u(4)];

    v = x0(Index.vv);
    w = x0(Index.ww);
    q = x0(Index.qq);

    qo = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);

    norm_q  = qo^2 + q1^2 + q2^2 + q3^2;

    Qu = [(qo^2 + q1^2 - q2^2 - q3^2), 2*(q1*q2 - qo*q3), 2*(q1*q3 + qo*q2);...
          2*(q1*q2 + qo*q3), (qo^2 - q1^2 + q2^2 - q3^2), 2*(q2*q3 - qo*q1);...
          2*(q1*q3 - qo*q2), 2*(q2*q3 + qo*q1), (qo^2 - q1^2 - q2^2 + q3^2)]/norm_q;
    
    Om = [0,    -w(1), -w(2), -w(3);
          w(1),     0,  w(3), -w(2);
          w(2), -w(3),     0,  w(1);
          w(3),  w(2), -w(1),    0];

    dxdt = [v;
            Qu*F_b./const.m;% - [0;0;const.g];
            -const.invII*(cross(w,const.II*w) - const.T*u + const.D*w);
            0.5*Om*q];
end




