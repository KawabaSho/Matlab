clear, close all
% STLファイルを読み込みます

%  Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fig   = FigureManager("Name","Sample");
Phys  = PhysicsManeger;

%  Setting Actors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynamics = @Dynamics;
x0 = [1;0; 0; 0; 1/sqrt(2); 0];
Phys.addActor(RigidBody(x0, dynamics, [], []))

%  Setting Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TR0 = stlread("sample.stl");
model.Faces    = TR0.ConnectivityList;
xx = Phys.getValue(Phys.ActorNum);
model.Vertices = TR0.Points + xx(1:3)';
Fig.addList(...
    Patch([0;0;0], TR0.Points, model.Faces,...
        model,'FaceColor',[0.7 0.7 1.0],'EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.5));

grid on
axis equal
lightangle(-45,70)
view(45,45)
xlim([-2 2])
ylim([-2 2])
zlim([-2 2])

%  Run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : 100000
    Phys.Integration([0;0], 0.01);
    for k = 1 : Fig.ListNum
        x_curr = Phys.getValue(k);
        Fig.PlotList{k}.Refresh(RodoriguesRotation([0;0;1],i*0.01),x_curr(1:3));
    end

    pause(0);
end


function dxdt = Dynamics(x0,u)
    xx = x0(1);
    xy = x0(2);
    vx = x0(4);
    vy = x0(5);
    L = xx^2 + xy^2;
    dxdt = [vx;
            vy;
            0;
            -xx/L + u(1);
            -xy/L + u(2);
            0];
end

function Rotation = RodoriguesRotation(x, rad)
    cs = cos(rad); sn = sin(rad);
    Rotation = ...
    [cs+x(1)*x(1)*(1-cs), x(1)*x(2)*(1-cs)-x(3)*sn, x(1)*x(3)*(1-cs)+x(2)*sn;
     x(1)*x(2)*(1-cs)+x(3)*sn, cs+x(2)*x(2)*(1-cs), x(2)*x(3)*(1-cs)-x(1)*sn;
     x(1)*x(3)*(1-cs)-x(2)*sn, x(2)*x(3)*(1-cs)+x(1)*sn, cs+x(3)*x(3)*(1-cs)];
    Rotation = Rotation/det33(Rotation);
end
function det = det33(RR)
    det =  (RR(1,1)*RR(2,2)*RR(3, 3) + RR(1, 2)*RR(2, 3)*RR(3, 1) + RR(1, 3)*RR(2, 1)*RR(3, 2))...
           -(RR(1, 3)*RR(2, 2)*RR(3, 1) + RR(1, 2)*RR(2, 1)*RR(3, 3) + RR(1, 1)*RR(2, 3)*RR(3, 2));
end

