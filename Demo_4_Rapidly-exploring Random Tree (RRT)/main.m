clear,close all



q_init = [-1;-1;-1];
Area = [-1,1;-1,1;-1,1];

rrt = class_RRT(Area,0.05,1000);

% obstacle
dX        = 0.75;
Pos_obj_c = -[dX;dX;dX]*0.5;
NN = [2,2,2];
ObsNum = NN(1)*NN(2)*NN(3);
Obstacles = class_Check_Collision();

for k = 1 : NN(3)
    k_th = NN(1)*NN(2)*(k-1);
    for j = 1 : NN(2)
        j_th = NN(1)*(j-1);
        for i = 1 : NN(1)
            id =  i + j_th + k_th;
            model = generateCube(0.1);
            model.Vertices = model.Vertices + Pos_obj_c' ...
                + [dX*(i-1), dX*(j-1), dX*(k-1)];
            rrt.addObstacle(model);
        end
    end
end


tic
rrt.Run(q_init,0);
toc
% Graph

[Ver, Num_v] = rrt.Graph.getVertex;
[Edge, Num_e] = rrt.Graph.getEdge;



figure()
hold on
for i = 1 : Num_e
    q_1 = [Ver(1,Edge(1,i)), Ver(1,Edge(2,i))];
    q_2 = [Ver(2,Edge(1,i)), Ver(2,Edge(2,i))];
    q_3 = [Ver(3,Edge(1,i)), Ver(3,Edge(2,i))];
    plot3(q_1,q_2,q_3,'Color','k')
end

plot3(Ver(1,:),Ver(2,:),Ver(3,:),'LineStyle','none','Marker','.','MarkerSize',5,'Color','b');
for R_th = 1 : ObsNum
    patch(rrt.Collision.Actors{R_th},'FaceColor', [0.8 0.8 1.0], ...
        'EdgeColor', 'none', ...
        'FaceLighting', 'gouraud', ...
        'AmbientStrength', 0.15)
end
hold off

axis equal
grid on
xlim(Area(1,:))
ylim(Area(2,:))
zlim(Area(3,:))
view([20,20])
lightangle(-45,70)
Record_RotatingGraph(gcf,"Output.avi");




















