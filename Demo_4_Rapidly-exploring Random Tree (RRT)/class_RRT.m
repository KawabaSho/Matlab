classdef class_RRT < handle

    properties
        Graph  class_Graph% Graph
        Collision class_Check_Collision
        Area % [q_max, q_min];
        L
        MaxIteration
        Obstacles % {}
        counter

        Tree class_Octree
    end

    methods
        function obj = class_RRT(Area,StepSize,MaxIteration)
            % Area:[min,max], StepSize:L
            obj.Graph = class_Graph;
            obj.MaxIteration = MaxIteration;
            obj.Area = Area;
            obj.L    = StepSize;
            obj.counter = 0;
            obj.Collision = class_Check_Collision;
            
            AreaL         = Area(:,2) - Area(:,1);
            obj.Tree      = class_Octree(3,Area,AreaL/10);

        end
        function addObstacle(obj, model)
            % model.Vertices model.Faces
            q_max = model.Vertices(1,:);
            q_min = q_max;
            for i = 1 : size(model.Vertices,1)
                for k = 1 : 3
                    if q_min(k) > model.Vertices(i,k)
                        q_min(k) = model.Vertices(i,k);
                    end
                end
                for k = 1 : 3
                    if q_max(k) < model.Vertices(i,k)
                        q_max(k) = model.Vertices(i,k);
                    end
                end 
            end
            obj.Tree.AddObject(q_min',q_max');

            obj.Collision.AddObstacle(model)
        end
        function obj = Run(obj,q_init,goal)
            
            dim      = size(q_init,1);
            MaxItre  = obj.MaxIteration;
            box      = zeros(dim,MaxItre+1);
            box(:,1) = q_init;

            Ver   = box;  % Position 
            Edg   = zeros(2,MaxItre);
            A_min = obj.Area(:,1);
            A     = (obj.Area(:,2) - A_min);
            d     = obj.L;
            for n = 1 : MaxItre
                q_rand = A.*rand(dim,1) + A_min;
                ne     = q_rand - Ver(:,1:n);
                ee     = obj.norm2(ne);
                [distance,Id] = min(ee);
                q_near = Ver(:,Id);
                rr     =  q_rand - q_near;
                q_new  =  q_near + rr/sqrt(rr'*rr)*d;

                flag_col = 0;
                id = obj.Tree.getIndex(q_new,2); % all obstacles level is 2
                for i = 1 : obj.Collision.Num
                    collision_tree = obj.Tree.Object2Tree{i};
                    if id == collision_tree(2)
                        flag_col = flag_col + obj.Collision.Check(i,q_new);
                    end
                end
                if flag_col
                    Ver(:,n+1) = NaN(size(q_init));
                    Edg(:,n)   = [1,1];
                else
                    Ver(:,n+1) = q_new;
                    Edg(:,n)   = [Id,n+1];
                end
                
            end
            obj.Graph.Vertex = Ver;  % Position 
            obj.Graph.Edge   = Edg;
        end
    end
    methods (Static)
        function ne = norm2(ee)
            % [q1,q2,q3,...]
            S = size(ee,2);
            ne = zeros(1,S);
            for i = 1 : S
                e = ee(:,i);
                ne(i) = e'*e;
            end
        end
    end

end