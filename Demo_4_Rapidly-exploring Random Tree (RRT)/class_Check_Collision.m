classdef class_Check_Collision < handle
    properties
        Actors
        Num
    end

    methods
        function obj = class_Check_Collision()
            obj.Num = 0;
            obj.Actors = {};
        end
        function AddObstacle(obj,model)
            % model.Vertices model.Faces
            obj.Num = obj.Num + 1;
            obj.Actors{obj.Num} = model;
        end
        function flag = Check(obj,id,vec_q)
            flag = 0;
            Obs = obj.Actors{id};
            face_num = size(Obs.Faces,1);
            for i = 1 : face_num
                p0 = Obs.Vertices(Obs.Faces(i,1),:)';
                p1 = Obs.Vertices(Obs.Faces(i,2),:)';
                p2 = Obs.Vertices(Obs.Faces(i,3),:)';

                nn = cross(p1-p0,p2-p0);
                if dot(nn,vec_q-p0) < 0
                    flag = flag + 1;
                end
            end
            if flag == face_num
                % pause;
                flag = 1;
            else
                flag = 0;
            end
        end
    end
end