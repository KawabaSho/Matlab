classdef class_Octree < handle
    % 2D

    properties
        Tree2Object
        Object2Tree
        Lenge
        PartitionLevel
        ObjNum

        O_axis
        WW
    end

    methods

        function obj = class_Octree(PartitionLevel,AreaSize,WW)
            % d = x_max(2) - x_min(1) 
            arguments
                PartitionLevel;
                AreaSize;
                WW;
            end
            obj.PartitionLevel  = PartitionLevel;
            obj.Lenge  = AreaSize(:,2) - AreaSize(:,1);
            
            obj.O_axis = AreaSize(:,1);

            obj.Tree2Object = cell(1,(8^(PartitionLevel+1) - 1)/7); % Tree id -> object id
            obj.Object2Tree = {};% Tree id <- object id 
            obj.ObjNum = 0;
            obj.WW = WW;
        end
        function id = getIndex(obj,xyz,userLevel)% vector
            [xyz] = obj.Transform(xyz);
            Box  = fix(xyz./obj.WW); % num
            N = obj.GetMortonNumber(Box(1,:),Box(2,:),Box(3,:));
            sf = 3*(obj.PartitionLevel-userLevel);
            index = bitshift(N,-sf);
            id = int32(index); % liner QSP
        end
        function AddObject(obj,AABB_L,AABB_R)% vector
            [AABB_L] = obj.Transform(AABB_L);
            [AABB_R] = obj.Transform(AABB_R);
            Level    = obj.PartitionLevel;
            LevelNum = Level + 1;

            BoxL  = fix(AABB_L./obj.WW); % num
            BoxR  = fix(AABB_R./obj.WW); % num
            NL = obj.GetMortonNumber(BoxL(1,:),BoxL(2,:),BoxL(3,:));
            NR = obj.GetMortonNumber(BoxR(1,:),BoxR(2,:),BoxR(3,:));
            NG = bitxor(NL,NR);
            for k = 1 : LevelNum % from top
                Itre = k-1;
                sf = 3*(Level-Itre);
                if sf == 1
                    flag = 0;
                else
                    flag = bitget(NG,sf)||bitget(NG,sf-1)||bitget(NG,sf-2);
                end
                if flag
                    break;
                end
            end
            obj.ObjNum = obj.ObjNum + 1;
            index = bitshift(NL,-sf);
            obj.Object2Tree{obj.ObjNum} = [Itre,index];
            id = (8^Itre - 1)/7 + int32(index); % liner QSP
            obj.Tree2Object{id} = obj.ObjNum;

        end
        function xx = Transform(obj,xx) % Vector
            xx = xx - obj.O_axis;
        end

    end
    methods (Static)
        function n = GetMortonNumber(x,y,z)
            n = bitor((bitor(BitSeparate3D(x),...
                (bitshift(BitSeparate3D(y),1)))),...
                  bitshift(BitSeparate3D(z),2));
            function n = BitSeparate3D(n)
                n = bitand(bitor(n,bitshift(n,8)), 0x0000f00f);
                n = bitand(bitor(n,bitshift(n,4)), 0x000c30c3);
                n = bitand(bitor(n,bitshift(n,2)), 0x00249249);
            end
        end
    end
end