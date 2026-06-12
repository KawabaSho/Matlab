classdef class_Quaternary_Spatial_Partitiong < handle
    % 2D

    properties
        Tree2Object
        Object2Tree
        Lenge
        ParticleSize
        PartitionLevel
        ObjNum
    end

    methods

        function obj = class_Quaternary_Spatial_Partitiong(PartitionLevel,AreaSize,ParticleSize)
            % d = x_max(2) - x_min(1) 
            arguments
                PartitionLevel;
                AreaSize;
                ParticleSize;
            end
            obj.PartitionLevel  = PartitionLevel;
            obj.Lenge  = AreaSize(1,2) - AreaSize(1,1);
            obj.Tree2Object = cell(1,(4^(PartitionLevel+1) - 1)/3); % Tree id -> object id
            obj.Object2Tree = {};% Tree id <- object id 
            obj.ObjNum = 0;
            obj.ParticleSize = ParticleSize;
        end

        function AddObject(obj,x,y)
            Level    = obj.PartitionLevel;
            LevelNum = Level + 1;
            L        = obj.Lenge;

            BoxPos = obj.GetBoxPosition(x,y,obj.ParticleSize); % pos
            U = L/2^Level; % unit
            BoxL  = fix(BoxPos{1}./U); % num
            BoxR  = fix(BoxPos{2}./U); % num
            NL = obj.GetMortonNumber(BoxL(1,:),BoxL(2,:));
            NR = obj.GetMortonNumber(BoxR(1,:),BoxR(2,:));
            NG = bitxor(NL,NR);
            for k = 1 : LevelNum % from top
                Itre = k-1;
                sf = 2*(Level-Itre);
                if sf == 0
                    flag = 0;
                else
                    flag = bitget(NG,sf)||bitget(NG,sf-1);
                end
                if flag
                    break;
                end
            end
            obj.ObjNum = obj.ObjNum + 1;
            index = bitshift(NL,-sf);
            obj.Object2Tree{obj.ObjNum} = [Itre,index];
            id = (4^Itre - 1)/3 + int32(index); % liner QSP
            obj.Tree2Object{id} = obj.ObjNum;

        end

    end
    methods (Static)


        function n = GetMortonNumber(x,y)
            n = bitor(BitSeparate32(x),...
                (bitshift(BitSeparate32(y),1)));
            function n = BitSeparate32(n)
            % 16bit
            n = bitand(bitor(n,bitshift(n,8)), 0x00ff00ff);
            n = bitand(bitor(n,bitshift(n,4)), 0x0f0f0f0f);
            n = bitand(bitor(n,bitshift(n,2)), 0x33333333);
            n = bitand(bitor(n,bitshift(n,1)), 0x55555555);
            end
        end
        function Box = GetBoxPosition(center_x,center_y,L)
            r = 0.5*L;
            c = [center_x;center_y];
            Box{1} = c+[-r;r];
            Box{2} = c+[r;-r];
        end

    end
end