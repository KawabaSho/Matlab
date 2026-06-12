classdef class_Graph < handle

    properties
        Vertex % array [q,q,...]
        Edge   % array [[id_1; id_2],...]
    end

    methods
        function obj = class_Graph()
        end
        function [Ver, Num] = getVertex(obj)
            Ver = obj.Vertex;
            Num = size(obj.Vertex,2);
        end
        function [Edge, Num] = getEdge(obj)
            Edge = obj.Edge;
            Num = size(obj.Edge,2);
        end

    end
end