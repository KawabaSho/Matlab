classdef Patch < handle
    properties
        hPatch
        position
        coord
        idx
        idx_num

        Vertices
        VerticesSize
    end
    methods
        function obj = Patch(position, coord, idx, varargin)
            % position : Object's position from an local coordinate system perspective.
            % coord    : Point cloud
            % idx      : Point id
            % varargin : Arguments of patch function

            obj.hPatch = patch(varargin{:});
            obj.position = position;
            obj.coord = coord;
            obj.idx   = idx;
            obj.idx_num = numel(idx);
            obj.Vertices = obj.hPatch.Vertices;
            obj.VerticesSize = size(obj.hPatch.Vertices);
        end
        function Refresh(obj, RR, xx_)
            obj.hPatch.Vertices = (obj.coord + obj.position')*RR' + xx_';
        end
        function Move(obj, xx_)
            obj.hPatch.Vertices = obj.Vertices + xx_';
        end
        function Rotate(obj, RR, xx_, RR_rotate)
            obj.hPatch.Vertices = (obj.coord*RR_rotate' + obj.position')*RR' + xx_';
        end
        function setVertices(obj, Vertices)
            obj.hPatch.Vertices = Vertices;
        end
        function setCData(obj, Data)
            obj.hPatch.FaceColor = Data;
        end
    end
end