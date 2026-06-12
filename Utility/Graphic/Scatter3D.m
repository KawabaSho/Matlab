classdef Scatter3D < handle
    properties
        % plot3 handle
        scatterhandle
    end
    methods
        function obj = Scatter3D(varargin)
            obj.scatterhandle = scatter3(varargin{:});
        end
        function add(obj, inputArg)
            obj.scatterhandle.XData = [obj.scatterhandle.XData, inputArg(1)];
            obj.scatterhandle.YData = [obj.scatterhandle.YData, inputArg(2)];
            obj.scatterhandle.ZData = [obj.scatterhandle.ZData, inputArg(3)];
        end
        function set(obj, Data)
            obj.scatterhandle.XData = Data(1,:);
            obj.scatterhandle.YData = Data(2,:);
            obj.scatterhandle.ZData = Data(3,:);
        end
        function setColorData(obj,data)
            obj.scatterhandle.CData = data;
        end
    end
end