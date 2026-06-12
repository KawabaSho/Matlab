classdef Scatter < handle
    properties
        % plot3 handle
        scatterhandle
    end
    methods
        function obj = Scatter(varargin)
            obj.scatterhandle = scatter(varargin{:});
        end
        function add(obj, inputArg)
            obj.scatterhandle.XData = [obj.scatterhandle.XData, inputArg(1)];
            obj.scatterhandle.YData = [obj.scatterhandle.YData, inputArg(2)];
        end
        function set(obj, Data)
            obj.scatterhandle.XData = Data(1,:);
            obj.scatterhandle.YData = Data(2,:);
        end
        function setColorData(obj,data)
            obj.scatterhandle.CData = data;
        end
    end
end