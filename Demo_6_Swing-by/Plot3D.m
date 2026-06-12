classdef Plot3D < handle
    properties
        % plot3 handle
        plot3handle
    end
    methods
        function obj = Plot3D(varargin)
            obj.plot3handle = plot3(varargin{:});
        end
        function add(obj, inputArg)
            obj.plot3handle.XData = [obj.plot3handle.XData, inputArg(1)];
            obj.plot3handle.YData = [obj.plot3handle.YData, inputArg(2)];
            obj.plot3handle.ZData = [obj.plot3handle.ZData, inputArg(3)];
        end
        function set(obj, Data)
            obj.plot3handle.XData = Data(1,:);
            obj.plot3handle.YData = Data(2,:);
            obj.plot3handle.ZData = Data(3,:);
        end
    end
end


