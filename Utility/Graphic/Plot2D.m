classdef Plot2D < handle
    properties
        % plot2 handle
        plot2handle
    end
    methods
        function obj = Plot2D(varargin)
            obj.plot2handle = ...
                plot(varargin{:});
        end
        function add(obj, inputArg)
            obj.plot2handle.XData = [obj.plot2handle.XData, inputArg(1)];
            obj.plot2handle.YData = [obj.plot2handle.YData, inputArg(2)];
        end
        function set(obj, inputArg)
            obj.plot2handle.XData = inputArg(1,:);
            obj.plot2handle.YData = inputArg(2,:);
        end
    end
end


