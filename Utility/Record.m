classdef Record < handle

    properties
        Flames
        FigureHandle
        FlameCount
    end

    methods
        function obj = Record(FlameNum)
            obj.FlameCount = 0;
            frames(FlameNum) = struct('cdata', [], 'colormap', []);
            obj.Flames = frames;
        end
        function setFigureHandle(obj,pFig)
            obj.FigureHandle = pFig;
            obj.FlameCount = obj.FlameCount + 1;
            obj.Flames(obj.FlameCount) = getframe(obj.FigureHandle);
        end
        function getFlame(obj)
            obj.FlameCount = obj.FlameCount + 1;
            obj.Flames(obj.FlameCount) = getframe(obj.FigureHandle);
        end
        function videowrite(obj, varargin)
            obj.getFlame();
            video = VideoWriter(varargin{:});
            open(video);
            writeVideo(video, obj.Flames(1 : obj.FlameCount));
            close(video);
            % reset flame
            obj.FlameCount = 0;
            obj.Flames(:) = struct('cdata', [], 'colormap', []);
        end
    end
end