classdef FigureManager < handle
    %  
    properties
        FigureHandle
        PlotList
        ListNum
    end

    methods
        function obj = FigureManager(Name)
            obj.FigureHandle = figure("Name",Name);
            obj.PlotList = {};
            obj.ListNum = 0;
        end

        function addList(obj, Handle)
            % input class : Plot3D

            obj.ListNum = obj.ListNum + 1;
            obj.PlotList{1,obj.ListNum} = Handle;
            
        end
        function setData(obj, Num, Data)
            obj.PlotList{1,Num}.set(Data);
        end
        function addData(obj, Num, Data)
            % s = obj.PlotList{1,Num}.plot3handle
            obj.PlotList{1,Num}.add(Data);
        end

    end
    methods (Static)
        function setGraphic
            grid on
            axis equal
            xlim([-10,10])
            ylim([-10,10])
            zlim([-10,10])
        end
    end
end