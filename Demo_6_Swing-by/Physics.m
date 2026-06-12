classdef Physics < handle
    %  values : [xx_', vv_', ww_', qq_', ]'
    %  flag  : active or inactive [collision, calculation, ]
    %  constraints : joint information
    %  information : dynamic/static
    %  class : Index.m
    properties 
        values
        Integrator
    end
    
    methods
        function obj = Physics(values, Integrator)
            % values      :  state variables
            obj.values      = values;
            obj.Integrator = Integrator;
        end
        function Integrate(obj,uu,dt)
            obj.values = obj.Integrator(obj.values, uu, dt);
        end
        function value = getValue(obj)
            value = obj.values;
        end
    end
end







