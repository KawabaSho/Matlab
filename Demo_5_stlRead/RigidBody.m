classdef RigidBody < handle
    %  values : [xx_', vv_', ww_', qq_', ]'
    %  flag  : active or inactive [collision, calculation, ]
    %  constraints : joint information
    %  information : dynamic/static
    %  class : Index.m
    properties 
        dynamics
        values
        flag
        constraints
        information
    end
    
    methods
        function obj = RigidBody(values, dynamics, constraints, information)
            % values      :  state variables
            % dynamics    :  dxdt = dynamics(x,u)
            % constrains  :  Cooming soon!
            % imformation :  user interface
            obj.dynamics    = dynamics;
            obj.values      = values;
            obj.constraints = constraints;
            obj.information = information;
        end
        function collision_cheak(obj, collision)
            % Coming soon!
            % information->size
            % 1. section cheak
            % 2. select other obstacles
            % 3. calculation distance
            obj.flag.collision = collision;
        end
        function activity_cheak(obj, values, input)
            % calculation flag
            % input : force or torque
            if any((input(:).^2 > 1e-24))
                bo = true;
            else
                bo = false;
            end
            v = values(Index.vv);
            if any((v(:).^2 > 1e-24))
                bo = true;
            end
            w = values(Index.ww);
            if any((w(:).^2 > 1e-24))
                bo = true;
            end
            obj.flag.activitey = bo;
        end
        function Integration(obj, uu_, dt)
            % integration
            x0 = obj.values;
            u0 = uu_;
            k1 = obj.dynamics(x0,u0).*dt;
            k2 = obj.dynamics(x0 + 0.5*k1,u0).*dt;
            k3 = obj.dynamics(x0 + 0.5*k2,u0).*dt;
            k4 = obj.dynamics(x0 + k3,u0).*dt;

            % collision cheak
            % Coming soon!

            % dynamic/static cheak
            % Coming soon!

            obj.values = x0 + (k1+2*(k2+k3)+k4)/6;
        end

        function value = getValue(obj)
            value = obj.values;
        end
    end
end







