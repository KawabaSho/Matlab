classdef PhysicsManeger < handle
    properties
        PhysicsHandles
        ActorNum
    end

    methods
        function obj = PhysicsManeger()
            obj.PhysicsHandles = {};
            obj.ActorNum = 0;
        end

        function addActor(obj,RigitBody_handle)
            % class of argument : RigitBody
            obj.ActorNum = obj.ActorNum + 1;
            obj.PhysicsHandles{obj.ActorNum} = RigitBody_handle;
        end

        function Integration(obj, input, dt)
            for i = 1:obj.ActorNum
                obj.PhysicsHandles{i}.Integrate(input, dt);
            end
        end

        function value = getValue(obj,Num)
            value = obj.PhysicsHandles{Num}.getValue;
        end

    end
end