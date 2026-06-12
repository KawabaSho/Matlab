classdef PhysicsManeger < handle
    properties
        RigidBodies
        ActorNum
    end

    methods
        function obj = PhysicsManeger()
            obj.RigidBodies = {};
            obj.ActorNum = 0;
        end

        function addActor(obj,RigitBody_handle)
            % class of argument : RigitBody
            obj.ActorNum = obj.ActorNum + 1;
            obj.RigidBodies{obj.ActorNum} = RigitBody_handle;
        end

        function Integration(obj, input, dt)
            for i = 1:obj.ActorNum
                obj.RigidBodies{i}.Integration(input, dt);
            end
        end

        function value = getValue(obj,Num)
            value = obj.RigidBodies{Num}.getValue;
        end

    end
end