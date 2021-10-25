classdef FourBarLinkage
    properties
        Lengths(1,4) double {mustBeReal, mustBeFinite}
        changingLengths(1,:) {mustBeInteger}            % Bool 1:4 array
        Angles(1,4) double {mustBeReal, mustBeFinite    % Angle in radians
        changingAngles(1,:) {mustBeInteger}             % Bool 1:4 array
    end
    
    % A four bar linkage will have two unknowns, one driving the other.
    
    methods
        function obj = FourBarLinkage(Lengths,changingLengths,Angles,changingAngles)
            if (length(changingLengths(changingLengths > 0)) + length(changingLengths(changingLengths > 0)) > 2)
                disp("Too many unknowns")
                return
            end
            obj.Lengths = Lengths;
            obj.changingLengths = changingLengths;
            obj.Angles = Angles;
            obj.changingAngles = changingAngles;
        end
        function MoveLinkage(move1or2,value)
            if (move1or2 == 1)
                
            end
        end
        function plotElement(obj,edgeColor,faceColor,edgeAlpha,faceAlpha)
        end
        
    end
end