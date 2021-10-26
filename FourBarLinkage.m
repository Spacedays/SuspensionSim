classdef FourBarLinkage
    properties
        linkage(2,4)        {mustBeReal} % First row is lengths, second row is angle [radians]; NaN for solution variables
        drivingVar = [2,1]      % coordinate of driving variable
        unknownPos1 = [2 2]     % Row, Col of unknown 1. First row is length, second is angle [rads]
        unknownPos2 = [2 3]
        priorGuesses = [0 0]    % Guess variables for fzero
        opt                     % Fzero Options
    end
    % A four bar linkage will have two unknowns, one driving the other.
    
    methods
        function obj = FourBarLinkage(linkage,drivingVar,initGuesses)   % Constructor
            % The unknows are defined by NaN values
            
            % Setup options for fsolve
            obj.opt = optimset('Display','off');
            
            % Set obj.priorGuesses if provided
            if (nargin == 3)
               obj.priorGuesses = initGuesses;
            end
            
            % Determine the array row & col of unknowns
            unknowns = isnan(linkage);
            if (length(unknowns(unknowns==1)) > 2)
                disp("Too many unknowns in linkage matrix")
                return
            end
            obj.linkage = linkage;
            obj.drivingVar = drivingVar;
            
            % ID unknown coordinates
            if (isempty(find(unknowns(1,:))))       % ==> Both unknowns are angles
                unknownIndexes = find(unknowns(2,:));
                obj.unknownPos1 = [2, unknownIndexes(1)];
                obj.unknownPos2 = [2, unknownIndexes(2)];
            elseif (isEmpty(find(unknowns(2,:))))   % ==> Both unknowns are lengths
                unknownIndexes = find(unknowns(2,:));
                obj.unknownPos1 = [1, unknownIndexes(1)];
                obj.unknownPos2 = [1, unknownIndexes(2)];
            else                                    % ==> One unk Angle, one unk Length
                obj.unknownPos1 = [1, find(unknowns(1,:))];     %length
                obj.unknownPos2 = [2, find(unknowns(2,:))];     %angle
            end
        end
        
        function [unknown1, unknown2] = CalcLinkage(obj,drivingVarValue,drivingVar) % Use this function to calculate a linkage position
            if (nargin < 2)
                disp("ERROR: Not Enough Argiments")
                return
            elseif (nargin == 2)
                drivingVar = obj.drivingVar;
            else
                obj.drivingVar = drivingVar;
            end
            obj.linkage(drivingVar(1),drivingVar(2)) = drivingVarValue;     % Change driving Var Value as requested
            
%             fxn = ;
%             a = @(X) LinkageEqn(obj,X);
            Xtemp = fsolve(@(X) LinkageEqn(obj,X), obj.priorGuesses, obj.opt);
            unknown1 = Xtemp(1);
            unknown2 = Xtemp(2);
        end
        
        function out = LinkageEqn(obj, X)      % Function defined using the geometry values defined 
            % Updates the unknowns' position to the previous guess
            obj.linkage(obj.unknownPos1(1),obj.unknownPos1(2)) = X(1);%obj.priorGuesses(1);
            obj.linkage(obj.unknownPos2(1),obj.unknownPos2(2)) = X(2);%obj.priorGuesses(2);
            
            eqn = obj.linkage(1,1).*exp(1i.*obj.linkage(2,1)) + obj.linkage(1,2).*exp(1i.*obj.linkage(2,2)) ...
                 - obj.linkage(1,3).*exp(1i.*obj.linkage(2,3)) - obj.linkage(1,4).*exp(1i.*obj.linkage(2,4));
            out = [real(eqn); imag(eqn)];
        end
        
    end
end