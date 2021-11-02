classdef NBarLinkage
    properties
        linkage(2,:)        {mustBeReal} % First row is lengths, second row is angles [radians]; NaN for solution variables
        drivingVar = [2,1]      % coordinate of driving variable
        unknownPos1 = [2 2]     % Row, Col of unknown 1. First row is lengths, second is angles [rads]
        unknownPos2 = [2 3]
        priorGuesses = [0 0]    % Guess variables for fzero
        opt                     % Fzero Options
    end
    % A four bar linkage will have two unknowns, one driving the other.
    
    methods
        function obj = NBarLinkage(linkage,drivingVar,initGuesses)   % Constructor
            % The unknows are defined by NaN values
            % Example fxn call: linkage = NBarLinkage([r1 r2 r3 r4; Th1 NaN Th3 NaN], [2,3], [10*D2R, 90*D2R]);
                % Here, Th3 is the driving variable with Th2 and Th4 as the solution variables
            
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
            if (isempty(find(unknowns(1,:),1)))       % ==> Both unknowns are angles
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
        
        function out = LinkageEqn(obj, X)      % Function defined using geometry values from the ojb's linkage array
            % Updates the unknowns' position to the previous guess
            obj.linkage(obj.unknownPos1(1),obj.unknownPos1(2)) = X(1);%obj.priorGuesses(1);
            obj.linkage(obj.unknownPos2(1),obj.unknownPos2(2)) = X(2);%obj.priorGuesses(2);
            if (max(size(obj.linkage)) == 5)
            eqn = obj.linkage(1,1).*exp(1i.*obj.linkage(2,1)) + obj.linkage(1,2).*exp(1i.*obj.linkage(2,2)) ...
                 - obj.linkage(1,3).*exp(1i.*obj.linkage(2,3)) - obj.linkage(1,4).*exp(1i.*obj.linkage(2,4));
            end
            if (max(size(obj.linkage)) == 4)
            eqn = obj.linkage(1,1).*exp(1i.*obj.linkage(2,1)) + obj.linkage(1,2).*exp(1i.*obj.linkage(2,2)) ...
                 - obj.linkage(1,3).*exp(1i.*obj.linkage(2,3)) - obj.linkage(1,4).*exp(1i.*obj.linkage(2,4));
            end
            if (max(size(obj.linkage)) == 3)
            eqn = obj.linkage(1,1).*exp(1i.*obj.linkage(2,1)) + obj.linkage(1,2).*exp(1i.*obj.linkage(2,2)) ...
                 - obj.linkage(1,3).*exp(1i.*obj.linkage(2,3));
            end
            out = [real(eqn); imag(eqn)];
        end
        
        function [unknown1, unknown2] = CalcLinkage(obj,drivingVarValue,drivingVar) % Use this function to calculate a linkage position
            if (nargin < 2)
                disp("ERROR: Not Enough Argiments")
                return
            elseif (nargin > 2)
                obj.drivingVar = drivingVar;
            end
            obj.linkage(obj.drivingVar(1),obj.drivingVar(2)) = drivingVarValue;     % Change driving Var Value as requested
            
%             fxn = ;
%             a = @(X) LinkageEqn(obj,X);
            Xtemp = fsolve(@(X) LinkageEqn(obj,X), obj.priorGuesses, obj.opt);
            unknown1 = Xtemp(1);
            unknown2 = Xtemp(2);
        end
        
        function solVectors = CalcSoltnRange(obj,drivingVarVector,drivingVar)
            % Returns a 2xlen(drivingVarVector) vector for the two solution vars
            
            if (nargin < 2)
                disp("ERROR: Not Enough Argiments")
                return
            elseif (nargin > 2)
                obj.drivingVar = drivingVar;
            end
            
%             Setup initial vectors
            solVectors = ones(2,length(drivingVarVector));
            for k=1:(length(drivingVarVector))
                % Solving for Position
            	drivingVarValue = drivingVarVector(k);
            	[sol1,sol2] = CalcLinkage(obj,drivingVarValue);
            	solVectors(1,k) = sol1;                  % Store solutions to vector
            	solVectors(2,k) = sol2;

            	obj.priorGuesses(1) = sol1;    % Use last solution for next guess
            	obj.priorGuesses(2) = sol2;
            end
        end
        
        function [rVectors, thVectors] = CalcLinkageRange(obj,drivingVarVector,drivingVar)
            % Returns vectors for the entire linkage
            
            if (nargin < 2)
                disp("ERROR: Not Enough Argiments")
                return
            elseif (nargin > 2)
                obj.drivingVar = drivingVar;
            end
            
            % Setup initial vectors
            rVectors = ones(length(obj.linkage), length(drivingVarVector));
            thVectors = ones(length(obj.linkage), length(drivingVarVector));
            for k=1:length(obj.linkage)
                rVectors(k,:) = rVectors(k) .* obj.linkage(1,k);
                thVectors(k,:) = thVectors(k) .* obj.linkage(2,k);
            end
            
%             switch obj.drivingVar(1)
%                 case 1
%                     rVectors(obj.drivingVar(1),:) = drivingVarVector;
%                 case 2
%                     thVectors(obj.drivingVar(1),:) = drivingVarVector;
%             end
            
            solVectors = CalcSoltnRange(obj,drivingVarVector);
            
            % Set 
            switch obj.unknownPos1(1)
                case 1
                    rVectors(obj.unknownPos1(2),:) = solVectors(1,:);
                case 2
                    thVectors(obj.unknownPos1(2),:) = solVectors(1,:);
            end
            switch obj.unknownPos2(1)
                case 1
                    rVectors(obj.unknownPos2(2),:) = solVectors(2,:);
                case 2
                    thVectors(obj.unknownPos2(2),:) = solVectors(2,:);
            end
            
            
            
        end

    end
end