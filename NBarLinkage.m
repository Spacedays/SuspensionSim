classdef NBarLinkage
    properties
        Linkage    (2,:)        {mustBeReal} % First row is lengths, second row is angles [radians]; NaN for solution variables
        DrivingVar (1,2) = [2,1]        % coordinate of driving variable
        PosVectors(1,:) logical     % Vectors with a positive sign in the linkage eqn (( e.g. for r1*e^j*th1 + r2*e^j*th2 - r3*e^j*th3 = 0, it would be [1 1 0] ))
        NegVectors(1,:) logical     % Vectors w/ negative sign in linkage eqn ...^          (( for r1*e^j*th1 + r2*e^j*th2 - r3*e^j*th3 = 0, it would be [0 0 1] ))
        UnknownPos1(1,2) = [2 2]        % Row, Col of unknown 1. First row is lengths, second is angles [rads]
        UnknownPos2(1,2) = [2 3]
        LinkageRange  (2,:,:)               % used to store a range of linkage arrangements - row 1&2 are length& angle, columns are the variable numbers, 3 is the step index
        PriorGuesses = [0 0]    % Guess variables for fzero"
        Opt                     % Fzero Options
    end
    % A four bar linkage will have two unknowns, one driving the other.
    
    methods
        function Obj = NBarLinkage(Linkage,DrivingVar,InitGuesses,PosVectors,NegVectors)   % Constructor
            % The unknows are defined by NaN values
            % Example fxn call: linkage = NBarLinkage([r1 r2 r3 r4; Th1 NaN Th3 NaN], [2,3], [10*D2R, 90*D2R]);
                % Here, Th3 is the driving variable with Th2 and Th4 as the solution variables
            % Setup options for fsolve
            Obj.Opt = optimset('Display','off'); % Needed for non-square systems: ,'Algorithm', 'levenberg-marquardt'
            
            % Set obj.priorGuesses if provided
            if (nargin >= 3)
               Obj.PriorGuesses = InitGuesses;
            end
            
            % Determine the array row & col of unknowns
            unknowns = isnan(Linkage);
            if (length(unknowns(unknowns==1)) > 2)
                disp("Too many unknowns in linkage matrix")
                return
            end
            Obj.Linkage = Linkage;
            Obj.DrivingVar = DrivingVar;
            
            % ID unknown coordinates
            if (isempty(find(unknowns(1,:),1)))       % ==> Both unknowns are angles
                unknownIndexes = find(unknowns(2,:));
                Obj.UnknownPos1 = [2, unknownIndexes(1)];
                Obj.UnknownPos2 = [2, unknownIndexes(2)];
            elseif (isEmpty(find(unknowns(2,:))))   % ==> Both unknowns are lengths
                unknownIndexes = find(unknowns(2,:));
                Obj.UnknownPos1 = [1, unknownIndexes(1)];
                Obj.UnknownPos2 = [1, unknownIndexes(2)];
            else                                    % ==> One unk Angle, one unk Length
                Obj.UnknownPos1 = [1, find(unknowns(1,:))];     %length
                Obj.UnknownPos2 = [2, find(unknowns(2,:))];     %angle
            end
            
            if (exist('PosVectors','var') && exist('NegVectors','var'))
                Obj.PosVectors = PosVectors;
                Obj.NegVectors = NegVectors;
            else
                if (max(size(Obj.Linkage)) == 4)
                    Obj.PosVectors = logical([1 1 0 0]);
                elseif (max(size(Obj.Linkage)) == 3)
                    Obj.PosVectors = logical([1 1 0]);
                else
                    disp(max(size(Obj.Linkage)))
                end
                Obj.NegVectors = ~Obj.PosVectors;
            end
            
        end
        
        function out = LinkageEqn(Obj, X)      % Function defined using geometry values from the obj's linkage array
            % Updates the unknowns' position to the previous guess
            Obj.Linkage(Obj.UnknownPos1(1),Obj.UnknownPos1(2)) = X(1);%obj.priorGuesses(1);
            Obj.Linkage(Obj.UnknownPos2(1),Obj.UnknownPos2(2)) = X(2);%obj.priorGuesses(2);
            if (max(size(Obj.Linkage)) == 5)
            eqn = Obj.Linkage(1,1).*exp(1i.*Obj.Linkage(2,1)) + Obj.Linkage(1,2).*exp(1i.*Obj.Linkage(2,2)) ...
                 - Obj.Linkage(1,3).*exp(1i.*Obj.Linkage(2,3)) - Obj.Linkage(1,4).*exp(1i.*Obj.Linkage(2,4));
            end
            if (max(size(Obj.Linkage)) == 4)
            eqn = Obj.Linkage(1,1).*exp(1i.*Obj.Linkage(2,1)) + Obj.Linkage(1,2).*exp(1i.*Obj.Linkage(2,2)) ...
                 - Obj.Linkage(1,3).*exp(1i.*Obj.Linkage(2,3)) - Obj.Linkage(1,4).*exp(1i.*Obj.Linkage(2,4));
            end
            if (max(size(Obj.Linkage)) == 3)
            eqn = Obj.Linkage(1,1).*exp(1i.*Obj.Linkage(2,1)) + Obj.Linkage(1,2).*exp(1i.*Obj.Linkage(2,2)) ...
                 - Obj.Linkage(1,3).*exp(1i.*Obj.Linkage(2,3));
            end
            out = [real(eqn); imag(eqn)];
        end
        
        function out = LinkageEqn2(Obj, X)      % Function defined using geometry values from the obj's linkage array
            % Vectorized LinkageEqn
            
            % Updates the unknowns' position to the previous guess
            Obj.Linkage(Obj.UnknownPos1(1),Obj.UnknownPos1(2)) = X(1);%obj.priorGuesses(1);
            Obj.Linkage(Obj.UnknownPos2(1),Obj.UnknownPos2(2)) = X(2);%obj.priorGuesses(2);
            
            eqn = sum(Obj.Linkage(1,Obj.PosVectors).*exp(1i.*Obj.Linkage(2,Obj.PosVectors)) ) - sum(Obj.Linkage(1,Obj.NegVectors).*exp(1i.*Obj.Linkage(2,Obj.NegVectors)) );
            
            out = [real(eqn); imag(eqn)];
        end
        
        function out = LinkageBCEqn(Obj, X)      % Fxn using Obj.linkage and a boundary condition eqn
            out = linkageEqn(Obj,X);
            
        end
        
        function MakeLinkageBCEqn()
        
        end
        
        function [unknown1, unknown2] = CalcLinkage(Obj,drivingVarValue,DrivingVar) % Use this function to calculate a linkage position
            if (nargin < 2)
                disp("ERROR: Not Enough Argiments")
                return
            elseif (nargin > 2)
                Obj.DrivingVar = DrivingVar;
            end
            Obj.Linkage(Obj.DrivingVar(1),Obj.DrivingVar(2)) = drivingVarValue;     % Change driving Var Value as requested
            
%             Xtemp = fsolve(@(X) LinkageEqn(Obj,X), Obj.PriorGuesses, Obj.Opt);
            Xtemp = fsolve(@(X) LinkageEqn2(Obj,X), Obj.PriorGuesses, Obj.Opt);
            unknown1 = Xtemp(1);
            unknown2 = Xtemp(2);
        end
        
        function [solV1, solV2] = CalcLinkageRange(Obj,drivingVarVector,options)
            % Returns an mxn vector of the entire linkage for the analysis range
            arguments
                Obj
                drivingVarVector (1,:) {mustBeNumeric}      
                options.DrivingVar (1,2) {mustBeNumeric}    
                options.FullSoltn (1,1) {mustBeNumeric} = 0
            end
            
            if (nargin < 2)
                disp("ERROR: Not Enough Argiments")
                return
            elseif (nargin > 2 && ~isempty(options.DrivingVar))
                Obj.DrivingVar = options.DrivingVar;
            end
            
            
            
            solVectors = ones(2,length(drivingVarVector));
                for k=1:(length(drivingVarVector))
                    % Solving for Position
                    DrivingVarValue = drivingVarVector(k);
                    [sol1,sol2] = CalcLinkage(Obj,DrivingVarValue);
                    solVectors(1,k) = sol1;                  % Store solutions to vector
                    solVectors(2,k) = sol2;

                    Obj.PriorGuesses(1) = sol1;    % Use last solution for next guess
                    Obj.PriorGuesses(2) = sol2;
                end
                
            if ~options.FullSoltn
                solV1 = solVectors(1,:);
                solV2 = solVectors(2,:);
            else
                % Setup initial vectors
                rVectors = ones(length(Obj.Linkage), length(drivingVarVector));
                thVectors = ones(length(Obj.Linkage), length(drivingVarVector));
                for k=1:length(Obj.Linkage)
                    rVectors(k,:) = rVectors(k) .* Obj.Linkage(1,k);
                    thVectors(k,:) = thVectors(k) .* Obj.Linkage(2,k);
                end
                
                %solVectors = CalcSoltnRange(obj,DrivingVarVector);
                % Put the solution vectors in the right place
                switch Obj.UnknownPos1(1)
                    case 1
                        rVectors(Obj.UnknownPos1(2),:) = solVectors(1,:);
                    case 2
                        thVectors(Obj.UnknownPos1(2),:) = solVectors(1,:);
                end
                switch Obj.UnknownPos2(1)
                    case 1
                        rVectors(Obj.UnknownPos2(2),:) = solVectors(2,:);
                    case 2
                        thVectors(Obj.UnknownPos2(2),:) = solVectors(2,:);
                end
                solV1 = rVectors;
                solV2 = thVectors;
            end
            
            
        end

        function [solV1, solV2] = CalcChangingLinkage(Obj,drivingLinkageVector,options)
            % Returns an mxn vector of the entire linkage for the analysis range
            arguments
                Obj
                drivingLinkageVector (:,:) {mustBeNumeric}      % MxN vector of M links over range N, with two vectors of NaN
                options.DrivingVar (1,2) {mustBeNumeric}        % Allows changing the driving variable on the object
                options.FullSoltn (1,1) {mustBeNumeric} = 0
            end
            
            if (nargin < 2)
                disp("ERROR: Not Enough Argiments")
                return
            elseif (nargin > 2 && ~isempty(options.DrivingVar))
                Obj.DrivingVar = options.DrivingVar;
            end
            
            Obj.LinkageRange = drivingLinkageVector;
            
            drivingVarVector = Obj.LinkageRange(Obj.DrivingVar(1), Obj.DrivingVar(2), :);
                for k=1:(length(drivingVarVector))
                    % Solving for Position
                    [sol1,sol2] = CalcLinkage(Obj,drivingVarVector(k));
                    Obj.LinkageRange(Obj.UnknownPos1(1),Obj.UnknownPos1(2),k) = sol1;                  % Store solutions to vector
                    Obj.LinkageRange(Obj.UnknownPos2(1),Obj.UnknownPos2(2),k) = sol2;   % solVectors(2,k) = sol2;

                    Obj.PriorGuesses(1) = sol1;    % Use last solution for next guess
                    Obj.PriorGuesses(2) = sol2;
                end
                
            if ~options.FullSoltn           %!!!
                solV1 = Obj.LinkageRange(Obj.UnknownPos1(1),Obj.UnknownPos1(2),:);
                solV2 = Obj.LinkageRange(Obj.UnknownPos2(1),Obj.UnknownPos2(2),:);
            else
                solV1 = Obj.LinkageRange(1,:,:);
                solV2 = Obj.LinkageRange(2,:,:);
%                 % Setup initial vectors
%                 rVectors = ones(length(Obj.Linkage), length(drivingLinkageVector));
%                 thVectors = ones(length(Obj.Linkage), length(drivingLinkageVector));
%                 for k=1:length(Obj.Linkage)
%                     rVectors(k,:) = rVectors(k) .* Obj.Linkage(1,k);
%                     thVectors(k,:) = thVectors(k) .* Obj.Linkage(2,k);
%                 end
%                 
%                 %solVectors = CalcSoltnRange(obj,DrivingVarVector);
%                 % Put the solution vectors in the right place
%                 switch Obj.UnknownPos1(1)
%                     case 1
%                         rVectors(Obj.UnknownPos1(2),:) = solVectors(1,:);
%                     case 2
%                         thVectors(Obj.UnknownPos1(2),:) = solVectors(1,:);
%                 end
%                 switch Obj.UnknownPos2(1)
%                     case 1
%                         rVectors(Obj.UnknownPos2(2),:) = solVectors(2,:);
%                     case 2
%                         thVectors(Obj.UnknownPos2(2),:) = solVectors(2,:);
%                 end
%                 solV1 = rVectors;
%                 solV2 = thVectors;
            end
            
            disp("Done")
        end
        
    end
end