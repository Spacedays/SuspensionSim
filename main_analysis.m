
clc;clear;%close all
%% misc utility
zci = @(v) find(diff(sign(v))); % zero crossing fcn - returns idx of zero cross
%% Plot setup
% 1-Camber, 2-Scrub, 3-Pullrod Motion Ratio, 4-Pullrod length
SPRows = 2; SPCols = 2;
sub_idx = 1;
figure(1); 

%% Roll Center Calculations
% TBD


%% Kinematics Setup
numpts = 40

D2R = pi/180;    %deg 2 radians
R2D = 180/pi;

% Suspension Linkage: Theta 3 drives the linkage, Theta 1 and all lengths are fixed.

Th1 = 63.0056*D2R; %63.01*D2R;% acos(80.83/137.5);
Th3 = 0;   % Driving Input Starting Value ((var will be iterated over))
r1 = 160.3174; %6.282*25.4;    % mm
r2 = 274.4; %12*25.4;
r3 = 343; %15*25.4;
r4 = 190.5; %7.5*25.4;

initGuesses = [10*D2R,90*D2R]; %Theta 2 and Theta 4 initial guesses
SL = NBarLinkage([r1 r2 r3 r4; Th1 NaN Th3 NaN], [2,3], initGuesses, PosVectors=[1 1 0 0]);   % Make SuspensionLinkage
VTh3 = linspace(-6,6,numpts)*D2R; % VTh3 is array of Theta3 angles to iterate over
bump = r3*sin(VTh3);

% Solve Kinematics using CalcLinkageRange
[~, thVectors] = CalcLinkageRange(SL,VTh3,fullSoltn=1);    % rVectors not used here
VTh1 = thVectors(1,:);
VTh2 = thVectors(2,:);
% VTh3 = thVectors(3,:);    % Th3 is already set so no need to re-set it, but it remains the same so it doesn't matter either way
VTh4 = thVectors(4,:);

%% Linkage Vectors
A = (216.23 + 1i*5*25.4)*ones(1,length(VTh3));
R1 = r1.*exp(1i.*VTh1.*ones(1,length(VTh3)) );  B = A+R1;
R2 = r2*ones(1,length(VTh3)).*exp(1i.*VTh2);    D = B+R2;
R3 = r3*ones(1,length(VTh3)).*exp(1i.*VTh3);    C = A+R3;
R4 = r4*ones(1,length(VTh3)).*exp(1i.*VTh4);    

%% Camber Calcs
figure(1)
subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
camber = atand((real(D)-real(C))./(imag(D)-imag(C)));
plot(bump,camber)
grid on
title("Camber vs Wheel Displacement")
xlabel("Wheel Displacement [mm]"), ylabel("Camber [deg]")
xticks(-40:10:40), yticks(-4:1:3)

%% Scrub Calcs
wheel_dx = (1+3.5)*25.4;
wheel_r = 10*25.4;
WM = C + .5 .* R4; % Wheel Mount Location
WC = WM + wheel_dx .* exp(1i.*(VTh4 - pi/2));    % Wheel centerline
WCP = WC + wheel_r .* exp(1i.* (VTh4 + pi));

% figure(2)
subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
% plot(VTh3*R2D, (real(WCP)-min(real(WCP))));
plot(bump, (real(WCP) - real(WCP(zci(bump))) ) )
% title("Wheel Scrub vs. \theta_3")
title("Wheel Scrub vs. Wheel Displacement")
% xlabel("\theta_3 [deg]"), ylabel("Wheel Scrub [mm]")
xlabel("Wheel Displacement [mm]"), ylabel("Wheel Scrub [mm]")
xticks(-60:10:60), yticks(-15:2.5:5)
grid on


%% Pull Rod Calculations
rocker_axis_pos = [15,27]; %[-15.4,40.8];   % x,y wrt lower right A-arm mount [mm] 
O_R = A + ones(1,length(VTh3))*(rocker_axis_pos(1) + 1i*rocker_axis_pos(2));    % Rocker axis of rotation

rocker_pull_radius = 30;   	% Distance from the rocker axis to the pullrod [mm]
l_pullrod = 270;            % Pullrod length [mm]
pullrod_upright_off = [-65 -20.64];    % pullrod mount x,y offset from upper upright suspension link (point D)
pullrod_Aarm_r = 65;            % Length offset from the upper A-arm mount along the plane of the upper A-Arm
pullrod_Aarm_vert = -20.64;          % Perpendicular offset of the pullrod A-Arm mount; measured perpendicular to the Upper A-Arm


% rocker_axis_pos, O_R, rocker_pull_radius, l_pullrod, pullrod_upright_off, pullrod_Aarm_r, pullrod_Aarm_vert


% E = D + pullrod_upright_off(1) + imag(pullrod_upright_off(2);      % Pull Rod Attachment Point (point E) based on upright
R5 = -pullrod_Aarm_r * exp(1i*VTh2) + pullrod_Aarm_vert * exp(1i*(VTh2+pi/4));
E = D + R5;

%% motion solver - New Fast & correct linkage vector
vers2 = true

if vers2
% Setup drivingLinkageVector Vers 2

V_rP3 = sqrt(real(E-O_R).^2 + imag(E-O_R).^2);
V_thP3 = atan2(imag(E-O_R), real(E-O_R));
V_rP2 = rocker_pull_radius .* ones(size(V_rP3));
V_thP2 = NaN(size(V_rP3));
V_rP1 = l_pullrod .* ones(size(V_rP3));
V_thP1 = NaN(size(V_rP3));

initGuesses2 = [20*D2R, 90*D2R];  % ((Th6 Th7))

drivingLinkageVector2 = zeros(2,3,length(V_rP3));
drivingLinkageVector2(1,:,:) = [V_rP1; V_rP2; V_rP3];
drivingLinkageVector2(2,:,:) = [V_thP1; V_thP2; V_thP3];

% Linkage Loop Eqn: F->E, F->O_R, O_R->E
RL2 = NBarLinkage([V_rP1(1) V_rP2(1) V_rP3(1); NaN NaN V_thP3(1)], [1,3], initGuesses2, PosVectors=[0 1 1]);         % Make RockerLinkage


[VTh6, VTh7] = CalcChangingLinkage(RL2,drivingLinkageVector2);

Th5 = atan2(pullrod_upright_off(2),pullrod_upright_off(1));  % Angle of the vector from upper upright susp. link to pullrod mount
VTh5 = Th5*ones(size(V_rP3));
VTh6 = squeeze(VTh6)';
VTh7 = squeeze(VTh7)';

r6 = V_rP1;
r7 = V_rP2;

R6 = r6.*exp(1i.*VTh6);
R7 = r7.*exp(1i*VTh7);
R8 = D - O_R;

%Two methods for finding point F - if one fails for an unknown reason try the other
F = E - R6;
F2 = O_R - R7;

else
%% -- Original (bad & slow) method using 4 vectors
% Setup drivingLinkageVector
Vr8 = ((real(D)-real(O_R)).^2+(imag(D)-imag(O_R)).^5).^0.5;   r8 = Vr8(1);  % Distance from rocker axis to upper upright susp. link
r5 = sum(pullrod_upright_off.^2).^0.5;  Vr5 = r5*ones(size(Vr8));           % Distance from upper upright suspension link to pullrod mount
r6 = l_pullrod; Vr6 = r6*ones(size(Vr8));
r7 = rocker_pull_radius;    Vr7 = r7*ones(size(Vr8));
Th5 = atan2(pullrod_upright_off(2),pullrod_upright_off(1));  % Angle of the vector from upper upright susp. link to pullrod mount
if (isnan(Th5)), Th5 = 0; elseif (Th5 < 0), Th5 = Th5 + 360*D2R; end
VTh5 = Th5*ones(size(Vr8));
VTh6 = NaN*ones(size(Vr8));
VTh7 = VTh6;
VTh8 = atan2((imag(D)-imag(O_R)),(real(D)-real(O_R)))+pi; Th8 = VTh8(1);
initGuesses = [20*D2R, -90*D2R];  % Th6 Th7

% Linkage Loop Eqn starting from point D
RL = NBarLinkage([r5 r6 r7 r8; Th5 NaN NaN Th8], [1,4], initGuesses, PosVectors=[1 0 1 1]);         % Make RockerLinkage
% [RLrVectors, RLthVectors] = CalcLinkageRange(SL,VTh3,fullSoltn=1);

drivingLinkageVector = zeros(2,4,length(Vr8));
drivingLinkageVector(1,:,:) = [Vr5; Vr6; Vr7; Vr8];
drivingLinkageVector(2,:,:) = [VTh5; VTh6; VTh7; VTh8];

[VTh6, VTh7] = CalcChangingLinkage(RL,drivingLinkageVector);
VTh6 = squeeze(VTh6)';
VTh7 = squeeze(VTh7)';

R6 = -r6.*exp(1i.*VTh6);
R7 = r7*exp(1i*VTh7);
R8 = D - O_R;

%Two methods for finding point F - if one fails for an unknown reason try the other
F = E - R6;%E - R6;
F2 = O_R - R7;

VTh7_2 = atan2(imag(O_R-F2),real(O_R-F2));
end

if (sum((real(F-F2).^2 + imag(F-F2).^2).^0.5 > 1) > 0 ) % If these deviate much, one is wrong
    disp("WARNING: F and F2 differ - at least one of them is incorrect")
    
    figure(6)
    R6B = F2-E;
    R7B = F2-O_R;
    plot(1:length(R5),(real(R5).^2 + imag(R5).^2).^.5);
    hold on
    plot(1:length(R5),(real(R6).^2 + imag(R6).^2).^.5);
    plot(1:length(R5),(real(R6B).^2 + imag(R6B).^2).^.5);
    plot(1:length(R5),(real(R7).^2 + imag(R7).^2).^.5);
    plot(1:length(R5),(real(R7B).^2 + imag(R7B).^2).^.5);

    title("Rocker Linkage length vs idx")
    legend("`Link5","Link6 (F->E)","Link6B (F2->E)","Link7 (F->O_R)","Link7B (F2->O_R)")
    hold off

end

%% Sanity check to ensure mechanisms arent't stretching - If it is correct, the pullrod length and rocker length should be constant
trigcalcs = false;
% plot_pullrod_length = true;
if (trigcalcs)
    % Calculate angle of rocker arm using lengths from O2-E, D-E, O2-D
    l_OrE = sqrt((real(E)-real(O_R)).^2 - (imag(E)-imag(O_R)).^2);    % distance from rocker origin to UPRIGHT pull rod attachment point
    %th_IF = acos((rocker_pull_radius.^2 - l_pullrod.^2 - l_OrF.^2)./(-2 .* l_pullrod .* l_OrF));   % angle EFO2
    th_E = atan2((imag(E)-imag(O_R)), (real(E)-real(O_R)));
    th_O = real(acos((l_pullrod.^2 - l_OrE.^2 - rocker_pull_radius.^2)./(-2*l_OrE*rocker_pull_radius)));
    th_E2 = pi/2 - (th_O + th_E);
    %th_HF = th_IF + th_F;

    %Ex = Fx - l_pullrod*cos(th_HF);
    %Ey = Fy - l_pullrod*sin(th_HF);
    Fx_trig = real(O_R) - rocker_pull_radius*cos(th_E2);
    Fy_trig = imag(O_R) + rocker_pull_radius*sin(th_E2);
end

% figure(3)
VL_pullrod = ((real(F) - real(E)).^2 + (imag(F) - imag(E)).^2).^0.5;          % Length of pullrod
VL_pullrod_rock = ((real(F) - real(O_R)).^2 + (imag(F) - imag(O_R)).^2).^0.5; % Length of pullrod-rocker link


plot_pullrod_length = true; plot_pullrod_rock_length = plot_pullrod_length;

if (max(VL_pullrod_rock) - min(VL_pullrod_rock) > 1 | max(VL_pullrod) - min(VL_pullrod) > 1)
    subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
    plot(1:length(VL_pullrod), VL_pullrod, 1:length(VL_pullrod_rock), VL_pullrod_rock)
    title("Pullrod / Rocker-Pull Link Length vs Idx")
    legend("Pullrod Len", "Pullrod-Rock Len", "location","best")
    xlabel("Index")
    grid on
end


%%
% if max(VL_pullrod) - min(VL_pullrod) > 1	% if the pullrod stretches by more than 1mm, plot length vs index
%     plot_pullrod_length = true;
% elseif 
%     plot_pullrod_rock_length = true;
%     plot_pullrod_length = false;
% else
%     plot_pullrod_length = false;
%     plot_pullrod_rock_length = false;
% end

% if (trigcalcs)
%     subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
%     VL_pullrod_trig = ((abs(Fx_trig - real(E))).^2 + (abs(Fy_trig - imag(E))).^2).^0.5;
%     plot(1:length(VL_pullrod), VL_pullrod, 1:length(VL_pullrod_trig), VL_pullrod_trig)
%     legend("kinematics","trig")
%     title("Trig vs Kinematics pullrod length")
%     grid on
% elseif plot_pullrod_length == true
%     subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
%     plot(1:length(VL_pullrod), VL_pullrod)
%     title("Pullrod Length vs Index")
%     grid on
% elseif plot_pullrod_rock_length == true
%     subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
%     plot(1:length(VL_pullrod_rock), VL_pullrod_rock)
%     title("Pullrod Rocker Length vs Index")
%     grid on
% end

%% 3D Forces

%S3D = Suspension_3D()

%% Motion Ratio Calculations
r_rocker_shockside = 80;   % mm
rocker_link_angle = 120;    % the angle between the two rocker links (ANIMIATION)
% shock_mount_pos = O_R +.75*r_rocker_shockside + 1i*(25.4*7.7);   %Shock mount pos wrt to rocker rotation axis; shock max len = 260mm
shock_mount_pos = O_R - 210 + 1i*(-20);   %Shock mount pos wrt to rocker rotation axis; shock max len = 260mm
% shock_mount_pos = O_R -.75.*r_rocker_shockside + 1i*(25.4*7.7);

% if ()

VTh9 = VTh7 + rocker_link_angle*D2R;
R9 = exp(1i*VTh9) .* r_rocker_shockside;  % Rocker Shock-Side link

% VTh9_2 = VTh7_2 - rocker_link_angle*D2R;
% R9_2 = exp(1i*VTh9_2) .* r_rocker_shockside;

G = O_R + R9;
% G2 = O_R + R9_2;
shock_length = ( real(shock_mount_pos - G).^2 + imag(shock_mount_pos - G).^2 ).^ 0.5;
spring_dx = diff(shock_length);

subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
plot(1:length(VL_pullrod),VTh7*R2D)
title("Rocker Angle vs Index")
xticks(0:10:100), yticks( (round(min(VTh7*R2D),2,"significant")-5):5:max(VTh7*R2D)+5  )
grid on
sa = min(shock_length); sb = max(shock_length);
fprintf("=== Shock Length ===\n\tmin: %.1f\n\tMax: %.1f\n\tDiff: %.1f \nover %.f mm of travel\n",sa, sb, sb-sa, max(bump)-min(bump));

%% Motion Ratio Optimization
% subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
SINGLE_MR_PLOT = true;
MR_RANGE = true;
MR_RANGE_3D = false;

if SINGLE_MR_PLOT
    figure(4)
    clf
%     shock_mount_pos = O_R -.75.*r_rocker_shockside + 1i*(25.4*7.7);
%     shock_mount_pos = O_R - 190 + 1i*(-20);
    bump_velocity = diff(bump);
    shock_length = ( real(shock_mount_pos - G).^2 + imag(shock_mount_pos - G).^2 ).^ 0.5;
    spring_dx = diff(shock_length);
    motion_ratio = abs(spring_dx) ./ abs(bump_velocity);
    plot(bump(2:end),motion_ratio)
    xticks(0-60:10:60)%, yticks(linspace( (round(min(spring_velocity),1,"significant")), max(spring_velocity) ,10 ) )
    ylim([.3,1.5])
    yticks(.3:.1:1.4)


    title("Motion Ratio vs Wheel Displacement")
    xlabel("Vertical Wheel Displacement [mm]")
    ylabel("Motion Ratio")
    grid on

end
if MR_RANGE    
    Rocker_Range(VTh7, O_R, shock_mount_pos, bump, rocker_link_angle, r_rocker_shockside);
end



%% Animation

% Determine the positions of the four 
% revolute joints at each iteration
plotlinkage = true;
plotslice = true;
silce_rainbowplot = true;

drawlabels = true;
drawVecLabels = true;
drawPtLabels = true;

if plotlinkage || plotslice
    
    %   link #:     1     2     3     4     5     6      7       8       9           10              11      12      13       14
    links = cat(3,[A;B],[B;D],[A;C],[C;D],[D;E],[F;E],[F;O_R],[WM;WC],[G;O_R],[G;shock_mount_pos],[WC;WCP],[E;F2],[F2;O_R]);%,[G2;O_R]);
    links = permute(links,[1 3 2]);     % Rearrange links to result in array: 2 x NumBars x NumIndexes
                                        % row 1 is the starting point, row 2 is end
    xlinks = real(links);
    ylinks = imag(links);
    midpts = mean(links,1); % used in vector labelling
    midpts(2,:,:) = imag(midpts);   midpts(1,:,:) = real(midpts(1,:,:));
    
    % Setup vector labelling
    labelledVecs = [6,12,7,13];%,9,14];
    vecLabels = ["R6","R6B","R7","R7B"]; %,"R9","R9B"];
    LTxtOps = ["VerticalAlignment","top","HorizontalAlignment","left"];     % Options for left text placement
    RTxtOps = ["VerticalAlignment","top","HorizontalAlignment","right"];    % Options for right text placement
    LVecs = logical([1 0 1 0]); %1:2:length(labelledVecs);     % left text
    RVecs = ~LVecs;             %2:2:length(labelledVecs);    % right text
    
    vecXL = squeeze(midpts(1,labelledVecs(LVecs),:));
    vecXR = squeeze(midpts(1,labelledVecs(RVecs),:));
    vecYL = squeeze(midpts(2,labelledVecs(LVecs),:));
    vecYR = squeeze(midpts(2,labelledVecs(RVecs),:));
    
%     text(vecXl(1,:,t), vecYl(1,:,t), vecLabels(labelledVecs(lVecs)),LTxtOps{:});
    
    % Setup Point labeling
    labelledPts = [F' F2' O_R']';
    ptLabels = ["F","F2","O_R"];
    LPts = logical([1 0 1]);
    RPts = ~LPts;
    
    xLPts = real(labelledPts(LPts,:));
    xRPts = real(labelledPts(RPts,:));
    yLPts = imag(labelledPts(LPts,:));
    yRPts = imag(labelledPts(RPts,:));
    
    
    if drawVecLabels && drawPtLabels
%         labeledStuff = [labelledVecs labelledPts];
        LLabels = [vecLabels(LVecs) ptLabels(LPts)];
        RLabels = [vecLabels(RVecs) ptLabels(RPts)];
        
        % X left & right pts
        XL = [vecXL' xLPts']';
        XR = [vecXR' xRPts']';
        % Y left & right pts
        YL = [vecYL' yLPts']';
        YR = [vecYR' yRPts']';
        
    elseif drawVecLabels
%         labeledStuff = labelledVecs;
        LLabels = LVecs;
        RLabels = RVecs;
        
        % X left & right pts
        XL = vecXL;
        XR = vecXR;
        % Y left & right pts
        YL = vecYL;
        YR = vecYR;
    elseif drawPtLabels
%         labeledStuff = labelledPts;
        LLabels = ptLabels(LPts);
        RLabels = ptLabels(RPts);
        
        % X left & right pts
        XL = xLPts;
        XR = xRPts;
        % Y left & right pts
        YL = yLPts;
        YR = yRPts;
    else
%         labeledStuff = [];
        LLabels = [];
        RLabels = [];
    end
    
    
    % Store plot handles to color-code lines via group (and maybe for speed?)
    h = cell(1,size(links,3));   % create cell array for line handles
    pullrod = 6;
    rocker_links = [7,9];
    secondary_pullrod_linkage = [11,12]; %,14];
                                
    % Find Plot Limits
    xmin = min(xlinks, [], [1 2 3]);
    xmax = max(xlinks, [], [1 2 3]);
    ymin = min(ylinks, [], [1 2 3]);
    ymax = max(ylinks, [], [1 2 3]);
    
    % Calculate plot gap space
    Space=0.03*max([abs(xmin) abs(xmax) ...
                    abs(ymin) abs(ymax)]);
    % Put equal gap around mechanism plot
    xmin = xmin-Space; 
    xmax = xmax+Space;
    ymin = ymin-Space;
    ymax = ymax+Space;

    % Mechanism Animation
%     figure(3) %2)%"WindowState","maximized") % Large Figure
    figure(3)    
    grid on
    if plotlinkage
        % TODO: ideas for later:
        %   1. update line positions instead of redrawing (if not fast enough)
        %   2. plot tubing sizes, etc (especially static stuff like chassis)
        %   3. manually ctrl draw timing (drawnow + setting to not autoupdate?)
        %   4. remove/find alt. to plotting w/ *, they slow things down
        for t=1:length(VTh3)        
            h{t} = plot(xlinks(:,:,t), ylinks(:,:,t),"b","LineWidth",2);
            set(h{t}(rocker_links),"Color",[0,.75,.5]);% = repmat([1 0 0],length(pullrod_linkage),1);    % Color the pullrod linkage
%             set(h{t}(secondary_pullrod_linkage),"Color",[.75,0,.75]);
            set(h{t}(pullrod),"Color",[1,0,0])
            
            if drawlabels
            % Label vectors at midpoint
                text(XL(:,t), YL(:,t), LLabels, LTxtOps{:});
                text(XR(:,t), YR(:,t), RLabels, RTxtOps{:});
%                 text(vecXL(1,:,t), vecYL(1,:,t), vecLabels(LVecs),LTxtOps{:});
%                 text(vecXR(1,:,t), vecYR(1,:,t), vecLabels(RVecs),RTxtOps{:});
            end
            
            axis equal   % Equal scale of x and y-axis 
            axis([xmin xmax ymin ymax]);
            M(t)=getframe; % For assembling movie frames
        end
    
    elseif plotslice
        t = round(length(links)/2);
        if silce_rainbowplot
            plot(real(links(:,:,t)), imag(links(:,:,t)),"LineWidth",2);   % ,bar5x,bar5y not plotted
        else
        h{t} = plot(xlinks(:,:,t), ylinks(:,:,t),"b*-","LineWidth",2);
        set(h{t}(rocker_links),"Color",[0,.75,.5]);% = repmat([1 0 0],length(pullrod_linkage),1);    % Color the pullrod linkage
%             set(h{t}(secondary_pullrod_linkage),"Color",[0,.5,1]);
        set(h{t}(pullrod),"Color",[1,0,0])
        end
        
        if drawlabels
        % Label vectors at midpoint
            text(midpts(1,labelledVecs(LVecs),t),midpts(2,labelledVecs(LVecs),t),...
                vecLabels(labelledVecs(LVecs)),LTxtOps{:});   % odd vectors
            text(midpts(1,labelledVecs(RVecs),t),midpts(2,labelledVecs(RVecs),t),...
                vecLabels(labelledVecs(RVecs)),RTxtOps{:});   % odd vectors
        end

        axis equal   % Equal scale of x and y-axis 
        axis([xmin xmax ymin ymax]);
    end
    %%% End Animation
    
end

%%
figure(7)
title("Rocker Linkage Angle vs idx ")
hold on
plot(1:length(R5),VTh5*R2D);
plot(1:length(R5),VTh6*R2D);
plot(1:length(R5),VTh7*R2D);
plot(1:length(R5),(VTh7-VTh9)*R2D);
legend("Th5","Th6","Th7","Th7->Th9")
hold off
% xlim([1, length(R5)*3.5])

%% Fxns
function [] = Rocker_Range(VTh7, O_R, shock_mount_pos, bump, rocker_link_angle, r_rocker_shockside)
    D2R = pi/180;
    num_changes = 5;
    
    figure(5)
    clf
    untouched_rocker_link_angle = rocker_link_angle;
    untouched_r_rocker_shockside = r_rocker_shockside;

    rock_angle = linspace(abs(rocker_link_angle) - 10, abs(rocker_link_angle) + 10, num_changes);     % The angle between the pullrod hole, the rocker axis, and the rocker shock mounting hole
    r_shock = linspace(r_rocker_shockside - 10,r_rocker_shockside + 10, num_changes);        % The radius from the rocker axis to the shock hole
    
    
%   This section iterates over an array change_vals which contains values for each optimization variable. Each "slice" 
%   in the third dimension represents a different set of changes. 
%   For example, using the array below will plot the motion ratio vs. a varying rocker angle in the 1st iter, then a
%   varying r_shock in the 2nd iter.
%       change_vals = [mean(r_shock)*ones(size(rock_angle)); rock_angle];     % constant shock radius, varying rocker angle
%       change_vals(:,:,2) = [mean(rock_angle)*ones(size(r_shock)); r_shock];

%       optimization variables: [r_rocker_shockside, rocker angle]
    change_vals =        [untouched_r_rocker_shockside*ones(size(rock_angle)); rock_angle];     % constant shock radius, varying rocker angle
    change_vals(:,:,2) = [r_shock;                              untouched_rocker_link_angle*ones(size(r_shock))];   % const. 
    
    clf % clear current figure
    changevars = 2; % number of variables to plot over
    for iter = 1:changevars
        subplot(2,1,iter)
        hold on
        for i = 1:length(change_vals)
            % Iterate variables
            r_rocker_shockside = change_vals(1,i,iter);
            rocker_link_angle = change_vals(2,i,iter);    % the angle between the two rocker links
            
            % Recalculate results
            R9 = exp(1i*(VTh7 + D2R.*rocker_link_angle)) .* r_rocker_shockside;     % Rocker Shock-Side link
            G = O_R + R9;
            shock_length = ( real(shock_mount_pos - G).^2 + imag(shock_mount_pos - G).^2 ).^ 0.5;
            spring_dx = diff(shock_length);
            bump_velocity = diff(bump);
            motion_ratio = abs(spring_dx) ./ abs(bump_velocity);
            
            % plot
            plot(bump(2:end),motion_ratio)
            grid on
        end
        xticks(0-60:10:60)%, yticks(linspace( (round(min(spring_velocity),1,"significant")), max(spring_velocity) ,10 ) )
        
        if iter == 1
            legend(string(round(change_vals(2,:,iter),3,"significant") ) + " deg")
            title("MR v. Wheel Disp. - Varying rocker angle")
        elseif iter == 2
%             change_vals = [r_shock; 110*ones(size(r_shock))];
            legend(string(round(change_vals(1,:,iter),3,"significant") ) + " mm")
            title("MR v. Wheel Disp. - Varying shock rocker radius")
%         elseif iter == 3
        end
    end
    
    rocker_link_angle = untouched_rocker_link_angle;
    r_rocker_shockside = untouched_r_rocker_shockside;
end


% function [MRvBump] = calc_MR()

