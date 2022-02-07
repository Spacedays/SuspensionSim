
clc;clear;%close all
%% Plot setup
% 1-Camber, 2-Scrub, 3-Pullrod Motion Ratio, 4-Pullrod length
SPRows = 2; SPCols = 2;
sub_idx = 1;
figure(1); 

%% Roll Center Calculations
% TBD


%% Kinematics Setup

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
VTh3 = linspace(-10,10,101)*D2R; % VTh3 is array of Theta3 angles to iterate over

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
% figure(1)
subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
camber = atand((real(D)-real(C))./(imag(D)-imag(C)));
bump = r3*sin(VTh3);
plot(bump,camber)
grid on
title("Camber vs Wheel Displacement")
xlabel("Wheel Displacement [mm]"), ylabel("Camber [deg]")
xticks(-60:10:60), yticks(-4:1:3)

%% Scrub Calcs
wheel_dx = (1+3.5)*25.4;
wheel_r = 10*25.4;
WM = C + .5 .* R4; % Wheel Mount Location
WC = WM + wheel_dx .* exp(1i.*(VTh4 - pi/2));    % Wheel centerline
WCP = WC + wheel_r .* exp(1i.* (VTh4 + pi));

% figure(2)
subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
% plot(VTh3*R2D, (real(WCP)-min(real(WCP))));
plot(bump, (real(WCP) - real(WCP(51)) ) )
% title("Wheel Scrub vs. \theta_3")
title("Wheel Scrub vs. Wheel Displacement")
% xlabel("\theta_3 [deg]"), ylabel("Wheel Scrub [mm]")
xlabel("Wheel Displacement [mm]"), ylabel("Wheel Scrub [mm]")
xticks(-60:10:60), yticks(-15:2.5:5)
grid on


%% Pull Rod Calculations
rocker_axis_pos = [-20,50]; %[-15.4,40.8];   % x,y wrt lower right A-arm mount [mm] 
O_R = A + ones(1,length(VTh3))*(rocker_axis_pos(1) + 1i*rocker_axis_pos(2));    % Rocker axis of rotation

rocker_pull_radius = 35; %1.5*25.4;   	% [mm]
l_pullrod = 390;% 16.5*25.4;  
pullrod_upright_off = [0 0];    % x,y offset from upper upright suspension link (point D)
pullrod_Aarm_r = 25;            % Radial offset from upper A-arm mount
pullrod_Aarm_vert = 0;

% E = D + pullrod_upright_off(1) + imag(pullrod_upright_off(2);      % Pull Rod Attachment Point (point E) based on upright

E = D - pullrod_Aarm_r * exp(1i*VTh2) + pullrod_Aarm_vert * exp(1i*(VTh2+pi/4));

% Setup drivingLinkageVector

Vr8 = ((real(D)-real(O_R)).^2+(imag(D)-imag(O_R)).^2).^0.5;   r8 = Vr8(1);  % Distance from rocker axis to upper upright susp. link
r5 = sum(pullrod_upright_off.^2).^0.5;  Vr5 = r5*ones(size(Vr8));   % Distance from upper upright suspension link to pullrod mount
r6 = l_pullrod; Vr6 = r6*ones(size(Vr8));
r7 = rocker_pull_radius;    Vr7 = r7*ones(size(Vr8));
Th5 = atan(pullrod_upright_off(2)/pullrod_upright_off(1));  % Angle of the vector from upper upright susp. link to pullrod mount
if (isnan(Th5)), Th5 = 0; end
VTh5 = Th5*ones(size(Vr8));
VTh6 = NaN*ones(size(Vr8));
VTh7 = VTh6;
VTh8 = atan((imag(D)-imag(O_R))./(real(D)-real(O_R)))+pi; Th8 = VTh8(1);
initGuesses = [110*D2R, 100*D2R];  % Th6 Th7

% Linkage Loop Eqn starting from point D
RL = NBarLinkage([r5 r6 r7 r8; Th5 NaN NaN Th8], [1,4], initGuesses, NegVectors=[1 0 0 0]);         % Make RockerLinkage
% [RLrVectors, RLthVectors] = CalcLinkageRange(SL,VTh3,fullSoltn=1);

drivingLinkageVector = zeros(2,4,length(Vr8));
drivingLinkageVector(1,:,:) = [Vr5; Vr6; Vr7; Vr8];
drivingLinkageVector(2,:,:) = [VTh5; VTh6; VTh7; VTh8];

[VTh6, VTh7] = CalcChangingLinkage(RL,drivingLinkageVector);
VTh6 = squeeze(VTh6)';
VTh7 = squeeze(VTh7)';

R6 = r6.*exp(1i.*VTh6);
R7 = r7*exp(1i*VTh7);

%Two methods for finding point F - if one fails for an unknown reason try the other
F2 = E - R6;
F = O_R + R7;

if (length((real(F-F2).^2 + imag(F-F2).^2).^0.5 > 1) > 0 ) % If these deviate much, one is wrong
    disp("WARNING: F and F2 differ - at least one of them is incorrect")
end

%% Sanity check to ensure mechanisms arent't stretching - If it is correct, the pullrod length and rocker length should be constant
trigcalcs = false;
% plot_pullrod_length = true;
if (trigcalcs)
    % Calculate angle of rocker arm using lengths from O2-E, D-E, O2-D
    l_OrE = sqrt((real(E)-real(O_R)).^2 - (imag(E)-imag(O_R)).^2);    % distance from rocker origin to UPRIGHT pull rod attachment point
    %th_IF = acos((rocker_pull_radius.^2 - l_pullrod.^2 - l_OrF.^2)./(-2 .* l_pullrod .* l_OrF));   % angle EFO2
    th_E = atan((imag(E)-imag(O_R))./(real(E)-real(O_R)));
    th_O = real(acos((l_pullrod.^2 - l_OrE.^2 - rocker_pull_radius.^2)./(-2*l_OrE*rocker_pull_radius)));
    th_E2 = pi/2 - (th_O + th_E);
    %th_HF = th_IF + th_F;

    %Ex = Fx - l_pullrod*cos(th_HF);
    %Ey = Fy - l_pullrod*sin(th_HF);
    Fx_trig = real(O_R) - rocker_pull_radius*cos(th_E2);
    Fy_trig = imag(O_R) + rocker_pull_radius*sin(th_E2);
end

% figure(3)
VL_pullrod = ((abs(real(F) - real(E))).^2 + (abs(imag(F) - imag(E))).^2).^0.5;    % Length of pullrod
VL_pullrod_rock = ((abs(real(F) - real(O_R))).^2 + (abs(imag(F) - imag(O_R))).^2).^0.5;    % Length of pullrod-rocker link


plot_pullrod_length = true; plot_pullrod_rock_length = plot_pullrod_length;
if max(VL_pullrod) - min(VL_pullrod) > 1    % if the pullrod stretches by more than 1mm, plot length vs index
    plot_pullrod_length = true;
elseif max(VL_pullrod_rock) - min(VL_pullrod_rock) > 1
    plot_pullrod_rock_length = true;
    plot_pullrod_length = false;
else
    plot_pullrod_length = false;
    plot_pullrod_rock_length = false;
end

if (trigcalcs)
    subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
    VL_pullrod_trig = ((abs(Fx_trig - real(E))).^2 + (abs(Fy_trig - imag(E))).^2).^0.5;
    plot(1:length(VL_pullrod), VL_pullrod, 1:length(VL_pullrod_trig), VL_pullrod_trig)
    legend('kinematics','trig')
    title("Trig vs Kinematics pullrod length")
    grid on
elseif plot_pullrod_length == true
    subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
    plot(1:length(VL_pullrod), VL_pullrod)
    title("Pullrod Length vs Index")
    grid on
elseif plot_pullrod_rock_length == true
    subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
    plot(1:length(VL_pullrod_rock), VL_pullrod_rock)
    title("Pullrod Rocker Length vs Index")
    grid on
end

%% 3D Forces

S3D = Suspension_3D()

%% Motion Ratio Calculations
r_rocker_shockside = 128;
rocker_link_angle = 135;    % the angle between the two rocker links (ANIMIATION)
shock_mount_pos = O_R -.75*r_rocker_shockside + 1i*(25.4*7.7);   %Shock mount pos wrt to rocker rotation axis; shock max len = 260mm
R9 = exp(1i*(VTh7 - D2R*rocker_link_angle)) * r_rocker_shockside;  % Rocker Shock-Side link
G = O_R + R9;
shock_length = ( real(shock_mount_pos - G).^2 + imag(shock_mount_pos - G).^2 ).^ 0.5;
spring_dx = diff(shock_length);

subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
plot(1:length(VL_pullrod),VTh7*R2D)
title("Rocker Angle vs Index")
xticks(0:10:100), yticks( (round(min(VTh7*R2D),2,'significant')-5):5:max(VTh7*R2D)+5  )
grid on

%% Animation

% Determine the positions of the four 
% revolute joints at each iteration
plotlinkage = false;
plotslice = true;

if plotlinkage || plotslice
    links = cat(3,[A;B],[B;D],[A;C],[C;D],[D;E],[E;F],[F;O_R],[WM;WC],[WC;WCP],[G;O_R],[G;shock_mount_pos]);
    links = permute(links,[1 3 2]);     % Rearrange links to result in array: 2 x NumBars x NumIndexes
                                        % row 1 is the starting point, row 2 is end
    
    % Find Plot Limits
    xmin = min(real(links(:,:,:)), [], [1 2 3]);
    xmax = max(real(links(:,:,:)), [], [1 2 3]);
    ymin = min(imag(links(:,:,:)), [], [1 2 3]);
    ymax = max(imag(links(:,:,:)), [], [1 2 3]);
    
    % Calculate plot gap space
    Space=0.03*max([abs(xmin) abs(xmax) ...
                    abs(ymin) abs(ymax)]);
    % Put equal gap around mechanism plot
    xmin = xmin-Space; 
    xmax = xmax+Space;
    ymin = ymin-Space;
    ymax = ymax+Space;

    % Mechanism Animation
    figure(2) %2)%'WindowState','maximized') % Large Figure
    grid on
    
    if plotlinkage
        for t=1:length(VTh3)        
            plot(real(links(:,:,t)), imag(links(:,:,t)),'LineWidth',2);   % ,bar5x,bar5y not plotted
            axis equal   % Equal scale of x and y-axis 
            axis([xmin xmax ymin ymax]);
            M(t)=getframe; % For assembling movie frames
        end
    
    elseif plotslice
        t = 51;
        plot(real(links(:,:,t)), imag(links(:,:,t)),'LineWidth',2);   % ,bar5x,bar5y not plotted
            axis equal   % Equal scale of x and y-axis 
            axis([xmin xmax ymin ymax]);
    end
    %%% End Animation
end


%% Motion Ratio Optimization
% subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
SINGLE_MR_PLOT = true;
MR_RANGE = true;
MR_RANGE_3D = false;


if SINGLE_MR_PLOT
    figure(3)
    clf
    shock_mount_pos = O_R -.75.*r_rocker_shockside + 1i*(25.4*7.7);
    bump_velocity = diff(bump);
    shock_length = ( real(shock_mount_pos - G).^2 + imag(shock_mount_pos - G).^2 ).^ 0.5;
    spring_dx = diff(shock_length);
    motion_ratio = spring_dx ./ bump_velocity;
    plot(bump(2:end),motion_ratio)
    xticks(0-60:10:60)%, yticks(linspace( (round(min(spring_velocity),1,'significant')), max(spring_velocity) ,10 ) )
    ylim([.3,1.5])
    yticks(.3:.1:1.5)


    title("Motion Ratio vs Wheel Displacement")
    xlabel("Vertical Wheel Displacement [mm]")
    ylabel("Motion Ratio")
    grid on
end

if MR_RANGE
    num_changes = 5
    
    figure(4)
    clf
    untouched_rocker_link_angle = rocker_link_angle;
    untouched_r_rocker_shockside = r_rocker_shockside;

    rock_angle = linspace(rocker_link_angle - 10, rocker_link_angle + 7, num_changes);     % The angle between the pullrod hole, the rocker axis, and the rocker shock mounting hole
    r_shock = linspace(r_rocker_shockside - 10,r_rocker_shockside + 10, num_changes);        % The radius from the rocker axis to the shock hole
    
    
%   This section iterates over an array change_vals which contains values for each optimization variable. Each 'slice' 
%   in the third dimension represents a different set of changes. 
%   For example, using the array below will plot the motion ratio vs. a varying rocker angle in the 1st iter, then a
%   varying r_shock in the 2nd iter.
%       change_vals = [mean(r_shock)*ones(size(rock_angle)); rock_angle];     % constant shock radius, varying rocker angle
%       change_vals(:,:,2) = [mean(rock_angle)*ones(size(r_shock)); r_shock];

%       optimization variables: [r_rocker_shockside, rocker angle]
    change_vals =        [untouched_r_rocker_shockside*ones(size(rock_angle)); rock_angle];     % constant shock radius, varying rocker angle
    change_vals(:,:,2) = [r_shock;                              untouched_rocker_link_angle*ones(size(r_shock))];
    
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
            shock_mount_pos = O_R -.75.*r_rocker_shockside + 1i*(25.4*7.7);   %Shock mount pos wrt to rocker rotation axis; shock max len = 260mm
            R9 = exp(1i*(VTh7 - D2R.*rocker_link_angle)) .* r_rocker_shockside;  % Rocker Shock-Side link
            G = O_R + R9;
            shock_length = ( real(shock_mount_pos - G).^2 + imag(shock_mount_pos - G).^2 ).^ 0.5;
            spring_dx = diff(shock_length);
            bump_velocity = diff(bump);
            motion_ratio = spring_dx ./ bump_velocity;
            
            % plot
            plot(bump(2:end),motion_ratio)
            grid on
        end
        xticks(0-60:10:60)%, yticks(linspace( (round(min(spring_velocity),1,'significant')), max(spring_velocity) ,10 ) )
        
        if iter == 1
            legend(string(round(change_vals(2,:,iter),3,'significant') ) + " deg")
            title("MR v. Wheel Disp. - Varying rocker angle")
        elseif iter == 2
%             change_vals = [r_shock; 110*ones(size(r_shock))];
            legend(string(round(change_vals(1,:,iter),3,'significant') ) + " mm")
            title("MR v. Wheel Disp. - Varying shock rocker radius")
%         elseif iter == 3
        end
    end
    
    
% This is an attempt at 3d optimization, but it is flawed 
elseif MR_RANGE_3D
    surf_anim = true;
    
    r_shock = 15:15:150;
    rock_angle = 90:10:160;
    [X,Y] = meshgrid(rock_angle,r_shock);
    ZZ = zeros(length(r_shock),length(rock_angle),length(bump)-1);
    
%     change_vals = [75*ones(size(a)); a];  % rocker link angle and rocker shockside radius

    result_vals = zeros(length(r_shock),length(rock_angle));
    for i = 1:length(r_shock)
        for j = 1:length(rock_angle)
            r_rocker_shockside = r_shock(i);
            rocker_link_angle = rock_angle(j);    % the angle between the two rocker links
            shock_mount_pos = O_R -.75.*r_rocker_shockside + 1i*(25.4*7.7);   %Shock mount pos wrt to rocker rotation axis; shock max len = 260mm
            R9 = exp(1i*(VTh7 - D2R.*rocker_link_angle)) .* r_rocker_shockside;  % Rocker Shock-Side link
            G = O_R + R9;
            spring_dx = diff(( real(shock_mount_pos - G).^2 + imag(shock_mount_pos - G).^2 ).^ 0.5);  % spring length vs index
            bump_velocity = diff(bump);
            motion_ratio = spring_dx ./ bump_velocity;
            ZZ(i,j,:) = motion_ratio;
        end
    end
    figure(5)
    clf
    surf(X,Y,ZZ(:,:,51))
    xlabel("Rocker Link Angle"), ylabel("Shock Side Length"), zlabel("Motion Ratio")
%     if surf_anim  
%     else
%     end

    title("Motion Ratio vs Wheel Displacement")
    grid on
end


% re-set vars overwritten by optimization (INCOMPLETE) (Animation runs after the fact, so vars used by it must be reset)
rocker_link_angle = 90;    % the angle between the two rocker links 
shock_mount_pos = O_R -.75*r_rocker_shockside + 1i*(25.4*7.7);   %Shock mount pos wrt to rocker rotation axis; shock max len = 260mm
R9 = exp(1i*(VTh7 - D2R*rocker_link_angle)) * r_rocker_shockside;  % Rocker Shock-Side link

