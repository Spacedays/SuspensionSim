
clc;clear;%close all
%% Roll Center Calculations
% TBD


%% Kinematics Setup

D2R = pi/180;    %deg 2 radians
R2D = 180/pi;

% Suspension Linkage: Theta 3 drives the linkage, Theta 1 and all lengths are fixed.

Th1 = 59.5636*D2R;% acos(80.83/137.5);
Th3 = 0;   % Driving Input Starting Value ((var will be iterated over))
r1 = 6.282*25.4;%159.56;    % mm
r2 = 12*25.4;%304.8;
r3 = 15*25.4;%381;
r4 = 7.5*25.4;%190.5;

initGuesses = [10*D2R,90*D2R]; %Theta 2 and Theta 4 initial guesses
SL = NBarLinkage([r1 r2 r3 r4; Th1 NaN Th3 NaN], [2,3], initGuesses);   % Make SuspensionLinkage
VTh3 = linspace(-10,10)*D2R; % VTh3 is array of Theta3 angles to iterate over

% Solve Kinematics using CalcLinkageRange
[~, thVectors] = CalcLinkageRange(SL,VTh3,fullSoltn=1);    % rVectors not used here
VTh1 = thVectors(1,:);
VTh2 = thVectors(2,:);
% VTh3 = thVectors(3,:);    % Th3 is already set so no need to re-set it, but it remains the same so it doesn't matter either way
VTh4 = thVectors(4,:);

%% Four Bar Linkage Coordinates

% These vars are also used for animation
Ax = 2.5*25.4*ones(1,length(VTh3));     % x Position of the origin
Ay = 5*25.4*ones(1,length(VTh3));       % y position of the origin
Bx = Ax + real(r1*exp(1i*Th1*ones(length(VTh3)) ));% x coordinate for point B
By = Ay + imag(r1*exp(1i*Th1*ones(length(VTh3)) ));% y coordinate for point B
Cx = Ax + real(r3*exp(1i*VTh3));    % x coordinate for point C
Cy = Ay + imag(r3*exp(1i*VTh3));    % y coordinate for point C
Dx = Cx + real(r4*exp(1i*VTh4));    % x coordinate for for point D
Dy = Cy + imag(r4*exp(1i*VTh4));    % y coordinate for for point D

%% Pull Rod Calculations
rocker_offset = [-30,40]; %[-15.4,40.8];   % x,y wrt lower right A-arm mount [mm] 
rocker_pull_radius = 1.5*25.4;   	% [mm]
l_pullrod = 17*25.4;  
pullrod_upright_off = [0 0]; % x,y offset from upper upright suspension link (point D)

Orx = Ax + ones(1,length(VTh3))*rocker_offset(1);   % Rocker Origin x
Ory = Ay + ones(1,length(VTh3))*rocker_offset(2);   % Rocker Origin x
Ex = Dx;    % x coordinate for UPRIGHT Pull Rod Attachment Point (point E)
Ey = Dy;    % y coordinate for point E

% Setup drivingLinkageVector

Vr8 = ((Dx-Orx).^2+(Dy-Ory).^2).^0.5;                % Distance from rocker axis to upper upright susp. link
r8 = Vr8(1);

r5 = sum(pullrod_upright_off.^2).^0.5;           % Distance from upper upright suspension link to pullrod mount
Vr5 = r5*ones(size(Vr8));

r6 = l_pullrod;
Vr6 = r6*ones(size(Vr8));

r7 = rocker_pull_radius;
Vr7 = r7*ones(size(Vr8));

Th5 = atan(pullrod_upright_off(2)/pullrod_upright_off(1));  % Angle of the vector from upper upright susp. link to pullrod mount
if (isnan(Th5)), Th5 = 0; end
VTh5 = Th5*ones(size(Vr8));

VTh6 = NaN*ones(size(Vr8));
VTh7 = VTh6;

VTh8 = atan((Dy-Ory)./(Dx-Orx))+pi;
Th8 = VTh8(1);

initGuesses = [30*pi/180, 270*pi/180];

RL = NBarLinkage([r5 r6 r7 r8; Th5 NaN NaN Th8], [1,4], initGuesses, NegVectors=[0 -1 0 0]);         % Make RockerLinkage
% [RLrVectors, RLthVectors] = CalcLinkageRange(SL,VTh3,fullSoltn=1);

drivingLinkageVector = zeros(2,4,length(Vr8));
% drivingLinkageVector(1,1,:) = r5*ones(size(Vr8));
drivingLinkageVector(1,:,:) = [Vr5; Vr6; Vr7; Vr8];
drivingLinkageVector(2,:,:) = [VTh5; VTh6; VTh7; VTh8];

[VTh6, VTh7] = CalcChangingLinkage(RL,drivingLinkageVector);
VTh6 = squeeze(VTh6)';
VTh7 = squeeze(VTh7)';

%% Plotting and Pull Rod Calculations -- Pull Rod calcs TBD
motion_ratio = 1.0;

trigcalcs = true;
kincalcs = true;
if (trigcalcs)
    % Calculate angle of rocker arm using lengths from O2-E, D-E, O2-D
    l_OrE = sqrt((Ex-Orx).^2 - (Ey-Ory).^2);    % distance from rocker origin to UPRIGHT pull rod attachment point
    %th_IF = acos((rocker_pull_radius.^2 - l_pullrod.^2 - l_OrF.^2)./(-2 .* l_pullrod .* l_OrF));   % angle EFO2
    th_E = atan((Ey-Ory)./(Ex-Orx));
    th_O = real(acos((l_pullrod.^2 - l_OrE.^2 - rocker_pull_radius.^2)./(-2*l_OrE*rocker_pull_radius)));
    th_E2 = pi/2 - (th_O + th_E);
    %th_HF = th_IF + th_F;

    %Ex = Fx - l_pullrod*cos(th_HF);
    %Ey = Fy - l_pullrod*sin(th_HF);
    Fx_trig = Orx - rocker_pull_radius*cos(th_E2);
    Fy_trig = Ory + rocker_pull_radius*sin(th_E2);
end
if kincalcs
    % Orx, Ory for rocker rotation axis
    Fx = Orx + real(Vr7.*exp(1i.*VTh7));
    Fy = Ory + imag(Vr7.*exp(1i.*VTh7));
end

% Ex = Orx ;         % x coordinate for point E
% Ey = Oy + rocker_origin(2);         % y coordinate for point E


%% Camber Calcs
figure(1)
camber = atand((Dx-Cx)./(Dy-Cy));
bump = r3*sin(VTh3);
plot(bump,camber)
grid on
title("Camber vs Bump")

%% Motion Ratio Calcs

VL_pullrod = ((abs(Ex)-abs(Fx)).^2 + (abs(Ey)-abs(Fy)).^2).^0.5;
figure(3)
plot(1:length(VL_pullrod), VL_pullrod)
title("Pullrod Length vs Index")

if (trigcalcs && kincalcs)
    figure(4)
    VL_pullrod_trig = ((abs(Ex)-abs(Fx_trig)).^2 + (abs(Ey)-abs(Fy_trig)).^2).^0.5;
    plot(1:length(VL_pullrod), VL_pullrod, 1:length(VL_pullrod_trig), VL_pullrod_trig)
    legend('kinematics','trig')
    title("Trig vs Kinematics pullrod length")
end

%% Animation
% Determine the positions of the four 
% revolute joints at each iteration

% Find Plot Limits
% xmin = min([min(Ox) min(Bx) min(Cx) min(Dx)]);
% xmax = max([max(Ox) max(Bx) max(Cx) max(Dx)]);
% ymin = min([min(Oy) min(By) min(Cy) min(Dy)]);
% ymax = max([max(Oy) max(By) max(Cy) max(Dy)]);

xmin = min([min(Ax) min(Bx) min(Cx) min(Dx) min(Ex) min(Ex) min(Orx) min(min(Fx))]);
xmax = max([max(Ax) max(Bx) max(Cx) max(Dx) max(Ex) max(Ex) max(Orx) max(max(Fx))]);
ymin = min([min(Ay) min(By) min(Cy) min(Dy) min(Ey) min(Ey) min(Ory) min(min(Fy))]);
ymax = max([max(Ay) max(By) max(Cy) max(Dy) max(Ey) max(Ey) max(Ory) max(max(Fy))]);

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
for t=1:length(VTh3)
    bar1x=[Ax(t) Bx(t)]; % Coordinates Link 1
    bar1y=[Ay(t) By(t)];
    bar2x=[Bx(t) Dx(t)]; % Coordinates Link 2
    bar2y=[By(t) Dy(t)];
    bar3x=[Ax(t) Cx(t)]; % Coordinates Link 3
    bar3y=[Ay(t) Cy(t)];
    bar4x=[Cx(t) Dx(t)]; % Coordinates Link 4
    bar4y=[Cy(t) Dy(t)]; 
    bar5x=[Dx(t) Ex(t)]; % Coordinates Link 5 - Upright Pullrod Offset
    bar5y=[Dy(t) Ey(t)]; 
    bar6x=[Ex(t) Fx(t)]; % Coordinates Link 6 - pull rod
    bar6y=[Ey(t) Fy(t)]; 
    bar7x=[Orx(t) Fx(t)]; % Coordinates Link 7 - rocker radius
    bar7y=[Ory(t) Fy(t)];
    plot(bar1x,bar1y,bar2x,bar2y,bar3x,bar3y,bar4x,bar4y,bar6x,bar6y,bar7x,bar7y,'LineWidth',2);   % ,bar5x,bar5y not plotted
    axis equal   % Equal scale of x and y-axis 
    axis([xmin xmax ymin ymax]);
    M(t)=getframe; % For assembling movie frames
end
%%% End Animation

