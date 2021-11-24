
clc;clear;%close all
%% Plot setup
NSubplots = 3;  % 1-Camber, 2-Scrub, 3-Pullrod Motion Ratio, 4-Pullrod length
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

%% Making some vectors
A = (216.23 + 1i*5*25.4)*ones(1,length(VTh3));
R1 = r1.*exp(1i.*VTh1.*ones(1,length(VTh3)) );  B = A+R1;
R2 = r2*ones(1,length(VTh3)).*exp(1i.*VTh2);    D = B+R2;
R3 = r3*ones(1,length(VTh3)).*exp(1i.*VTh3);    C = A+R3;
R4 = r4*ones(1,length(VTh3)).*exp(1i.*VTh4);    

%% Four Bar Linkage Coordinates

% These vars are also used for animation
% Ax = 2.5*25.4*ones(1,length(VTh3));     % x Position of the origin
% Ay = 5*25.4*ones(1,length(VTh3));       % y position of the origin
% Bx = Ax + real(r1*exp(1i*Th1*ones(length(VTh3)) ));% x coordinate for point B
% By = Ay + imag(r1*exp(1i*Th1*ones(length(VTh3)) ));% y coordinate for point B
% Cx = Ax + real(r3*exp(1i*VTh3));    % x coordinate for point C
% Cy = Ay + imag(r3*exp(1i*VTh3));    % y coordinate for point C
% Dx = Cx + real(r4*exp(1i*VTh4));    % x coordinate for for point D
% Dy = Cy + imag(r4*exp(1i*VTh4));    % y coordinate for for point D

Ax = real(A); Ay = imag(A);
% Bx = real(A+R1); By = imag(A+R1);
Bx = real(B); By = imag(B);  
% Cx = real(A+R3); Cy = imag(A+R3);
Cx = real(C); Cy = imag(C);
% Dx = A + real(R2); Dy = imag(Bx+R2);
Dx = real(D); Dy = imag(D);




%% Pull Rod Calculations
rocker_offset = [-25.4,3*25.4]; %[-15.4,40.8];   % x,y wrt lower right A-arm mount [mm] 
rocker_pull_radius = 1.5*25.4; %1.5*25.4;   	% [mm]
l_pullrod = 405;% 16.5*25.4;  
pullrod_upright_off = [0 0]; % x,y offset from upper upright suspension link (point D)

Orx = Ax + ones(1,length(VTh3))*rocker_offset(1);   % Rocker Origin x
Ory = Ay + ones(1,length(VTh3))*rocker_offset(2);   % Rocker Origin x
E = D;      % UPRIGHT Pull Rod Attachment Point (point E)
Ex = real(E);
Ey = imag(E);

% Setup drivingLinkageVector

Vr8 = ((Dx-Orx).^2+(Dy-Ory).^2).^0.5;   r8 = Vr8(1);  % Distance from rocker axis to upper upright susp. link
r5 = sum(pullrod_upright_off.^2).^0.5;  Vr5 = r5*ones(size(Vr8));   % Distance from upper upright suspension link to pullrod mount
r6 = l_pullrod; Vr6 = r6*ones(size(Vr8));
r7 = rocker_pull_radius;    Vr7 = r7*ones(size(Vr8));
Th5 = atan(pullrod_upright_off(2)/pullrod_upright_off(1));  % Angle of the vector from upper upright susp. link to pullrod mount
if (isnan(Th5)), Th5 = 0; end
VTh5 = Th5*ones(size(Vr8));
VTh6 = NaN*ones(size(Vr8));
VTh7 = VTh6;
VTh8 = atan((Dy-Ory)./(Dx-Orx))+pi; Th8 = VTh8(1);
initGuesses = [90*D2R, 120*D2R];  % Th6 Th7

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

F = E - R6;
F2 = Orx + 1i*Ory + real(R7) + imag(R7)*1i;
% F2 = 

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
%     Fx = Orx - real(Vr7.*exp(1i.*VTh7));
%     Fy = Ory - imag(Vr7.*exp(1i.*VTh7));
    Fx = real(F);
    Fy = imag(F);
end

% Ex = Orx ;         % x coordinate for point E
% Ey = Oy + rocker_origin(2);         % y coordinate for point E


%% Camber Calcs
% figure(1)
subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
camber = atand((Dx-Cx)./(Dy-Cy));
bump = r3*sin(VTh3);
plot(bump,camber)
grid on
title("Camber vs Bump")
xlabel("Bump [mm]"), ylabel("Camber [deg]")
xticks(-60:10:60), yticks(-4:1:3)

%% Scrub Calcs
wheel_dx = (1+3.5)*25.4;
wheel_r = 10*25.4;
WM = C + .5 .* R4; % Wheel Mount Location
WMx = real(WM); WMy = imag(WM);
WC = WM + wheel_dx .* exp(1i.*(VTh4 - pi/2));    % Wheel centerline
WCx = real(WC); WCy = imag(WC);
WCP = WC + wheel_r .* exp(1i.* (VTh4 + pi));
WCPx = real(WCP); WCPy = imag(WCP);

% figure(2)
subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
% plot(VTh3*R2D, (WCPx-min(WCPx)));
plot(bump, (WCPx-WCPx(51)))
% title("Wheel Scrub vs. \theta_3")
title("Wheel Scrub vs. Bump")
% xlabel("\theta_3 [deg]"), ylabel("Wheel Scrub [mm]")
xlabel("Bump [mm]"), ylabel("Wheel Scrub [mm]")
xticks(-60:10:60), yticks(-15:2.5:5)
grid on

%% Motion Ratio Calcs


% figure(3)
subplot(SPRows,SPCols,sub_idx);   sub_idx = sub_idx+1;
VL_pullrod = ((abs(Fx - Ex)).^2 + (abs(Fy - Ey)).^2).^0.5;    % Length of pullrod
if (trigcalcs && kincalcs)
    VL_pullrod_trig = ((abs(Fx_trig - Ex)).^2 + (abs(Fy_trig - Ey)).^2).^0.5;
    plot(1:length(VL_pullrod), VL_pullrod, 1:length(VL_pullrod_trig), VL_pullrod_trig)
    legend('kinematics','trig')
    title("Trig vs Kinematics pullrod length")
    grid on
else 
    plot(1:length(VL_pullrod), VL_pullrod)
    title("Pullrod Length vs Index")
    grid on
end

%% Animation
% Determine the positions of the four 
% revolute joints at each iteration
plotlinkage = true;

if plotlinkage
    
    % Find Plot Limits
    xmin = min([min(Ax) min(Bx) min(Cx) min(Dx) min(Ex) min(Ex) min(Orx) min(min(Fx)) min(WMx) min(WCx) min(WCPx)]);
    xmax = max([max(Ax) max(Bx) max(Cx) max(Dx) max(Ex) max(Ex) max(Orx) max(max(Fx)) max(WMx) max(WCx) max(WCPx)]);
    ymin = min([min(Ay) min(By) min(Cy) min(Dy) min(Ey) min(Ey) min(Ory) min(min(Fy)) min(WMy) min(WCy) min(WCPy)]);
    ymax = max([max(Ay) max(By) max(Cy) max(Dy) max(Ey) max(Ey) max(Ory) max(max(Fy)) max(WMy) max(WCy) max(WCPy)]);

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
        bar1x=real([A(t) B(t)]); % Coordinates Link 1
        bar1y=imag([A(t) B(t)]);
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
        bar7x=[Fx(t) Orx(t)]; % Coordinates Link 7 - rocker radius
        bar7y=[Fy(t) Ory(t)];
        wheelax_x = [WMx(t) WCx(t)];    wheelax_y = [WMy(t) WCy(t)];
        wheelr_x = [WCx(t) WCPx(t)];   wheelr_y = [WCy(t) WCPy(t)];
        plot(bar1x,bar1y,bar2x,bar2y,bar3x,bar3y,bar4x,bar4y,bar6x,bar6y,bar7x,bar7y,...
            wheelax_x,wheelax_y, wheelr_x, wheelr_y,'LineWidth',2);   % ,bar5x,bar5y not plotted
        axis equal   % Equal scale of x and y-axis 
        axis([xmin xmax ymin ymax]);
        M(t)=getframe; % For assembling movie frames
    end
    %%% End Animation
end
