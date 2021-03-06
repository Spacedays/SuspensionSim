
clc;clear;close all
%% Roll Center Calculations
% TBD


%% Kinematics Setup
% Next goal: do stuff to compare different suspension setups
% clc;clear;close all
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
linkage = NBarLinkage([r1 r2 r3 r4; Th1 NaN Th3 NaN], [2,3], initGuesses);

VTh3 = linspace(-10,10)*D2R; % VTh3 is array of Theta3 angles to iterate over



%% Kinematics setup options -- Use only one of the three options

%%% Kinematics using CalcLinkageRange
% [rVectors, thVectors] = CalcLinkageRange(linkage,VTh3,fullSoltn=1);
% VTh1 = thVectors(1,:);
% VTh2 = thVectors(2,:);
% % VTh3 = thVectors(3,:);    % Th3 is already set so no need to re-set it, but it remains the same so it doesn't matter either way
% VTh4 = thVectors(4,:);


%%% Kinematics using CalcLinkage only returning the two unknowns
[VTh2,VTh4] = CalcLinkageRange(linkage,VTh3); % optional argument: fullSoltn=0
VTh1 = ones(1,length(VTh3))*Th1;
%  = solVectors(1,:);
%  = solVectors(2,:);

%%% Kinematics using CalcLinkage
% VTh1 = ones(1,length(VTh3))*Th1;
% VTh2 = zeros(1,length(VTh3));
% VTh4 = zeros(1,length(VTh3));
% 
% %Iterate over driving variable
% for k=1:(length(VTh3))
% % Solving for Position
%   Th3 = VTh3(k);
%   [Th2,Th4] = CalcLinkage(linkage,Th3);
%   VTh2(k)=Th2;                  % Store solutions to vector
%   VTh4(k)=Th4;
%   obj.priorGuesses(1) = Th2;    % Use last solution for next guess
%   obj.priorGuesses(2) = Th4;
% end
% %%% End loops

%% Pull Rod Calculations
motion_ratio = 1.0;
rocker_origin = [-15.4,40.8];   % x,y wrt lower right A-arm mount [mm] 
rocker_pull_radius = 1.5*25.4;   	% [mm]
l_pullrod = 17*25.4;  


% These vars are also used for animation
Ox = 2.5*25.4*ones(1,length(VTh3));     % x Position of the origin
Oy = 5*25.4*ones(1,length(VTh3));     % y position of the origin
Bx = Ox + real(r1*exp(1i*Th1*ones(length(VTh3)) ));% x coordinate for point B
By = Oy + imag(r1*exp(1i*Th1*ones(length(VTh3)) ));% y coordinate for point B
Cx = Ox + real(r3*exp(1i*VTh3));    % x coordinate for point C
Cy = Oy + imag(r3*exp(1i*VTh3));    % y coordinate for point C
Dx = Cx + real(r4*exp(1i*VTh4));    % x coordinate for for point D
Dy = Cy + imag(r4*exp(1i*VTh4));    % y coordinate for for point D


Orx = Ox + ones(1,length(VTh3))*rocker_origin(1);
Ory = Oy + ones(1,length(VTh3))*rocker_origin(2);
Fx = Dx;    % x coordinate for for point Pull Rod Attachment Point (point F)
Fy = Dy;    % y coordinate for for point F

% Calculate angle of rocker arm using lengths from O2-E, D-E, O2-D
l_OrF = sqrt((Fx-Orx).^2 - (Fy-Ory).^2);    % distance from rocker origin to pull rod attachment point
%th_IF = acos((rocker_pull_radius.^2 - l_pullrod.^2 - l_OrF.^2)./(-2 .* l_pullrod .* l_OrF));   % angle EFO2
th_F = atan((Fy-Ory)./(Fx-Orx));
th_O = real(acos((l_pullrod.^2 - l_OrF.^2 - rocker_pull_radius.^2)./(-2*l_OrF*rocker_pull_radius)));
th_E = pi/2 - (th_O + th_F);
%th_HF = th_IF + th_F;

%Ex = Fx - l_pullrod*cos(th_HF);
%Ey = Fy - l_pullrod*sin(th_HF);
Ex = Orx - rocker_pull_radius*cos(th_E);
Ey = Ory + rocker_pull_radius*sin(th_E);


% Ex = Orx ;         % x coordinate for point E
% Ey = Oy + rocker_origin(2);         % y coordinate for point E


%% Camber Calcs
figure(1)
grid on
title("Camber vs Bump")
camber = atand((Dx-Cx)./(Dy-Cy));
bump = r3*sin(VTh3);
plot(bump,camber)

%% Motion Ratio Calcs


%% Animation
% Determine the positions of the four 
% revolute joints at each iteration

% Find Plot Limits
% xmin = min([min(Ox) min(Bx) min(Cx) min(Dx)]);
% xmax = max([max(Ox) max(Bx) max(Cx) max(Dx)]);
% ymin = min([min(Oy) min(By) min(Cy) min(Dy)]);
% ymax = max([max(Oy) max(By) max(Cy) max(Dy)]);

xmin = min([min(Ox) min(Bx) min(Cx) min(Dx) min(Ex) min(Fx) min(Orx)]);
xmax = max([max(Ox) max(Bx) max(Cx) max(Dx) max(Ex) max(Fx) max(Orx)]);
ymin = min([min(Oy) min(By) min(Cy) min(Dy) min(Ey) min(Fy) min(Ory)]);
ymax = max([max(Oy) max(By) max(Cy) max(Dy) max(Ey) max(Fy) max(Ory)]);

% Calculate plot gap space
Space=0.03*max([abs(xmin) abs(xmax) ...
                abs(ymin) abs(ymax)]);
% Put equal gap around mechanism plot
xmin = xmin-Space; 
xmax = xmax+Space;
ymin = ymin-Space;
ymax = ymax+Space;

% Mechanism Animation
figure(2)%'WindowState','maximized') % Large Figure
grid on
for t=1:length(VTh3)
    bar1x=[Ox(t) Bx(t)]; % Coordinates Link 1
    bar1y=[Oy(t) By(t)];
    bar2x=[Bx(t) Dx(t)]; % Coordinates Link 2
    bar2y=[By(t) Dy(t)];
    bar3x=[Ox(t) Cx(t)]; % Coordinates Link 3
    bar3y=[Oy(t) Cy(t)];
    bar4x=[Cx(t) Dx(t)]; % Coordinates Link 4
    bar4y=[Cy(t) Dy(t)]; 
    bar5x=[Ex(t) Fx(t)]; % Coordinates Link 5 - pull rod
    bar5y=[Ey(t) Fy(t)]; 
    bar6x=[Orx(t) Ex(t)]; % Coordinates Link 6 - rocker radius
    bar6y=[Ory(t) Ey(t)]; 
    plot(bar1x,bar1y,bar2x,bar2y,bar3x,bar3y,bar4x,bar4y,bar5x,bar5y,bar6x,bar6y,'LineWidth',2);   
    axis equal   % Equal scale of x and y-axis 
    axis([xmin xmax ymin ymax]);
    M(t)=getframe; % For assembling movie frames
end
%%% End Animation

