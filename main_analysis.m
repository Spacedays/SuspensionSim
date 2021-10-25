

%% Roll Center Calculations


%% Pull Rod Calculations

%% Setup
motion_ratio = 1.0

% linkage = FourBarLinkage(





%% 

%% Setup
% clc;clear;close all
D2R = pi/180;    %deg 2 radians
R2D = 180/pi;

% Suspension Linkage: Theta 3 drives the linkage, Theta 1 and all lengths are fixed.


Th1 = acos(80.83/137.5);
Th3 = 0   % Driving Input Starting Value ((var will be iterated over))
r1 = 159.56;    % mm
r2 = 304.8;
r3 = 381;
r4 = 190.5;

initGuesses = [10*D2R,90*D2R]; %Theta 2 and Theta 4 initial guesses
linkage = FourBarLinkage([r1 r2 r3 r4; Th1 NaN Th3 NaN], [2,3], initGuesses);

VTh3 = linspace(-7.5,7.5)*D2R; % VTh3 is array of Theta3 angles to iterate over
VTh1 = ones(1,length(VTh3))*Th1;
VTh2 = zeros(1,length(VTh3));
VTh4 = zeros(1,length(VTh3));


[Th2,Th4] = CalcLinkage(linkage,Th3);


%% Loop through every position

%Changing Theta 2 values
for k=1:(length(VTh3))/2
% Solving for Position
  Th3 = VTh3(k);
  [Xtemp, fval] = fsolve(@PosEq5bar,Xinit,opt);
  %Rem(k)=sqrt(fval(1)^2+fval(2)^2);  
  Th3=Xtemp(1); % Split off solutions
  Th4=Xtemp(2);
  VTh3(k)=Th3;  % Store solutions to vector
  VTh4(k)=Th4;
  % Use last solution for next guess
  Xinit=[Th3,Th4];
end  % End loop through every position
%%% End loops


%% Animation
% Determine the positions of the four 
% revolute joints at each iteration
% x Position of the origin

Ox = zeros(1,length(VTh3));
% y position of the origin
Oy = zeros(1,length(VTh3));
% x coordinate for point B
Bx = Ox + real(r1*exp(1i*Th1*ones(length(VTh3)) ));
% y coordinate for point B
By = Oy + imag(r1*exp(1i*Th1*ones(length(VTh3)) ));
% x coordinate for point C
Cx = Ox + real(r3*exp(1i*VTh3));
% y coordinate for point C
Cy = Oy + imag(r3*exp(1i*VTh3));
% x coordinate for for point D
Dx = Cx + real(r4*exp(1i*VTh2));
% y coordinate for for point D
Dy = Cy + imag(r4*exp(1i*VTh2));

% Find Plot Limits
xmin = min([min(Ox) min(Bx) min(Cx) min(Dx)]);
xmax = max([max(Ox) max(Bx) max(Cx) max(Dx)]);
ymin = min([min(Oy) min(By) min(Cy) min(Dy)]);
ymax = max([max(Oy) max(By) max(Cy) max(Dy)]);
% Calculate plot gap space
Space=0.03*max([abs(xmin) abs(xmax) ...
                abs(ymin) abs(ymax)]);
% Put equal gap around mechanism plot
xmin = xmin-Space; 
xmax = xmax+Space;
ymin = ymin-Space;
ymax = ymax+Space;

% Mechanism Animation
figure()%'WindowState','maximized') % Large Figure
for t=1:length(VTh3)
    bar1x=[Ox(t) Bx(t)]; % Coordinates Link 2
    bar1y=[Oy(t) By(t)];
    bar2x=[Bx(t) Dx(t)]; % Coordinates Link 2
    bar2y=[By(t) Dy(t)];
    bar3x=[Ox(t) Cx(t)]; % Coordinates Link 3
    bar3y=[Oy(t) Cy(t)];
    bar4x=[Cx(t) Dx(t)]; % Coordinates Link 4
    bar4y=[Cy(t) Dy(t)]; 
    plot(bar1x,bar1y,bar2x,bar2y,bar3x,bar3y,bar4x,bar4y,'LineWidth',2);   
    axis equal   % Equal scale of x and y-axis 
    axis([xmin xmax ymin ymax]);
    M(t)=getframe; % For assembling movie frames
end
%%% End Animation

