

%% Roll Center Calculations


%% Pull Rod Calculations

%% Setup
motion_ratio = 1.0

% linkage = FourBarLinkage(





%% 

%% Setup
% clc;clear;close all
D2R = pi/180;    %deg 2 radians
R2D = 180/p;

% Suspension Linkage: Theta 3 drives the linkage, Theta 1 and all lengths are fixed.


Th1 = acos(80.83/137.5);
Th3 = 0   % Driving Input Starting Value ((var will be iterated over))
r1 = 159.56;    % mm
r2 = 304.8;
r3 = 381;
r4 = 190.5;

initGuesses = [10*D2R,90*D2R]; %Theta 2 and Theta 4 initial guesses
linkage = FourBarLinkage([r1 r2 r3 r4 r5; Th1 NaN Th3 NaN], [2,3], initGuesses);

VTh3 = linspace(-7.5,7.5)*D2R; % VTh3 is array of Theta3 angles to iterate over
VTh1 = ones(1,length(VTh3))*Th1;
VTh2 = zeros(1,length(VTh3));
VTh4 = zeros(1,length(VTh3));


[Th2,Th4] = calcLinkage(linkage);


% Initial guesses for Theta3 and Theta4 
% when Theta2=0 
Th3 = 90*D2R;  Th4 = 270*D2R; 
Xinit=[Th3 Th4];  % Vector of guesses
VTh3 = zeros(1,length(VTh3)); % Predefine
VTh4 = zeros(1,length(VTh3)); % vectors
opt = optimset('Display', 'off'); 
%%% End Setup

%% Loop through every position

%Changing Theta 2 values
for k=1:(length(VTh3))/2
% Solving for Position
  Th3 = VTh3(k);
  [Xtemp, fval] = fsolve(@PosEq5bar,Xinit,opt);
  Rem(k)=sqrt(fval(1)^2+fval(2)^2);  
  Th3=Xtemp(1); % Split off solutions
  Th4=Xtemp(2);
  VTh3(k)=Th3;  % Store solutions to vector
  VTh4(k)=Th4;
  % Use last solution for next guess
  Xinit=[Th3,Th4];
end  % End loop through every position


%Changing Theta 5 Values
for k=(length(VTh5)/2):length(VTh5)
  Th5 = VTh5(k);
  [Xtemp, fval] = fsolve(@PosEq5bar,Xinit,opt);
  Rem(k)=sqrt(fval(1)^2+fval(2)^2);  
  Th3=Xtemp(1); % Split off solutions
  Th4=Xtemp(2);
  VTh3(k)=Th3;  % Store solutions to vector
  VTh4(k)=Th4;
  % Use last solution for next guess
  Xinit=[Th3,Th4];
end
%%% End loops


%% Animation
% Determine the positions of the four 
% revolute joints at each iteration
% x Position of the origin

Ox = zeros(1,length(VTh3));
% y position of the origin
Oy = zeros(1,length(VTh3));
% x coordinate for point B
Bx = real(r1*exp(1i*Th1*ones(length(VTh3)) ));
% y coordinate for point B
By = imag(r1*exp(1i*Th1*ones(length(VTh3)) ));
% x coordinate for point C
Cx = real(r2*exp(1i*VTh3));
% y coordinate for point C
Cy = imag(r2*exp(1i*VTh3));
% x coordinate for for point D
Dx = real(r3*exp(1i*VTh3))+Cx;
% y coordinate for for point D
Dy = imag(r3*exp(1i*VTh3))+Cy;
% x coordinate for point D
Ex = real(r4*exp(1i*VTh4))+Dx;
% x coordinate for point D
Ey = imag(r4*exp(1i*VTh4))+Dy; 
% x coordinate for point E
Gx = real(lr6*exp(1i*(VTh4+Th_GED)))+Ex;
% y coordinate for point E
Gy = imag(lr6*exp(1i*(VTh4+Th_GED)))+Ey; 

% Find Plot Limits
xmin = min([min(Ox) min(Bx) min(Cx) min(Dx) min(Ex) min(Gx)]);
xmax = max([max(Ox) max(Bx) max(Cx) max(Dx) max(Ex) max(Gx)]);
ymin = min([min(Oy) min(By) min(Cy) min(Dy) min(Ey) min(Gy)]);
ymax = max([max(Oy) max(By) max(Cy) max(Dy) max(Ey) max(Gy)]);
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
    bar2x=[Ox(t) Cx(t)]; % Coordinates Link 2
    bar2y=[Oy(t) Cy(t)];
    bar3x=[Cx(t) Dx(t)]; % Coordinates Link 3
    bar3y=[Cy(t) Dy(t)];
    bar4x=[Dx(t) Ex(t)]; % Coordinates Link 4
    bar4y=[Dy(t) Ey(t)]; 
    plot(bar2x,bar2y,bar3x,bar3y,bar4x,bar4y,'LineWidth',2);   
    axis equal   % Equal scale of x and y-axis 
    axis([xmin xmax ymin ymax]);
    M(t)=getframe; % For assembling movie frames
end
%%% End Animation

