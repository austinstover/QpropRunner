%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Austin Stover, Caitlind Walker
% Date:   October 2018 - January 2020
% Script for running QProp; Creates files and runs QProp from command line
%
% Notes
% QPROP usage:
% qprop propfile motorfile Vel Rpm [ Volt dBeta Thrust Torque Amps Pele ]   (single-point)
% qprop propfile motorfile Vel1,Vel2,dVel Rpm ["]              (multi-point 1-parameter sweep over Vel, Rpm set)
% qprop propfile motorfile Vel1,Vel2,dVel 0 Volt ["]           (multi-point 1-parameter sweep over Vel, Volt set)
% qprop propfile motorfile Vel1,Vel2,dVel Rpm1,Rpm2,dRpm ["]   (multi-point 2-parameter sweep over Vel and Rpm)
% qprop propfile motorfile runfile                             (multi-point, via file specification)
%
% Running Command Line:
% https://www.mathworks.com/help/matlab/ref/system.html
% status = system(command)
% [status,cmdout] = system(command)
%
% Creating/opening file:
% mkdir tests
% file = "tests/new_script2.m";
% edit file
%
% Writing data to file:
% https://www.mathworks.com/help/matlab/ref/fprintf.html
% Write a short table of the exponential function to a text file called exp.txt.
% x = 0:.1:1;
% A = [x; exp(x)];
% fileID = fopen('exp.txt','w');
% fprintf(fileID,formatSpec,A1,...,An)
% fprintf(fileID,'%6.2f %12.8f\n',A);
% fclose(fileID);
%
% View file contents
% type exp.txt
%
% TODO: Create fluid constants file
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Propeller Properties
airfoil = 'NACA4412';
nBlades = 2;
D_in = 9;       %in, Diameter
p_in = 7;       %in, Pitch
velocity = 0;%20;  %m/s, Approximate forward airspeed
rpm = 10000;    %Approximate RPM        For v=20, rpm=10000, D_in=9: J=0.524

%Propeller Airfoil Analysis Parameters
clMinBound = -0.3;
clMaxBound = 1.0; %0.9;
alphaMin = -10;
alphaMax = 20;
alphaStep = 0.4;


% Motor Properties
m_name = "Sunnysky_380Kv";
kv = 380; %RPM/V
R = 0.1; %Ohms, Motor internal resistance
I0 = 0.3; %A, No-lod current


%Subdirectories in which to put Qprop init files
propSubDirName = "Propellers";
motorSubDirName = "Motors";

%Fluid constants TODO: Implement fluid constants Qprop file
a = 340.0; %m/s, Speed of sound
nu = 1.78e-5; %Kinematic viscosity, m^2/s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

radius = D_in*0.0254 / 2; %m
radius_in = D_in / 2; %in

%Implement estimate for prop chord over radial fraction
x_root = 0.25;
x = linspace(x_root,1,10); %Radial fraction from root to tip
r = radius_in*x;

c_r = -0.3463*x.^2 + 0.235*x + 0.1845; %Chord/rad (unitless) for Graupner CAM 6x3 Qprop example
c_in = c_r*radius_in;
cRoot = (-0.3463*(x_root)^2 + 0.235*(x_root) + 0.1845)*radius; %And for the root (m)

beta_deg = rad2deg(atan(p_in./(D_in*pi*x)));

%The Re num with which to run XFOIL. Pick reasonably given prop speed, at
%about 0.75*radius. RE = V/v, V = airfoil freestream speed with unit chord,
%v = freestream kinematic viscosity (see xfoil documentation for more info)

%TODO: Calculate correct airspeed
uRoot = sqrt((rpm/60 * 2*pi*radius*x_root)^2 + velocity^2); %Airspeed into prop blade element at 3/4 rad
REref = cRoot*uRoot/nu; %2e4;
REexp = -0.5; %"Picking REexp = -0.5 is reasonable for most low Re airfoils."

Mach = uRoot/a;

%Conversion factors
Rfac = 0.0254; %m/in  r_SI (m) = r*Rfac + Radd
Cfac = 0.0254; %m/in  c_SI (m) = chord*Cfac + Cadd
Bfac = 1.0; %b_SI (rad) = (beta*Bfac + Badd) * pi/180
Radd = 0;
Cadd = 0;
Badd = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN XFOIL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate blade airfoil polars
clBounds = [clMinBound, clMaxBound];
alphas = alphaMin:alphaStep:alphaMax;
[cl0,cl_a,clMin,clMax,cd0,cd2u,cd2l,clcd0] = polarAnalyzer3(airfoil, clBounds, alphas, REref, Mach);%fopen(polarFileName));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WRITE FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: Write to fluid constants file for non-sea-level flight and a different
% kinematic viscosity, etc.

%Propeller file
propTitle = sprintf('%s_%gx%g_%gm_s_%grpm',airfoil,D_in,p_in,velocity,rpm);
propFileName = fullfile(propSubDirName,strcat(propTitle,'.txt'));
propFile = fopen(propFileName,'w');

fprintf(propFile,'%s\n',propTitle);

fprintf(propFile,'\n%g\t%g\t! Nblades\t[ R ]\n',nBlades,radius_in);

fprintf(propFile,'\n%g\t%g\t! CL0\t CL_a\n',cl0,cl_a);
fprintf(propFile,'%g\t%g\t! CLmin\t CLmax\n',clMin,clMax);

fprintf(propFile,'\n%g\t%g\t%g\t%g\t! CD0\t CD2u\t CD2l\t CLCD0\n',cd0,cd2u,cd2l,clcd0);
fprintf(propFile,'%g\t%g\t! REref\t REexp\n', REref, REexp);

fprintf(propFile,'\n%g\t%g\t%g\t! Rfac\t Cfac\t Bfac\n', Rfac, Cfac, Bfac);
fprintf(propFile,'%g\t%g\t%g\t! Radd\t Cadd\t Badd\n', Radd, Cadd, Badd);

fprintf(propFile,'\n# r\t chord\t beta');
for i=1:length(r)
    fprintf(propFile,'\n%g\t%g\t%g', r(i), c_in(i), beta_deg(i));
    if(i == 1)
        fprintf(propFile,'\t! root station');
    end
    if(i == length(r))
        fprintf(propFile,'\t! tip station, also gives R = r if R is omitted from Line 2\n');
    end
end

%Motor file
motorFileName = fullfile(motorSubDirName,strcat(m_name,'.txt'));
motorFile = fopen(motorFileName,'w');

motor_type = 1; %BLDC or brushed DC
fprintf(motorFile,'%s\t! name\n',m_name);

fprintf(motorFile,'\n%d\t! motor type  (1 = permanent-magnet brushed or brushless DC motor)\n',motor_type);

fprintf(motorFile,'\n%g\t! motor parameter 1, R (Ohms) for motor type 1\n',R);
fprintf(motorFile,'%g\t! motor parameter 2, Io (Amps) for motor type 1\n',I0);
fprintf(motorFile,'%g\t! motor parameter 3, Kv (rpm/Volt) for motor type 1\n',kv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN QPROP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TODO: Loop over different advance ratios for dynamic thrust, rpm for static

% qprop propfile motorfile Vel Rpm [ Volt dBeta Thrust Torque Amps Pele ]   (single-point)
commandStr = sprintf('qprop.exe %s %s %g %g',propFileName,motorFileName,velocity,rpm);
status = system(commandStr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTERPRET RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TODO: Use 10,11,12,13 I think

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CREATE PROPELLER COEFFICIENTS FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc