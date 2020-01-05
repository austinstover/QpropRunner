function [ cl0,cl_a,clMin,clMax,cd0,cd2u,cd2l,clcd0] = polarAnalyzer3(polarFile, clBounds, alphas, REref, Mach)
%Finds the QProp parameters from the exported XFLR5 polars
%   Detailed explanation goes here

%Res = linspace(2e4, 1e5, 3);
%RE_SCALING_FACTOR = 1e6; %Divides Re by this for storage in arrays/matrices

%polar = []; %Initialize the list of points
%REref = Res(ceil(end/2)); %The reference reynolds value
% for Re = Res %Iter over Reynold's nums
% 	[pol,foil] = xfoil(airfoil,alphas,Re,Mach,'oper iter 50'); %Run xfoil for 1 Re
% 	%polar = [pol.alpha pol.CL pol.CD pol.CDp pol.Cm, pol.Top_xtr, pol.Bot_xtr];
% 	
% 	if(Re == REref) %Get the polar of a representative (the middle) Re val  TODO: Equate indices instead
% 		polarRef = [pol.alpha, pol.CD, pol.CL];
% 		REref = Re;
% 	end
% 	
% 	%CDList columns = [Re, alpha, CD, CL]
% 	%polar = vertcat(polar, [Re/RE_SCALING_FACTOR*ones(size(pol.alpha)), pol.alpha, pol.CD, pol.CL]);
% end

% I moved alpha definition to outside of polarAnalyzer3 function

% % FOR TESTING
% polarFile = 'NACA4412';
% clMinBound = -0.3;
% clMaxBound = 0.9;
% alphaMin = -10;
% alphaMax = 50;
% alphaStep = 0.4;
% REref = 9.8011e+04;
% Mach = 0.2705;
% clBounds = [clMinBound, clMaxBound];
% alphas = alphaMin:alphaStep:alphaMax;


airfoil = polarFile;

Re = REref; %2e4;
[pol,foil] = xfoil(airfoil,alphas,Re,Mach,'oper iter 50'); %Run xfoil for 1 Re
polarRef = [pol.alpha, pol.CD, pol.CL];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot and Analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%finding beta
% p = pitch
% D = diameter
% x = fraction along propeller
%beta = atan(p/(D*pi*x));

%finding chord as function of radius
%???

% CL(alpha)  =  ( CL0  +  CL_a*alpha ) / beta    ,  clipped to  CLmin..CLmax  range
% beta = sqrt(1 - M^2)  is the local Prantdl-Meyer compressibility factor
% (1 for ideal gas)
% x = alpha, y = CL

%BUG: CL can hit max after it starts decreasing, then increasing again,
%fooling the algorithm into thinking the max/min happen at the wrong pts
[~,CLMinInd] = min(polarRef(:,3)); %TODO: Change from min to where CL starts increasing
[~,CLMaxInd] = max(polarRef(:,3)); %TODO: Change from max to where CL stops increasing

[startInd, stopInd] = longestTrue(diff(polarRef(CLMinInd:CLMaxInd,3)) > 0); %Find longest consecutive series of increasing vals btwn the min and max inds

CLStartInd = startInd + CLMinInd - 1; %CL starts its consecutive increase
CLStopInd = stopInd  + CLMinInd; %CL stops its consecutive increase

clMin = polarRef(CLStartInd,3);
clMax = polarRef(CLStopInd,3);
noSepInd = 1:length(polarRef(:,3))>=CLStartInd & 1:length(polarRef(:,3))<=CLStopInd; %Range clip; deletes values post stall
[polarRefNoSep1,polarRefNoSep2] = prepareCurveData(polarRef(noSepInd,1),polarRef(noSepInd,3));
CLAlphaFit = fit(polarRefNoSep1,polarRefNoSep2,'poly1');

figure;
hold on
plot(CLAlphaFit, polarRef(:,1), polarRef(:,3));
xlabel('\alpha');
ylabel('CL');
title(sprintf('%s at Re=%g',airfoil, Re));
CLAlphaLim = ylim;
plot(polarRef(find(noSepInd,1,'first'),1)*[1,1],CLAlphaLim,'-');
plot(polarRef(find(noSepInd,1,'last'),1)*[1,1],CLAlphaLim,'-');
hold off
legend('plot','lin fit','left stall','right stall','location','best');

beta = sqrt(1 - Mach*Mach); % The local Prantdl-Meyer compressibility factor
CLAlphaFitCoeffs = coeffvalues(CLAlphaFit);
clMinBound = min(clBounds); %I added a function input to specify max and min cl to fit to
clMaxBound = max(clBounds);

cl_a = CLAlphaFitCoeffs(1)*beta;
cl0 = CLAlphaFitCoeffs(2)*beta;

%CLCD0 is the CL at the lowest CD
[~,CDMinInd] = min(polarRef(noSepInd,2));
CLNoSep = polarRef(noSepInd,3);
clcd0 = CLNoSep(CDMinInd); %The CLCD0 derived from the CD plot
CLInBounds = polarRef(:,3) > clMinBound & polarRef(:,3) < clMaxBound;
polarRefULog = (polarRef(:,3) > clcd0); %The boolean logical vector for which indices have CL > CLCD0
polarRefU = polarRef(polarRefULog & noSepInd' & CLInBounds,:);
polarRefL = polarRef(~polarRefULog & noSepInd' & CLInBounds,:);

%Fit CD(CL) w/ middle Re polar
% CD(CL)  =  CD0 + CD2*(CL-CLCD0)^2 = (CD0 + CD2*CLCD0^2) + (2CLCD0*CD2)CL + (CD2)CL^2
CDCLFitU = fit(polarRefU(:,3), polarRefU(:,2),'poly2');
CDCLFitL = fit(polarRefL(:,3), polarRefL(:,2),'poly2');
figure;
hold on;
plot(CDCLFitU, polarRefU(:,3), polarRefU(:,2),'.r');
plot(CDCLFitL, polarRefL(:,3), polarRefL(:,2),'.b');
ylabel('CD'); xlabel('CL');
title(sprintf('%s at Re=%g',airfoil, Re));

yLimMax = max(max(polarRefU(:,2)),max(polarRefL(:,2)));
yLimMin = min(min(polarRefU(:,2)),min(polarRefL(:,2)));
plot([clcd0,clcd0],[yLimMin,yLimMax]);
ylim([yLimMin, yLimMax]);
legend('CDU','CDCLFitU','CDL','CDCLFitL','CLCD0');
hold off;
CDCLFitCoeffsU = coeffvalues(CDCLFitU);
CDCLFitCoeffsL = coeffvalues(CDCLFitL);

% CD(CL,Re)  =  [ CD0 + CD2*(CL-CLCD0)^2 ] * [Re/REref]^REexp
if(clcd0 <= 0) %y-intercept occurs on upper parabola
    CD2 = CDCLFitCoeffsU(1); %The x^2 term
    cd0 = CDCLFitCoeffsU(3) - CD2*clcd0*clcd0; %cd0 is y-intercept on CL-CD plot (drag at 0 lift)
else %y-intercept occurs on lower parabola
    CD2 = CDCLFitCoeffsL(1);
    cd0 = CDCLFitCoeffsL(3) - CD2*clcd0*clcd0;
end

cd2u = CDCLFitCoeffsU(1); %cd2 = cd2u where cl>clcd0
cd2l = CDCLFitCoeffsL(1); %cd2 = cd2l where cl<clcd0


% c = 0;
% g1 = fittype( @(a,b,x) a*x.^2+b*x+c )

% Params to find: CD0    CD2u   CD2l   CLCD0 REref  REexp
% ft = fittype( @(CD2,REexp) )

% CD(CL,Re)  =  [ CD0 + CD2*(CL-CLCD0)^2 ] * [Re/REref]^REexp
% ft = fittype('(CD0 + CD2*(CL-CLCD0)^2) * (Re/REref)^REexp',...
%             'problem','CD0','CD2','CLCD0',...
% 			'dependent','CL','Re','REref','REexp',...
%             'independent','CD');
% f2 = fit( [CDList(:,1), CDList(:,4)], CDList(:,3), 'poly23' )
% figure;
% plot(f2, [CDList(:,1), CDList(:,4)], CDList(:,3));