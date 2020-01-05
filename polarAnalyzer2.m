%function [ cl0,cl_a,clMin,clMax,cd0,cd2u,cd2l,clcd0 ] = polarAnalyzer( polarFile )
%Finds the QProp parameters from the exported XFLR5 polars
%   Detailed explanation goes here

airfoil = 'NACA2412';
Mach = 0;
alphas = -12:12;
Res = linspace(1e5, 1e6, 5);

LEN_ROW_POLAR = 7;
polarList = zeros(length(alphas),LEN_ROW_POLAR,length(Res)); %Initialize the 3D matrix of polar data
nanRow = NaN(1,LEN_ROW_POLAR);
for i = 1:length(Res) %Iter over Reynold's nums
	[pol,foil] = xfoil(airfoil,alphas,Res(i),Mach,'oper iter 50'); %Run xfoil for 1 Re
	polar = [pol.alpha pol.CL pol.CD pol.CDp pol.Cm, pol.Top_xtr, pol.Bot_xtr];
		
	if(length(alphas') > length(pol.alpha)) %If one analysis didn't converge
		for j = 1:length(alphas)
			if(j > size(polar, 1) || polar(j,1) ~= alphas(j))
				polar = insertrows(polar, nanRow, j-1); %Fill in gaps with NaNs
			end
		end
	end
	polarList(:,:,i) = polar;
end

