%function [ cl0,cl_a,clMin,clMax,cd0,cd2u,cd2l,clcd0 ] = polarAnalyzer( polarFile )
%Finds the QProp parameters from the exported XFLR5 polars
%   Detailed explanation goes here

polarDirName = 'Analysis_Polars';
files = dir(fullfile(polarDirName, '*.txt'));
polarData3D = [];
reData = [];
for polarFile = transpose(files)
	try %TODO: Loop through polarData2D and interpolate to fill in missing rows/alphas before concatenating with polarData3D
		%	   Use F = fillmissing(A,method)
		[ReStartInd, ReEndInd] = regexp(polarFile.name, '(?<=Re).*?(?=_)'); %Match substring by regexp with reynolds number and output indices
		re = str2double(string(polarFile.name(ReStartInd : ReEndInd))); %Get reynolds num associated with file, in millions
		re = re*1000000;

		fopen(fullfile(polarDirName, polarFile.name));

		polarFileText = fileread(polarFile.name); %Convert file to text
		%Get all the text after the last '-\n' i.e. all the data
		findLastHdrChar = strfind(polarFileText,'---') + 2;
		polarDataString = extractAfter(polarFileText, findLastHdrChar(end));
		%Convert datastring to matrix after stripping text of newlines
		numCols = 10; %The number of columns in the data
		%Replace newlines with spaces, convert str to cell array to matrix, reshape, transpose to match original data
		polarData2D = transpose(reshape(cell2mat(textscan(strrep(polarDataString,newline,''),'%f')),numCols,[]));

		polarData3D = cat(3, polarData3D, polarData2D); %Append data matrix for this Reynolds num to multidim data matrix
		reData = [reData, re]; %#ok<AGROW> %Append Reynolds num to Reynolds num vector
	catch ME
		if (strcmp(ME.identifier,'MATLAB:catenate:dimensionMismatch'))
			msg = ['Dimension mismatch occurred: main polars data matrix has '...
				num2str(size(polarData3D,1)),' rows while ', polarFile.name, ' has ', ...
				num2str(size(polarData2D,1)),' rows.'];
			causeException = MException('MATLAB:myCode:dimensions',msg);
			ME = addCause(ME,causeException);
		end
		rethrow(ME)
	end
end
	polarFileName = 'NACA 2412_T1_Re0.500_M0.00_N9.0.txt';
	
	
	cl0 = 0; cl_a = 0; clMin = 0; clMax = 0; cd0 = 0; cd2u = 0; cd2l = 0; 
	clcd0 = 0;
%end

