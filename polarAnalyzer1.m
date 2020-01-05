polarFile = 'Analysis_Polars/NACA 2412_T1_Re0.015_M0.00_N9.0.txt';

fopen(fullfile(polarFile, polarFile));

polarFileText = fileread(polarFile); %Convert file to text
%Get all the text after the last '-\n' i.e. all the data
findLastHdrChar = strfind(polarFileText,'---') + 2;
polarDataString = extractAfter(polarFileText, findLastHdrChar(end));
%Convert datastring to matrix after stripping text of newlines
numCols = 10; %The number of columns in the data
%Replace newlines with spaces, convert str to cell array to matrix, reshape, transpose to match original data
polarData2D = transpose(reshape(cell2mat(textscan(strrep(polarDataString,newline,''),'%f')),numCols,[]));

alpha = polarData2D(:,1);
cl = polarData2D(:,2);

cl_alpha_fit = fit(alpha, cl,'poly2');
plot(cl_alpha_fit, alpha, cl)