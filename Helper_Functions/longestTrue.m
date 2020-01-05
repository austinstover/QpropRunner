% Author: Austin Stover
% Date: January 2020
% Finds longest consecutive string of 1s in logical array (or nonzeros in
% numerical array). If no 1s found, returns 0s for start and stop index.
function [startInd, stopInd] = longestTrue(lArr)
    longest = 0;
    numCons = 0;
    potentialStartInd = 1;
    startInd = 0;
    stopInd = 0;
    for i = 1:length(lArr)
        if(lArr(i))
            if(numCons == 0)
                potentialStartInd = i;
            end
            numCons = numCons + 1;
            if(numCons > longest)
                longest = numCons;
                startInd = potentialStartInd;
                stopInd = i;
            end
        else
            numCons = 0;
        end
    end
end