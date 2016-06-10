function PrintNearestNeighbors(embeddingsToTry,embeddingNames,vocab,inds,indsToIgnore,matrices)
if(nargin ~= 5)
    indsToIgnore = [];
end
useMatrices = nargin == 6;
for wordIdx = inds
    for i = 1:length(embeddingsToTry)
        if(useMatrices)
                             nn = NearestNeighborPredictionsFromQuadForm(wordIdx,embeddingsToTry{i},5,matrices{i});    
        else
                 nn = NearestNeighborPredictions2(wordIdx,embeddingsToTry{i},5,indsToIgnore);    
        end
        fprintf('%s --> %s %s %s %s %s   (%s)\n',vocab{wordIdx}, vocab{nn(1)}, vocab{nn(2)},vocab{nn(3)}, vocab{nn(4)}, vocab{nn(5)},embeddingNames{i});
    end
end

