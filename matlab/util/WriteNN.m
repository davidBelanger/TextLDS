function WriteNN(embeddings_norm,k,workingDir,vocab,range,info)
if(nargin == 5)
    info = '';
end
outname = sprintf('../%s/nn-%d-%s.txt',workingDir,k,info);
out = fopen(outname,'w');
for wordIdx = range
    nn = NearestNeighborPredictions2(wordIdx,embeddings_norm,5);
    fprintf(out,'%s --> %s %s %s %s %s\n',vocab{wordIdx}, vocab{nn(1)}, vocab{nn(2)},vocab{nn(3)}, vocab{nn(4)}, vocab{nn(5)});
end
fprintf('wrote nearest neighbors to %s\n',outname);

