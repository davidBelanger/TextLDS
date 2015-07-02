function WriteVocab(vocab,outFile)
out = fopen(outFile,'w');
for i = 1:length(vocab)
    fprintf(out,'%s\n',vocab{i});
end
