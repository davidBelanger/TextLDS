function WriteBinaryVocab(vocab,file)
o = fopen(file,'w');
for i = 1:length(vocab)
   fprintf(o,'%s ',vocab{i});
end

fclose(o);