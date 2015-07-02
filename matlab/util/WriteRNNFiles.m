function WriteRNNFiles(vocabIdxFile,base,params,ape)
Ey_x_0 = ape.unWhiten*params.Ey_x_0; %%bc the LDS (and thus the posterior, might be w.r.t whitened coordinates, but we want one-hots). 
Ex_x_0 = params.Ex_x_0;

%construct the class conditional matrix
%word_idx count class_idx
a = dlmread(vocabIdxFile,'\t');
%%account for 0 indexing
a(:,1) = a(:,1) + 1;
a(:,3) = a(:,3) + 1;
classes = a(:,3);
freqs = a(:,2)/sum(a(:,2));
numWords = size(a,1);
numClasses = max(classes);

classConditionalProbs = zeros(numClasses,numWords);
for i = 1:numClasses
   inds_in_class = find(classes == i);
   classConditionalProbs(i,inds_in_class) = freqs(inds_in_class); 
   classConditionalProbs(i,:) = classConditionalProbs(i,:)/sum(classConditionalProbs(i,:));
end

%classMatrix = zeros(200,200); 

%wordIdx classId cnt

Eclass_x_0 = classConditionalProbs*Ey_x_0;
classMatrix = (Eclass_x_0)*inv((Ex_x_0));

size(classMatrix)
   classMatrix(1:10)
fprintf('max = %f, sum = %f, min = %f\n',max(classMatrix(:)),sum(classMatrix(:)),min(classMatrix(:)));

for whiten = [true false]
    for hidScale = [1 2]
        for obsScaleExp = [0 1 2 3 4]
            obsScale = 10^(obsScaleExp);
            outfile = sprintf('%s-whiten-%d-hid-%f-obs-%f.lds',base,whiten,hidScale,obsScaleExp);
            fprintf('writing %s\n',outfile);
            PrintRNNFormat(outfile,params,whiten,hidScale,obsScale,classMatrix)
        end
    end
end
