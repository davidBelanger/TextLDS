function getCoocurrenceAtLags(workingDir,id,V,filenameList)
%V = 150001;

%filenameList = 'matFileList';
fid = fopen(filenameList);
filenames = textscan(fid,'%s','delimiter','\n');
filenames = filenames{1};
len = length(filenames);
numLags = 10;

dataAtLags = cell(1,len);
 for i = 1:(numLags+1)
    dataAtLags{i} = sparse(V,V); 
 end
 
totalWords = 0;

for fi = 1:len
    fi
    s = load(filenames{fi},'ints');
    data = double(s.ints)+1;
    data(find(data == 150002)) = V;
    ld = length(data);
    totalWords = totalWords + ld;

    %%this is dumb, since adding to a sparse matrix is slow. should just
    %%get a massive index list once and then call the constructor once. 
    for lagP1 = 1:(numLags + 1)
        lag = lagP1 - 1
     dataAtLags{lagP1} = dataAtLags{lagP1} + sparse(data(1:end-lag),data(1+lag:end),ones(1,ld-lag),V,V);
     
        %  dataAtLags{lagP1} = dataAtLags{lagP1} + sparse(data(1+lag:end),data(1:end-lag),ones(1,ld-lag),V,V);
         end       
end

disp('saving');
save(sprintf('%s/DataAtLags_%s',workingDir,id),'totalWords','dataAtLags','-v7.3');

