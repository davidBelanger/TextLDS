function AA = ConstructHankel(ape,order,whiten,numObs)
n = size(ape.y_y_{1},1);

klim = 2*order;
is = cell(1,klim);
js = cell(1,klim);
vs = cell(1,klim);

for i = 1:klim %%todo: you might be able to make this slightly smaller (i.e. you don't use all of them later)
    [ii jj vv] = find(whiten*ape.y_y_{i}*whiten');
    
    is{i} = ii';
    js{i} = jj';
    vs{i} = vv';
end

disp('got elements');
len = 0;
for row = 1:order
    for col = 1:order
        covIdx = order -(col-1) + (row-1)+1;
        len = len + length(is{covIdx});
    end
end
disp(['got length of nz vector: ' num2str(len)]);

rowInds = zeros(1,len,'int32');
colInds = zeros(1,len,'int32');
values = zeros(1,len);
offset = 0;

start = tic;
tic;
for row = 1:order
    for col = 1:order
        rowIdxShift = (row-1)*n;
        colIdxShift = (col-1)*n;
        covIdx = order -(col-1) + (row-1)+1;
        len = length(is{covIdx});
        rowInds(offset+1:(offset + len)) = is{covIdx} + rowIdxShift;
        colInds(offset+1:(offset + len)) = js{covIdx} + colIdxShift;
        values(offset+1:(offset + len)) =  vs{covIdx};
        offset = offset + len;
        fprintf('%02d ',covIdx-1);
    end
    fprintf('\n');
end

disp('got indices');
disp('constructing matrix');
AA = sparse(double(rowInds),double(colInds),values/numObs,n*order,n*order);
fprintf('time to construct matrix: %f\n',toc);
