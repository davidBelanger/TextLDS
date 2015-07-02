function getWhitenedCoocurrenceAtLags(dataDir,kappa)
doSqrt = 0;
numLagsToTake = 10;
doSin = 0;

if(doSin)
    secondTransform = sinTransform;
else
    secondTransform = @(x,T) x;
end

CooccurrenceAtLagsFile = [dataDir '/DataAtLags.mat'];
fprintf('loading from %s', CooccurrenceAtLagsFile);
load(CooccurrenceAtLagsFile);
disp('loaded');

n = size(dataAtLags{1},1);
dataAtLags = dataAtLags(1:numLagsToTake);
sqrtTransform = @(x) sqrt(x);
%sqrtTransform = @(x) x.^(1/3);
if(doSqrt)
    for i = 1:numLagsToTake
       dataAtLags{i} =  sqrtTransform(dataAtLags{i});
    end
end

dataAtLags{1} = dataAtLags{1} + kappa*speye(n);
rawDataAtLags = dataAtLags;
counts = spdiags(dataAtLags{1});

T = sum(counts);
frequencies = counts/T;
whiteningTransform = spdiags(1./sqrt(frequencies),0,n,n);
whitenedMean = whiteningTransform*frequencies;
disp('whitening');
%todo: check that everything lines up

dataAtLags = dataAtLags(1:numLagsToTake);
for i = 1:numLagsToTake
    if(sum(dataAtLags{i}(:)) == 0)
        break;
    end
   dataAtLags{i} =  secondTransform(whiteningTransform*rawDataAtLags{i}*whiteningTransform,T);
end

load('news/vocab');
disp('perfectly correlated words');
z = dataAtLags{2}/T;
[is js] = find(z == 1.0);
t = length(is);

for tt = 1:t
    tt
   fprintf('%s %s %d %d\n',vocab{is(tt)},vocab{js(tt)},full(dataAtLags{1}(is(tt),is(tt))),full(dataAtLags{1}(js(tt),js(tt)))); 
end

disp('checking');
for i = 1:numLagsToTake
    if(sum(dataAtLags{i}(:)) == 0)
        break;
    end
    whiteCov = (dataAtLags{i}/T);
%    vec = whitenedMean'*(whiteCov - whitenedMean*whitenedMean');
    vec = whitenedMean'*whiteCov  - (whitenedMean'*whitenedMean)*whitenedMean';
    norm(vec)
 %   assert(norm(vec) < 0.0001);
end

unWhiten = inv(whiteningTransform);
whiten = whiteningTransform;
meanShift = whitenedMean;

save(sprintf('%s/DataAtLags_Whitened_k%d.mat',dataDir,kappa,doSqrt,doSin),'dataAtLags','unWhiten','whiten','meanShift','T','-v7.3');

