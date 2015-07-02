function params = learn_word_lds(varargin)

addpath(genpath('.'));
addpath('../machineInfo/');

p = inputParser;
addOptional(p,'powerIterations',3,@isnumeric);
addOptional(p,'datasetName','news',@ischar);
addOptional(p,'iparams','',@ischar);
covStyles = {'full','diagonal','spherical','bohning'};
checkMember = @(x,y) any(validatestring(x,y));
checkCov = @(x) checkMember(x,covStyles);
addOptional(p,'covStyle','full',checkCov);
addOptional(p,'klimForASOS',10,@isnumeric);
addOptional(p,'klimForSSID',2,@isnumeric);
addOptional(p,'kappa',100,@isnumeric);
addOptional(p,'doSSID',true,@islogical);
addOptional(p,'hidsize',200,@isnumeric);
addOptional(p,'expt','lds',@ischar);
addOptional(p,'ssidTransformPower',0,@isnumeric);
addOptional(p,'randSeed',0,@isnumeric);
addOptional(p,'doPOS',true,@islogical);
addOptional(p,'randomInit',false,@islogical);
addOptional(p,'saveParams',true,@islogical);
addOptional(p,'writeRNNFiles',false,@islogical);
addOptional(p,'emLoops',20,@isnumeric);
addOptional(p,'rnnParamDir','',@ischar);
addOptional(p,'workingDir','',@ischar);



parse(p,varargin{:});
p.Results

numEMLoops = p.Results.emLoops;

rng(p.Results.randSeed);
saveParams = p.Results.saveParams;
doPOS = p.Results.doPOS;
randomInit = p.Results.randomInit;
updateInit = true;
numIters = 5; %num EM iters between when we evaluate
klim = p.Results.klimForASOS;
klimForSSID = p.Results.klimForSSID;
machineInfo = GetMachineInfo();

insize = 0;
kappa = p.Results.kappa;
hidsize = p.Results.hidsize;
simpleMethodForC = false;
ssidTransformPower = p.Results.ssidTransformPower;
useTropp = true;
powerIterations = p.Results.powerIterations;
datasetName = p.Results.datasetName;
iparamsFile = p.Results.iparams;
dataDir = ['./' datasetName];
covStyle = p.Results.covStyle;

updateTransMat = true;
updateHiddenCovMat = false;
updateObservedCovMat = true;
firstIterToUpdateEmbeddings = 2;
firstIterToUpdateTransitions = 2;
doSubspaceId = p.Results.doSSID;
loadFromFile = ~doSubspaceId && ~randomInit;
workingDir = p.Results.workingDir;
fprintf('working in %s\n',workingDir);
mkdir(['../' workingDir]);

CooccurrenceAtLagsFile = sprintf('%s/DataAtLags_Whitened_k%d.mat',dataDir,kappa);
fprintf('loading statistics from %s\n',CooccurrenceAtLagsFile);
vocabFile = [dataDir '/vocab.mat'];
outFile = ['../' workingDir '/params'];
load(vocabFile);
edgesize = 0;
ape = word_GetWordStatisticsFromFile(CooccurrenceAtLagsFile,klim, 2*klim+1, edgesize, insize);
ape.frequencies = ape.meanShift.^2;
T = ape.T;
disp('finished loading data');

%%%%%%%%%%%%%%%%%%%%%%%%
%%parameter initialization (either SSID or from a file of init values)
if(doSubspaceId)
    fsvdOverSampling = 10;
    n = length(ape.frequencies);
    transformForSVD = spdiags(ape.frequencies.^ssidTransformPower,0,n,n);
    iparams = WordSubspaceId(ape,hidsize,useTropp,fsvdOverSampling,powerIterations,ape.T,klimForSSID, simpleMethodForC,transformForSVD);
elseif(loadFromFile)
    fprintf('loading iparams from %s\n',iparamsFile);
    load(iparamsFile);
    iparams = params;
else
    assert(1==0,'you must either do SSID or initialize with parameters from a file');
end
%%%%%%%%%%%%%%%%%%%%%%%%


params = iparams;
disp('starting EM');
LLhist = [];
totalIters = 0;
doPOSOnFirst = ~loadFromFile;
for k = 1:numEMLoops
    if(doPOSOnFirst || k > 1)
        params = word_GetSteadyStateParams(params,ape);
        params = word_GetParamsForFiltering(params,ape);
        if(saveParams)
            save([outFile '-' num2str(k)],'params','LLhist','totalIters','-v7.3');       
        end
        
        modes = {'kalman'};
        if(doPOS)
            for mi = 1:length(modes);
                mode = modes{i};
                POS('dev',mode,num2str(k),dataDir,workingDir,params,vocab,[],useUniversalTags,posDebugMode,posConfigs,machineInfo,true,false,false,true,true,'',-1);
            end
        end
        
        if(k == 1)
            embeddings_model = normalizeColumns(orthogonalizeColumns(ape.whiten * params.C*inv(params.V_0_1))')';
        else
            embeddings_model = normalizeColumns(orthogonalizeColumns(ape.whiten * params.C*inv(params.V_0_1))')';
            % embeddings_model = normalizeColumns(orthogonalizeColumns(ape.whiten * params.Ey_x_0*inv(params.V_0_1))'   )';
        end
        embeddings_C = normalizeColumns(params.C')';
        %the indices in the next two lines are arbitrary. this is just so that you don't do Nearest Neighbors for the whole vocabulary
        PrintNearestNeighbors({embeddings_model,embeddings_C},{'model','C'},vocab,1000:1200);
        WriteNN(embeddings_model,k,workingDir,vocab,1000:2500); 
    end
    
    if(k >= 1)
        firstIterToUpdateEmbeddings = 0;
        firstIterToUpdateTransitions = 0;
    end
    fprintf('starting EM stage %d\n',k)
    [params, LLhist1] = word_ApproxEMforLDS( params,[], numIters,updateTransMat, updateHiddenCovMat, updateObservedCovMat,updateInit,firstIterToUpdateEmbeddings,firstIterToUpdateTransitions , ape , T,covStyle);
    fprintf('finished EM stage %d\n',k)
    totalIters = totalIters + numIters;
    
    LLhist = [LLhist LLhist1];
    oo = fopen(sprintf('../%s/nn-%d-ll.txt',workingDir,k),'w');
    fprintf(oo,'%f ',LLhist(:));
    fclose(oo);
end


if(p.Results.writeRNNFiles)
    params.T = T;
   if(numEMLoops == 0)
      params.Ey_x_0 = zeros(size(params.C));
      params.Ex_x_0 = eye(size(params.A,1));
      
   end
   params = word_GetSteadyStateParams(params,ape);
   params = word_GetParamsForFiltering(params,ape);
   fprintf('writing %s\n',[outFile '-' 'final']);
    save([outFile '-' 'final'],'params','LLhist','totalIters','-v7.3');       

   outBase =  [p.Results.rnnParamDir '/rnn-' p.Results.expt];
   WriteRNNFiles([dataDir '/class_vocabIndices'],outBase,params,ape); 
end


