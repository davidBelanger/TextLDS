function params = word_GetParamsForFiltering(params,ape)
n = size(params.C,1);

filteringParams.pi_1 = params.pi_1;
filteringParams.A = params.A;

wordSpace2NormalizedSpace = ape.whiten; 
normalizedSpace2WordSpace = ape.unWhiten;

filteringParams.K_noshift = params.K*wordSpace2NormalizedSpace;
Km = params.K*ape.meanShift(:);
filteringParams.K = filteringParams.K_noshift - repmat(Km,[1 size(params.K,2)]);

filteringParams.C = normalizedSpace2WordSpace*params.C;
filteringParams.J = params.J;
filteringParams.KC = filteringParams.K_noshift*filteringParams.C; %%todo: is it right to use K_noshift here?

filteringParams.SinvC  =  wordSpace2NormalizedSpace*MIL_rightMultiply(params.invS,filteringParams.C);
filteringParams.Ctranspose_Sinv_C = MIL_innerMultiply(params.invS,params.C',params.C);
filteringParams.likTerm = MIL_diag(params.invS); %%todo: account for meanshift in this (it's just a constant factor, though)

filteringParams.A_KCA = params.A - filteringParams.KC*params.A;
filteringParams.invS_logDet = params.invS.logDet + 2*sparseDiagLogdet(normalizedSpace2WordSpace);

params.filteringParams = filteringParams;

fprintf('max ev of A - KCA: %f\n',max(abs(eig(filteringParams.A_KCA))));
fprintf('max ev of J: %f\n',max(abs(eig(params.filteringParams.J))));