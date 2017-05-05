function idxs = NearestNeighborPredictions(queryIdx,params,num,bigram)

%%todo: this is completely ignoring the effect of x_\infty

if(~bigram)
    pred = params.C*params.K(:,queryIdx);
else
    pred = params.C*params.A*params.K(:,queryIdx);
end

[ps pi] = sort(pred,'descend');
idxs = pi(2:(num+1));
