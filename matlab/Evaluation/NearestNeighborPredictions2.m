function idxs = NearestNeighborPredictions2(queryIdx,C,num,indsToIgnore)

pred = C*C(queryIdx,:)';
%    pred = params.Cscaled*params.A*params.steadyStateParams.Kscaled(:,queryIdx);


if(nargin == 4)
    pred(indsToIgnore) = min(pred)-100000000;
end

[ps pi] = sort(pred,'descend');

idxs = pi(2:(num+1));
