function o = word_GetWordStatisticsFromFile(fname,klim, klag, edgesize, insize)

%%we're assuming the input data being loaded has been whitened
load(fname,'dataAtLags','unWhiten','whiten','meanShift','T');
o = struct;

 total = klim+1+1;
 outsize = size(dataAtLags{1},1);
 fu_u_ =zeros( insize, insize, total);
 fu_y_ = zeros( insize, outsize, total );
 fy_u_ = zeros( outsize, insize, total );
 
 n = size(dataAtLags{1},1);

o.T = T;
 
o.y_y_0_with_counts_scaled =  dataAtLags{1}/T;

fy_y_ = dataAtLags(1:(klim+1));

o.y_y_0 = fy_y_{1};
o.klag = klag;
o.edgesize = edgesize;

%this is just a dummy to get the indexing to line up
 T = 1000;

o.y_lead_ = sparse(outsize,edgesize);
o.u_lead_ = sparse(insize,edgesize);
len = length(T-edgesize+1:T);
o.y_trail_ = sparse(outsize,len);
o.u_trail_ = sparse(insize,len);

len2      = length(edgesize+1:edgesize+klag);
len3      = length(T-klag+1-edgesize:T-edgesize);
o.y_pre_  = sparse(outsize,len2);
o.u_pre_  = sparse(insize,len2);
o.y_post_ = sparse(outsize,len3);
o.u_post_ = sparse(insize,len3);

o.klim = klim;
o.y_y_ = fy_y_;
o.u_u_ = fu_u_;
o.u_y_ = fu_y_;
o.y_u_ = fy_u_;

 
o.y_u_0 = sparse(outsize,insize);
o.u_u_0 = sparse(insize,insize);
o.u_u_0_2 = sparse(insize,insize);
o.inv_u_u_0 = inv(o.u_u_0);
o.inv_u_u_0_2 = inv(o.u_u_0_2);

o.unWhiten = unWhiten;
o.whiten = whiten;
o.meanShift = meanShift;


