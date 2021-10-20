function [err,nI, resid]=findSimilarity(s,y,reduceOrder, ni, na)
%
% Usage : [err,nI, resid]=findSimilarity(s,y)
%
% Find the most similar vector to s in y without scaling.
% DFT approaches exist which are based on the expansion of (s-y)^2 = 2^2 - 2sy + y^2. These approaches
% are most effective when the vectors s and y are large.
% This implementation does a stepped search for a global minimum of s in a reduced order y, then does 
% a refined search to find the global minimum.
% Author : Matt Flax <flatmax@>
% Copyright October 2021
% BSD License.
% inputs :
%          s : The search vector
%          y : The longer data vector
%reduceOrder : Whether to do the first stage approximation or the refinement
%         ni : The mInimum n variable to start a refined search from
%         na : The mAximum n variable to stop searching at
% ouptuts :
%        err : The RMS error between s and y at nI
%         nI : The first index of the matched vector
%      resid : The residual vector from the match

if nargin<3 % run in reduced order mode to find the likely global min, then recurse one more time to do a local search
	reduceOrder=1;
end
Ns=size(s,1);
Ny=size(y,1);

Noverlap=round(5*Ns/6);
if ~reduceOrder
	%     ni=1;
	if nargin>=4 & ~isempty(ni)
		Noverlap=Ns-1;
		%          yOrig=y;
		if na>length(y)
			na=min(na,length(y));
			ni=round(length(y)-1.5*Ns);
		end
	else
		ni=1; na=length(y);
	end
	res=buffer(y(ni:na),Ns,Ns-1,'nodelay');
	res=s-res;
	%     errM=rms(res); % this is more expensive to compute then mean square
	%     [err,nI]=min(errM);
	errM3=mean(res.^2);
	[err,nI]=min(errM3);
	resid=res(:,nI);
	err=sqrt(err);
	nI=nI+ni-1;
	%     if nargin>=4 % use this for testing
	%         [err2,nI2, resid2]=findSimilarity(s,yOrig,0);
	%     end
	%
else % reduced order estimate
	Nskip=Ns-Noverlap;
	
	res=buffer(y,Ns,Noverlap,'nodelay');
	res=s-res;
	%     errM=rms(res); % too expensive, use mean square
	% 	[~,nIEst]=min(errM);
	errM3=mean(res.^2);
	[~,nIEst]=min(errM3);
	reduceOrder=0;
	[err,nI, resid]=findSimilarity(s,y,reduceOrder, (nIEst-1)*Nskip+1, nIEst*Nskip+1+Ns);
	
end
end
