function [curAllX,clusterIndex]=CoherentFilter(trks,curTime,d,K,lamda)
%% Algorithm 1 of Coherent Filtering: Detecting coherent motion patterns at each frame
% Last revised by Bolei Zhou Jan.24,2013

%% step1: find K nearest neighbor set at each time

[allXset,allVset] = CF_trk2XV(trks,curTime,d);
curAllX=allXset{1,1};
nPoint=size(curAllX,2);
% display('find K nearest neighbor set at each time...');
[neighborSet,correlationSet] = CF_neighbor(allXset,allVset,K);

%% step2: find the invariant neighbor and pairwise connection set

% display('search invariant neighbor and construct the pairwise connections...');
[pairSet,correSet] = CF_neighbor2pair(neighborSet,correlationSet,d);

%% step3: threshold pairwise connection set by the averaged correlation values, then generate cluster components 

% display('threshold the average velocity correlations, and get the clustering...');
pairIndex=find(correSet>lamda); % included pairwise connection
clusterIndex=pair2cluster(pairSet(:,pairIndex),nPoint);


end

