function [clusteredLabels, Q, W, affin, rk] = graphAgglomerativeClustering (distance_matrix, groupNumber, strDescr, K, p, a, usingKcCluster, X, normalize)
%% Graph Agglomerative Clustering toolbox
% Input: 
%   - distance_matrix: pairwise distances, d_{i -> j}
%   - groupNumber: the final number of clusters
%   - strDescr: structural descriptor. The choice can be
%                 - 'gdl': graph degree linkage algorithm
%                 - others to be added
%   - K: the number of nearest neighbors for KNN graph, default: 20
%   - p: merging (p+1)-links in l-links algorithm, default: 1
%   - a: for covariance estimation, default: 1 
%       sigma^2 = (\sum_{i=1}^n \sum_{j \in N_i^3} d_{ij}^2) * a
% Output:
%   - clusteredLabels: clustering results
% by Wei Zhang (wzhang009 at gmail.com), June, 8, 2011
%
% Please cite the following paper, if you find the code is helpful
%
% W. Zhang, X. Wang, D. Zhao, and X. Tang. 
% Graph Degree Linkage: Agglomerative Clustering on a Directed Graph.
% in Proceedings of European Conference on Computer Vision (ECCV), 2012.

%% parse inputs
% disp('--------------- Graph Agglomerative Clustering ---------------------');

if nargin < 2, error('GAC: input arguments are not enough!'); end
if nargin < 3, strDescr = 'gdl'; end
if nargin < 4, K = 20;   end
if nargin < 5, p = 1;    end
if nargin < 6, a = 1; end
if nargin < 7, usingKcCluster = true; end
if nargin < 9, normalize = true; end

%% initialization
%     disp('---------- Building graph and forming initial clusters with l-links ---------');
    [graphW, NNIndex] = gacBuildDigraph_c(distance_matrix, K, a);
    initialClusters = gacBuildLlinks_cwarpper(distance_matrix, p, NNIndex);
    clear distance_matrix NNIndex

%     disp('-------------------------- Zeta merging --------------------------');
        if usingKcCluster
            clusteredLabels = gdlMergingKNN_c(graphW, initialClusters, groupNumber);
        else
            [clusteredLabels, Q, W, affin, rk] = gdlMerging_c(graphW, initialClusters, groupNumber, X, normalize);
        end

end