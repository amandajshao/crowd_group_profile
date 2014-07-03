function [center, label] = litekmeans(X, k)
% Perform k-means clustering.
%   X: d x n data matrix
%   k: number of seeds
% Written by Michael Chen (sth4nth@gmail.com).
n = size(X,2);
last = 0;
label = ceil(k*rand(1,n));  % random initialization
kk = 0;
while any(label ~= last)
    [u,~,label] = unique(label);   % remove empty clusters
    k = length(u);
    E = sparse(1:n,label,1,n,k,n);  % transform label into indicator matrix
    m = X*(E*spdiags(1./sum(E,1)',0,k,k));    % compute m of each cluster
    last = label'; % modification
    [~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1); % assign samples to the nearest centers
    fprintf('.');
    kk = kk+1;
    if mod(kk, 50) == 0, fprintf('\n'); end
end
[~,~,label] = unique(label);
center = m;
disp('k means clustering has been done!');