function X_inP = hnormalize( X )
% x_inS = hnormalize( x )
% Returns the normalized vector so that X(end,:)=1
% Works with column vectors

nh = size(X,1);
X_inP = X ./ repmat(X(end,:),nh,1);