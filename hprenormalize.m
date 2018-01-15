function X_inP = hprenormalize( X )
% x_inS = hprenormalize( x )
% Returns the normalized vector so that X(1,:)=1
% Works with column vectors

nh = size(X,1);
X_inP = X ./ repmat(X(1,:),nh,1);