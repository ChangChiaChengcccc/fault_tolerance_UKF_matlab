function [X0,X_plus_vec,X_minus_vec] = sigma_fun(L,lambda,P,X0)
    % convert dimension
    X0 = X0.';

    % bias
    bias = sqrt(L+lambda);
    bias = bias*chol(P);
    % X_plus_vec, X_minus_vec(12*12), X0(1*12)
    X_plus_vec = X0(ones(L,1),:)+bias;
    X_minus_vec = X0(ones(L,1),:)-bias;
    
    % transpose X0(12*1) X_plus_vec(12*12) X_minus_vec(12*12)
    X0 = X0.';
    X_plus_vec = X_plus_vec.';
    X_minus_vec = X_minus_vec.';
end