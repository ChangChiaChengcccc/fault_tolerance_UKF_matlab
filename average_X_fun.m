function X_k_minus=average_X_fun(L,lambda,X)
    % next state of X0, X_plus_vec, X_minus_vec
    X0 = X(:,1);
    X_plus_vec = X(:,2:13);
    X_minus_vec = X(:,14:25);
    Wm0 = lambda/(L+lambda);
    Wm = 1/(2*(L+lambda));
    X_k_minus = Wm0*X0;
    
    %X_k_minus 25 times (12*1)matrix plus 
    for i =1:L
    X_k_minus = X_k_minus + Wm*X_plus_vec(:,i) + Wm*X_minus_vec(:,i);
    end
end 