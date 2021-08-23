function Pxy = average_Pxy_fun(L,lambda,X_k_minus,Y_k_minus,X,Y)
    % weight
    alpha=1e-3;
    beta=2;  
    Wc0 = lambda/(L+lambda)+(1-alpha^2+beta);
    Wc  = 1/(2*(L+lambda));
    Wc  = [Wc0 Wc+zeros(1,2*L)];
            
    % Pxy(12*8)
    Pxy = (X-X_k_minus(:,ones(1,2*L+1)))*diag(Wc)*(Y-Y_k_minus(:,ones(1,2*L+1))).';
end 