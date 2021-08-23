function Pxx = average_Pxx_fun(L,lambda,X_k_minus,X)
    % Q
    q=0.01;          %std of process 
    Q=q^2*eye(L);   % covariance of process
    
    % weight
    alpha=1e-3;
    beta=2;  
    Wc0 = lambda/((L+lambda)+(1-alpha^2+beta));
    Wc  = 1/(2*(L+lambda));
    Wc  = [Wc0 Wc+zeros(1,2*L)];
         
    % Pxx(12*12)
    Pxx = (X-X_k_minus(:,ones(1,2*L+1)))*diag(Wc)*(X-X_k_minus(:,ones(1,2*L+1))).';
    
    % process noise
    Pxx = Pxx + Q;
end 