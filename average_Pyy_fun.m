function Pyy = average_Pyy_fun(L,lambda,Y_k_minus,Y)    
    % weight
    alpha=1e-3;
    beta=2;  
    Wc0 = lambda/(L+lambda)+(1-alpha^2+beta);
    Wc  = 1/(2*(L+lambda));
    Wc  = [Wc0 Wc+zeros(1,2*L)];
    
    % Pyy(12*12)
    Pyy = (Y-Y_k_minus(:,ones(1,2*L+1)))*diag(Wc)*(Y-Y_k_minus(:,ones(1,2*L+1))).';
end 