function y_k_minus=average_Y_fun(L,lambda,Y)
    % next state of Y0, Y_plus_vec, Y_minus_vec
    Y0 = Y(:,1);
    Y_plus_vec = Y(:,2:13);
    Y_minus_vec = Y(:,14:25);
    Wm0 = lambda/(L+lambda);
    Wm = 1/(2*(L+lambda));
    y_k_minus = Wm0*Y0;
    
    %x_k_minus 25 times (12*1)matrix plus 
    for i =1:L
    y_k_minus = y_k_minus + Wm*Y_plus_vec(:,i) + Wm*Y_minus_vec(:,i);
    end
end 