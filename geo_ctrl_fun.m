function buffer=geo_ctrl_fun(ctrl_state)

% convert parameter
x = ctrl_state(1:3);
v = ctrl_state(4:6);
O = ctrl_state(7:9);
W = ctrl_state(10:12);
E = ctrl_state(13:16);
R = ctrl_state(17:25);

% from main function
ma = 1.15;
dt = 0.01;

% reshape R
R = reshape(R,3,3);

% controller elements 
% state x, v, W, R, xd_dot_dot, Wd_dot, 
% xd, b1d, O and desired state
vd = [0 0 0]';
Wd = [0 0 0]';
xd = [0 0 1]';
b1d = [1 0 0]';
% Rd = eye(3);
xd_dot_dot = [0 0 0]';
Wd_dot = [0 0 0]';


% E,E_diag
E_diag = diag(E);

% allocation matrix
d = 0.225; 
c_tauf = 1.347e-2;
L = [1          1          1       1; ...
     0         -d          0       d; ...
     d          0         -d       0; ...
     -c_tauf    c_tauf    -c_tauf  c_tauf];
L_inv = inv(L);

% other coefficients
g = 9.81; e_3 = [0 0 1]';
J = [0.0131 0 0; 0 0.0131 0; 0 0 0.0244];

% gain kx, kv, kW, kR 
kx = 2*ma;  
kv = 1.5*ma; 
kR = 8.81;
kW = 2.54;

%for plot
buffer = zeros(25,1);
bufferR = zeros(3,3);

for i = 1:1
    % error
    ex = x-xd;
    ev = v-vd;

    % get b3_d & update Rd
    tmp =  (-kx*ex - kv*ev + ma*xd_dot_dot - ma*g*e_3 );
    b3d = -tmp / norm(tmp);
    b2d = cross(b3d,b1d) / norm(cross(b3d,b1d));
    b1d_proj = cross(b2d, b3d) / norm(cross(b2d, b3d)) ;
    Rd = [b1d_proj b2d b3d];

    %eR, eW
    eR = 0.5*vee((Rd'*R -R'*Rd));
    eW = W - R'*Rd*Wd;

    % controller(f,M) 
    f =  (kx*ex + kv*ev + ma*g*e_3 - ma*xd_dot_dot)'*(R*e_3);
    M = - kR*eR - kW*eW + cross(W, J*W) ...
        - J*(hat(W)*R'*Rd*Wd - R'*Rd*Wd_dot);
    
    % controller(ff,Mf)
    U  = [f;M];
    F  = inv(L)*U;
    Ff = E_diag*F;
    Uf = L*Ff;
    ff = Uf(1);
    Mf = Uf(2:4);
    
    % The geometric equation of motion 
    %x_dot = v;
    v_dot = g*e_3 - ff*R*e_3/ma;
    R_dot = R*hat(W);
    W_dot = J\(Mf - cross(W, J*W));

    %update v, x, W, R, O
    v = v + v_dot*dt;
    x = x + v*dt;
    W = W + W_dot*dt;
    R = R + R_dot*dt;
    O = O + W*dt;
    
    % get the E elements from state 
    % compare elements
    tmp1 = R*e_3;
    tmp2 = ma*g*e_3 - ma*v_dot;
    % average ff
    state_ff = tmp2./tmp1;
    state_ff = state_ff(3);%mean(state_ff,'all');
    % Mf
    state_Mf = J*W_dot + cross(W,J*W);
    % state_Uf
    state_Uf = [state_ff;state_Mf];
    tmp1 = L_inv*state_Uf;
    tmp2 = F;
    state_E = tmp1./tmp2;

    % for plot and pass the value
    buffer(1:3,i) = x;
    buffer(4:6,i) = v;
    buffer(7:9,i) = O;
    buffer(10:12,i) = W;
    buffer(13:16,i) = state_E;
    buffer(17:25,i) = reshape(R,9,1);

end
end

% % plot buffer
% time =  linspace(1,iteration,iteration);
% for j = 1:16
% subplot(6,3,j)
% plot(time,buffer(j,:),'--rs','LineWidth',1,...
%                        'MarkerEdgeColor','b',...
%                        'MarkerFaceColor','b',...
%                        'MarkerSize',1)
% end