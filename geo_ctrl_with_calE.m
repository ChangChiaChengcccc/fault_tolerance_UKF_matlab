clear all 
clc

% from main function
ma = 1.15;
time = 10;
dt = 0.001;
iteration = time/dt;

% controller elements 
% state x, v, W, R, xd_dot_dot, Wd_dot, 
% xd, b1d, O and desired state
x = [0 0 0]'; xd = [0 0 1]';
v = [0 0 0]'; vd = [0 0 0]';
O = [0 0 0]';
W = [0 0 0]'; Wd = [0 0 0]';
R = [1 0 1e-10;
     0 1 1e-10;
     0 0 1 ];
b1d = [1 0 0]';

% E,E_diag
E = [0.7 1 0.7 1]';
E_diag = diag(E);

% Rd = eye(3);
xd_dot_dot = [0 0 0]';
Wd_dot = [0 0 0]';


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
buffer = zeros(25,iteration);
bufferR = zeros(3,3);

%state_E_vec
state_E_vec = zeros(4,iteration);

for i = 1:iteration
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
    state_E_vec(:,i) = state_E;
      
    % for plot and pass the value
    buffer(1:3,i) = x;
    buffer(4:6,i) = v;
    buffer(7:9,i) = O;
    buffer(10:12,i) = W;
    buffer(13:21,i) = reshape(R,9,1);
    buffer(22:25,i) = state_E;
    
    % for check X value
    X = [x(3); v(3); O(1); W(1); O(2); W(2); O(3); W(3); state_E];
end
figure(1)
% plot buffer
time =  linspace(1,iteration,iteration);
for j = 1:12
subplot(4,3,j);
plot(time,buffer(j,:));
end
figure(2)
subplot(4,1,1)
plot(1:iteration,buffer(22,:))
axis([0 iteration -1.5 1.5])

subplot(4,1,2)
plot(1:iteration,buffer(23,:))
axis([0 iteration -1.5 1.5])

subplot(4,1,3)
plot(1:iteration,buffer(24,:))
axis([0 iteration -1.5 1.5])

subplot(4,1,4)
plot(1:iteration,buffer(25,:))
axis([0 iteration -1.5 1.5])