close all

%% initialize
% initial ctrl_state = [x; v; O; W; R; E] (25*1)
x = [0 0 0]'; 
v = [0 0 0]'; 
O = [0 0 0]';
W = [0 0 0]'; 
E = [1 1 1 1]';
R = [1 0 0 0 1 0 1e-10 1e-10 1]';
ctrl_state = [x; v; O; W; E; R];

% initial lambda, P
L=12;      %number of state
m=8;       %number of measurement state                          
lambda=0;  %scaling factor
P = 0.01*eye(L);

% declare next_ctrl_state, next_ctrl_state_vec
next_ctrl_state = zeros(25,1);
next_ctrl_state_vec = ctrl_state(:,ones(1,25));

% declare next_X_vec, X_k_vec_ukf, Y_sensor 
next_X_vec = zeros(12,25);
X_k_vec_ukf = zeros(12,1000);
Y_sensor = zeros(8,1);

% declare next_ctrl_state_vec_without_ukf, X_k_vec_without_ukf
next_ctrl_state_vec_without_ukf = zeros(25,1000);
X_k_vec_without_ukf = zeros(12,1000);

%% loop
for j = 1:1000
%% start Predict 
% X_sigma_points(12*25)
X0 = ctrl_state2state(ctrl_state); %(25*1->12*1)
[X0,X_plus_vec,X_minus_vec] = sigma_fun(L,lambda,P,X0);
X_sigma_points = [X0 X_plus_vec X_minus_vec];

    for i = 1:25
        % tmp_X(12*1)
        tmp_X = X_sigma_points(:,i);

        % revise ctrl_state by X_sigma_points   
        ctrl_state = state2ctrl_state(tmp_X,ctrl_state); %(12*1->25*1)
        
        % get next_ctrl_state from geo_ctrl_fun
        % ctrl_state = [x; v; O; W; R; E] (25*1)
        next_ctrl_state = geo_ctrl_fun(ctrl_state);
        next_ctrl_state_vec(:,i) = next_ctrl_state;

        % get next_X by next_ctrl_state
        % store to next_X_vec (memorize the state)
        next_X = ctrl_state2state(next_ctrl_state);
        next_X_vec(:,i) = next_X;
    end

    % get x_k_minus by averaging sigma points
    X_k_minus = average_X_fun(L,lambda,next_X_vec);

    % get Pxx with noise Q
    Pxx = average_Pxx_fun(L,lambda,X_k_minus,next_X_vec);

    %% start Correct 
    % get Y_sensor(with noise), Y_k_minus
    Y_sensor(1) = next_X_vec(1,1)+wgn(1, 1, -60);  % z
    Y_sensor(2) = next_X_vec(2,1)+wgn(1, 1, -50);  % z_dot
    Y_sensor([4 6 8]) = next_X_vec([4 6 8],1);     % angular velocity
    Y_sensor([3 5 7]) = next_X_vec([3 5 7],1)+wgn(3, 1, -42); % orientation
    
    Y = next_X_vec(1:8,:); %(8*25)
    Y_k_minus = average_Y_fun(L,lambda,Y); %(8*1)

    % get Pyy(8*8)
    Pyy = average_Pyy_fun(L,lambda,Y_k_minus,Y);

    % get Pxy(12*8)
    Pxy = average_Pxy_fun(L,lambda,X_k_minus,Y_k_minus,next_X_vec,Y);

    % R(noise)
    r=0.01;                 %std of process 
    R_cov=r^2*eye(m);       %covariance of process
    % gain K(12*8)
    K = Pxy*inv(Pyy);       %+R_cov

    % X_k_ukf
    X_k_ukf = X_k_minus + K*(Y_sensor-Y_k_minus); % Y is from sensor

    % P
    %P = Pxx - K*Pyy*K.';

    % update state
    X_k_vec_ukf(:,j) = X_k_ukf;
    ctrl_state = next_ctrl_state_vec(:,1);
    ctrl_state = state2ctrl_state(X_k_ukf,ctrl_state);
    ctrl_state(13:16) = E';
        
    % without ukf
    next_ctrl_state_vec_without_ukf(:,j) = next_ctrl_state_vec(:,1);
    X_k_vec_without_ukf(:,j) = ctrl_state2state(next_ctrl_state_vec_without_ukf(:,j));
end

%% plot position state with/without ukf
figure(1)
title('Position')
subplot(4,1,1)
plot(1:1000,X_k_vec_ukf(1,:),'-','LineWidth',2); hold on;
plot(1:1000,X_k_vec_without_ukf(1,:),'--','LineWidth',2);
title('Z height with iteration')
xlabel('iteration')
ylabel('z height')
legend('model with noise','after ukf')

subplot(4,1,2)
plot(1:1000,X_k_vec_ukf(3,:),'-','LineWidth',2); hold on;
plot(1:1000,X_k_vec_without_ukf(3,:),'--','LineWidth',2);
title('roll with iteration')
xlabel('iteration')
ylabel('roll(radian)')

subplot(4,1,3)
plot(1:1000,X_k_vec_ukf(5,:),'-','LineWidth',2); hold on;
plot(1:1000,X_k_vec_without_ukf(5,:),'--','LineWidth',2);
title('pitch with iteration')
xlabel('iteration')
ylabel('pitch(radian)')

subplot(4,1,4)
plot(1:1000,X_k_vec_ukf(7,:),'-','LineWidth',2); hold on;
plot(1:1000,X_k_vec_without_ukf(7,:),'--','LineWidth',2);
title('yaw with iteration')
xlabel('iteration')
ylabel('yaw(radian)')

%% plot velocity state with/without ukf
figure(2)
title('Velocity')
subplot(4,1,1)
plot(1:1000,X_k_vec_ukf(2,:),'-','LineWidth',2); hold on;
plot(1:1000,X_k_vec_without_ukf(2,:),'--','LineWidth',2);
title('Z velocity with iteration')
xlabel('iteration')
ylabel('Z velocity(m/s)')
legend('model with noise','after ukf')

subplot(4,1,2)
plot(1:1000,X_k_vec_ukf(4,:),'-','LineWidth',2); hold on;
plot(1:1000,X_k_vec_without_ukf(4,:),'--','LineWidth',2);
title('roll velocity with iteration')
xlabel('iteration')
ylabel('roll_velocity(radian/s)')

subplot(4,1,3)
plot(1:1000,X_k_vec_ukf(6,:),'-','LineWidth',2); hold on;
plot(1:1000,X_k_vec_without_ukf(6,:),'--','LineWidth',2);
title('pitch velocity with iteration')
xlabel('iteration')
ylabel('pitch_velocity(radian/s)')

subplot(4,1,4)
plot(1:1000,X_k_vec_ukf(8,:),'-','LineWidth',2); hold on;
plot(1:1000,X_k_vec_without_ukf(8,:),'--','LineWidth',2);
title('yaw velocity with iteration')
xlabel('iteration')
ylabel('yaw_velocity(radian/s)')

%% plot e1~e4 with/without ukf
figure(3)
title('Efficiency')
subplot(4,1,1)
plot(1:1000,X_k_vec_ukf(9,:),'-','LineWidth',1.3); hold on;
plot(1:1000,X_k_vec_without_ukf(9,:),'--','LineWidth',1.3);
title('e1 percentage with iteration')
xlabel('iteration')
ylabel('e1 percentage')
axis([0 1000 -5 5])
legend('model with noise','after ukf')

subplot(4,1,2)
plot(1:1000,X_k_vec_ukf(10,:),'-','LineWidth',1.3); hold on;
plot(1:1000,X_k_vec_without_ukf(10,:),'--','LineWidth',1.3);
title('e2 percentage with iteration')
xlabel('iteration')
ylabel('e2 percentage')
axis([0 1000 -5 5])

subplot(4,1,3)
plot(1:1000,X_k_vec_ukf(11,:),'-','LineWidth',1.3); hold on;
plot(1:1000,X_k_vec_without_ukf(11,:),'--','LineWidth',1.3);
title('e3 percentage with iteration')
xlabel('iteration')
ylabel('e3 percentage')
axis([0 1000 -5 5])

subplot(4,1,4)
plot(1:1000,X_k_vec_ukf(12,:),'-','LineWidth',1.3); hold on;
plot(1:1000,X_k_vec_without_ukf(12,:),'--','LineWidth',1.3);
title('e4 percentage with iteration')
xlabel('iteration')
ylabel('e4 percentage')
axis([0 1000 -5 5])

%% test
% position
figure(4)
title('Position')
subplot(4,1,1)
plot(1:30,X_k_vec_ukf(1,1:30),'-'); hold on;
plot(1:30,X_k_vec_without_ukf(1,1:30),'--');
title('Z height with iteration')
xlabel('iteration')
ylabel('z height')
legend('model with noise','after ukf')

subplot(4,1,2)
plot(1:30,X_k_vec_ukf(3,1:30),'-'); hold on;
plot(1:30,X_k_vec_without_ukf(3,1:30),'--');
title('roll with iteration')
xlabel('iteration')
ylabel('roll(radian)')

subplot(4,1,3)
plot(1:30,X_k_vec_ukf(5,1:30),'-'); hold on;
plot(1:30,X_k_vec_without_ukf(5,1:30),'--');
title('pitch with iteration')
xlabel('iteration')
ylabel('pitch(radian)')

subplot(4,1,4)
plot(1:30,X_k_vec_ukf(7,1:30),'-'); hold on;
plot(1:30,X_k_vec_without_ukf(7,1:30),'--');
title('yaw with iteration')
xlabel('iteration')
ylabel('yaw(radian)')

% velocity
figure(5)
title('Velocity')
subplot(4,1,1)
plot(1:30,X_k_vec_ukf(2,1:30),'-','LineWidth',2); hold on;
plot(1:30,X_k_vec_without_ukf(2,1:30),'--','LineWidth',2);
title('Z velocity with iteration')
xlabel('iteration')
ylabel('Z velocity(m/s)')
legend('model with noise','after ukf')

subplot(4,1,2)
plot(1:30,X_k_vec_ukf(4,1:30),'-','LineWidth',2); hold on;
plot(1:30,X_k_vec_without_ukf(4,1:30),'--','LineWidth',2);
title('roll velocity with iteration')
xlabel('iteration')
ylabel('roll_velocity(radian/s)')

subplot(4,1,3)
plot(1:30,X_k_vec_ukf(6,1:30),'-','LineWidth',2); hold on;
plot(1:30,X_k_vec_without_ukf(6,1:30),'--','LineWidth',2);
title('pitch velocity with iteration')
xlabel('iteration')
ylabel('pitch_velocity(radian/s)')

subplot(4,1,4)
plot(1:30,X_k_vec_ukf(8,1:30),'-','LineWidth',2); hold on;
plot(1:30,X_k_vec_without_ukf(8,1:30),'--','LineWidth',2);
title('yaw velocity with iteration')
xlabel('iteration')
ylabel('yaw_velocity(radian/s)')

% efficiency
figure(6)
title('Efficiency')
subplot(4,1,1)
plot(1:30,X_k_vec_ukf(9,1:30),'-','LineWidth',1.3); hold on;
plot(1:30,X_k_vec_without_ukf(9,1:30),'--','LineWidth',1.3);
title('e1 percentage with iteration')
xlabel('iteration')
ylabel('e1 percentage')
axis([0 30 -5 5])
legend('model with noise','after ukf')

subplot(4,1,2)
plot(1:30,X_k_vec_ukf(10,1:30),'-','LineWidth',1.3); hold on;
plot(1:30,X_k_vec_without_ukf(10,1:30),'--','LineWidth',1.3);
title('e2 percentage with iteration')
xlabel('iteration')
ylabel('e2 percentage')
axis([0 30 -5 5])

subplot(4,1,3)
plot(1:30,X_k_vec_ukf(11,1:30),'-','LineWidth',1.3); hold on;
plot(1:30,X_k_vec_without_ukf(11,1:30),'--','LineWidth',1.3);
title('e3 percentage with iteration')
xlabel('iteration')
ylabel('e3 percentage')
axis([0 30 -5 5])

subplot(4,1,4)
plot(1:30,X_k_vec_ukf(12,1:30),'-','LineWidth',1.3); hold on;
plot(1:30,X_k_vec_without_ukf(12,1:30),'--','LineWidth',1.3);
title('e4 percentage with iteration')
xlabel('iteration')
ylabel('e4 percentage')
axis([0 30 -5 5])