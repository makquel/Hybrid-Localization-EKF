%% Odometric(Dead Reckoning) Sensor Fusion using kalman filter

clear all
close all
clc

%% System parameters

L = 0.255; % Length of longitudinal axis of vehicle
W = 0.16; % Width of vehicle
dT = 0.252; % Sampling time

%% Load Wheels ands Steering Measurements
% Z_i = importdata('UNO/Sensors.csv');
Z_i = importdata('DOS/Sensors.csv');
% Z_i = importdata('TRES/Sensors.csv');
% Z_i = importdata('QUATRO/Sensors.csv');
Z_i(:,1) = AWrap(Z_i(:,1).*0.1); %awgn(10*ones(1,length(Z_i)),20,'measured'); % Omega <Noisy input>
Z_i(:,2) = Z_i(:,2)*0.001;% FL
Z_i(:,3) = Z_i(:,3)*0.001;% FR
Z_i(:,4) = Z_i(:,4)*0.001;% RL
Z_i(:,5) = Z_i(:,5)*0.001;% RR

%% Load GPS Measurements
% Z_gps = importdata('UNO/gps.csv');
Z_gps = importdata('DOS/gps.csv');
% Z_gps = importdata('TRES/gps.csv');
% Z_gps = importdata('QUATRO/gps.csv');

%% Average Speed
Speed = mean(Z_i(:,2)*0.001)/dT;

%% Simulation Setup
N = length(Z_gps);% length(Z_i);% Length of recorded data

% Standard deviation of variables measured
SigmaR1 = 0.005;%deve ser maior da resolucao do encoder
SigmaR2 = 0.005;
SigmaR3 = 0.005;
SigmaR4 = 0.005;
SigmaSt = 0.03;%0.001
% Covariance Matrices
P = 1*diag([0.1^2 0.1^2]);
Q = 1*diag([0.1^2 0.1^2]);
R = diag([SigmaSt^2 SigmaR1^2 SigmaR2^2 SigmaR3^2 SigmaR4^2]);

% Initial guess (prior @ k = 0)
steeri=10*pi/180;
delta=0.1;
x = [delta; atan(L*steeri/delta)];% \delta \omega mean(Z_i(3,2:5)) 


SNR = 20; % Signal-to-noise ratio dB

XErrStore = NaN*zeros(2,N);
PStore = NaN*zeros(2,N);
XStore = NaN*zeros(2,N);
XStore_odometry = NaN*zeros(3,N);
%% Odometric EKF                    NOTE: omega = omega_l = omega_r
figure(1);axis([0 2.3 0 1.75]);hold on;grid on;
min_observation_index = 10;
xtrue = [Z_gps(min_observation_index,1); ...
         Z_gps(min_observation_index,2); ...
         Z_gps(min_observation_index,3)*pi/180];  
     
for k = min_observation_index:N
    
    H = [-L*x(2)/(x(1)^2) L/x(1); ...
        1 -(W/2); ...
        1 (W/2); ...
        1 -(W/2); ...
        1 (W/2)];
  % Prediction
    x = x; % + awgn(x,SNR,'measured'); % Add white Gaussian noise 
    P = P + Q;     
  % Update  
    Z = [tand(Z_i(k,1)); Z_i(k,4); Z_i(k,5); Z_i(k,2)*cosd(Z_i(k,1)); Z_i(k,3)*cosd(Z_i(k,1))]; 
    ze = H*x;
    ze(1) = L*x(2)/x(1);
    y = Z - ze;
    S = H*P*H'+ R;
    K = P*H'*inv(S);
    x = x + (K*y);
    P = (eye(2) - (K*H))*P;
    PStore(:,k) = sqrt(diag(P));
    XStore(:,k) = x;
    XErrStore(:,k) = x - [(Z_i(k,4)+Z_i(k,5))/2; Z_i(k,1)*pi/180];
    
    u = [XStore(1,k)/dT;-atan(L*XStore(2,k)./XStore(1,k))];% WARNING: Angle correction
    xtrue = AckermannModel(xtrue,u,dT,L);
    XStore_odometry(:,k) = xtrue;
    if(mod(k-1,4)==0)
        DrawSimpleCar(xtrue,u(2),'r');
    end
    
    if(k > min_observation_index)
        line([Z_gps(k-1,1),Z_gps(k,1)],[Z_gps(k-1,2),Z_gps(k,2)],'color','b');
        line([XStore_odometry(1,k-1),XStore_odometry(1,k)],[XStore_odometry(2,k-1),XStore_odometry(2,k)],'color','r');
    end;
    
    ylabel('y');
    xlabel('x');
 
%     pause(dT);
   
end
% matlab2tikz( 'odometry-animation.tikz' );
%% Delta
figure(2);
DStore = NaN*zeros(1,N);
for k = 1:N
     if(k > min_observation_index)
        DStore(1,k) = sqrt((Z_gps(k-1,1)-Z_gps(k,1))^2+(Z_gps(k-1,2)-Z_gps(k,2))^2);
    end;
end
% DStore(1,1) = 0;
% DStore(1,2) = 0;
plot(XStore(1,min_observation_index:N),'k');hold on; grid on;
%plot(Z_i(min_observation_index:N,2),'-.r');plot(Z_i(min_observation_index:N,3),'-.g');plot(Z_i(min_observation_index:N,4),'-.y');plot(Z_i(min_observation_index:N,5),'-.b');
% plot((Z_i(min_observation_index:N,4)+Z_i(min_observation_index:N,5))/2,'r');
plot(DStore(1,min_observation_index:N),'r');
ylabel('\Delta(k)');
xlabel('k');
% legend('\Delta','FL','FR','RR','RL');
legend('\Delta','\Delta real');
matlab2tikz( 'DR-deltaresponse.tikz' );
%% Omega 
% OmegaEst = atan(L*XStore(2,min_observation_index:N)./XStore(1,min_observation_index:N))*180/pi;
% figure(3);
% % plot(atan(L*XStore(2,:)./XStore(1,:))*180/pi,'k');hold on; grid on;
% plot(OmegaEst,'k');hold on; grid on;
% plot(Z_i(min_observation_index:N,1),'r');
% ylabel('\omega(k)');
% xlabel('k');
% legend('Estimado','\omega');
% % matlab2tikz( 'DR-omegaresponse.tikz' );
%% Error

% figure(4);
% subplot(2,1,1);hold on;
% plot(XErrStore(1,min_observation_index:N),'b');
% plot(3*PStore(1,min_observation_index:N),'r');
% plot(-3*PStore(1,min_observation_index:N),'r');
% ylabel('\Delta');
% title('Erro e covariancia');
% subplot(2,1,2);hold on;
% plot(XErrStore(1,min_observation_index:N),'b');
% plot(3*PStore(2,min_observation_index:N),'r');
% plot(-3*PStore(2,min_observation_index:N),'r');
% ylabel('\omega');
% % matlab2tikz( '3sigma-odometry.tikz' );
%% Plot Results

% figure(5);hold on;axis([0 2.3 0 1.75]);grid on;
% x_o = Z_gps(min_observation_index:N,1);
% y_o = Z_gps(min_observation_index:N,2);
% theta_o = Z_gps(min_observation_index:N,3);
% xtrue = [Z_gps(min_observation_index,1); Z_gps(min_observation_index,2);
% Z_gps(min_observation_index,3)*pi/180];
% 
% % xtrue = [Z_gps(k,1); Z_gps(k,2); Z_gps(k,3)*pi/180];
% % for(k = min_observation_index:N)
% %     plot(XStore_odometry(1,k),XStore_odometry(2,k),'r');
% %     plot(Z_gps(k,1),Z_gps(k,2),'b');
% %     if(k > min_observation_index)
% %         line([Z_gps(k-1,1),Z_gps(k,1)],[Z_gps(k-1,2),Z_gps(k,2)],'color','b');
% %         line([XStore_odometry(1,k-1),XStore_odometry(1,k)],[XStore_odometry(2,k-1),XStore_odometry(2,k)],'color','r');
% %     end;
% %     if(k==N)
% %         DrawCar(XStore_odometry(:,k),u(2),'r');
% %     end;  
% % 
% % end
% plot(XStore_odometry(1,min_observation_index:N),XStore_odometry(2,min_observation_index:N),'r');hold on;
% plot(Z_gps(min_observation_index:N,1),Z_gps(min_observation_index:N,2),'b');
% DrawCar(XStore_odometry(:,N),u(2),'r');
% ylabel('y');
% xlabel('x');
% hleg1 = legend('Odometria convencional','Trajetoria real');
% set(hleg1,'Location','NorthWest')
% set(hleg1,'Interpreter','none')
% % matlab2tikz( 'Dead-Reckoning.tikz' );

%% Animation
fig_pos = [0 0 533 400]; % position and size of the figure window
fillplot_ax_pos = [30 30 500 370];        % position and size of fill plot

fig_col = [1 1 1];  % figure background color
text_col = [0 0 0]; % text color


movieflag = 1;
moviefilename = 'dead-reckoning.avi';
 
% only if our flag is set to 1, will we open up a movie object:
if movieflag == 1
    aviobj = avifile(moviefilename, 'fps', 2, 'compression', 'none');
end

fh= figure('color', fig_col, 'name', 'Tutorial animation movie','Position', fig_pos);
grid on;

tic;
fprintf('\nWe are entering the loop...\n');
% ***************  START THE BIG LOOP AFTER CREATING FIG. ***************
for k = min_observation_index:N
 
%     titlestr = sprintf('GPS(%03d)', k); % our title will change
%     fprintf('Frame %d\n', k); % print loop counter
 
 
    % ***************  1) GPS  ***************
    fillplot_ax = axes;
    set(fillplot_ax, 'Units', 'pixels', 'Position', fillplot_ax_pos);
    u = [XStore(1,k)/dT;-atan(L*XStore(2,k)./XStore(1,k))];% WARNING: Angle correction
    
    hold on;   
    plot(Z_gps(min_observation_index:k,1), Z_gps(min_observation_index:k,2),'-.b');
    plot(XStore_odometry(1,min_observation_index:k), XStore_odometry(2,min_observation_index:k),'r');    
    DrawSimpleCar(XStore_odometry(:,k),u(2),'k');
    hold off;
    
    axis([-0 2.6 -0 1.75]);
    xlabel('x'); ylabel('y');
     
    % ***************  FINISH THE FRAME (AVI if selected)  ***************
    pause(dT);
    if movieflag == 1
        frame = getframe(gcf);           % capture current figure
        aviobj = addframe(aviobj,frame); % use addframe to append frame
    end
    if k < N
        clf; % clear figure except for very last frame
    end
 
end % of the big loop
fprintf('\nDone looping...\n');
 
if movieflag == 1
    fprintf('Saving movie...\n\n');
    aviobj = close(aviobj);
end
toc;
