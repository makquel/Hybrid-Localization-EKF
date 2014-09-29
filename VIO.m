%% Localization(GPS + Visual Odometry) Sensor Fusion using kalman filter

clc
clear all
close all

%% Load Transformations (Measurements from visual odometry algorithm)

% H_i = importdata('ZERO/transformation.csv');
H_i = importdata('TRES/transformation.csv');
% H_i = importdata('QUATRO/transformation.csv');
%% Load GPS Measurements

% Z_gps = importdata('ZERO/gps.csv');
% Z_gps = importdata('TRES/gps.csv');
Z_gps = importdata('TRES/gps.csv');
% Z_gps = importdata('QUATRO/gpsx.csv');
% Z_gps(:,3) = Z_gps(:,3)*180/pi;
%% Simulation Setup

Nh = length(H_i);
i = 1:4:Nh+1;
x_i = zeros(length(i)-1,3);
R_t = eye(4);
Explored_area = [0 2.3 0 1.75];
N = length(Z_gps);
%%%%%%%%%  storage  %%%%%%%%
InnovStore = NaN*zeros(3,N);
uStore = NaN*zeros(3,N-1);
XStore = NaN*zeros(3,N-1);
VOStore = NaN*zeros(3,N-1);
VOErrorStore = NaN*zeros(3,N-1);
%3*sigma criterion
XErrStore = NaN*zeros(3,N);
PStore = NaN*zeros(3,N-1);
%% Plot Observations GPS vs Visual Odometry

%Map Setup
figure(1);hold on;grid on;
axis(Explored_area); 
xlabel('x (m)');
ylabel('y (m)');

%Initial position and heading
beta = (Z_gps(1,3)-90)*pi/180;
Tgps = [cos(beta) 0 -sin(beta) Z_gps(1,1); ...
        0 1 0 0; ...
        sin(beta) 0 cos(beta) Z_gps(1,2); ...
        0 0 0 1];
P = [0 0.15 0 0 0; ...
     0 0 0 0 0; ...
     0 0 0 0.15 0; ...
     1 1 1 1 1];
T = Tgps;
x = T*P;

%Vision Coordinate System
p0 = [x(1,1) x(3,1)];
p1 = [x(1,2) x(3,2)];
p2 = [x(1,3) x(3,3)];
p3 = [x(1,4) x(3,4)];
vectarrow(p0,p1,'r');hold on
vectarrow(p2,p3,'b');hold on
plot(x(1),x(3),'go');hold on

for observation = 1:length(i)-1

    %GPS
    quiver(Z_gps(observation,1),Z_gps(observation,2),cosd(Z_gps(observation,3)),sind(Z_gps(observation,3)),0.2,'k');hold on
    plot(Z_gps(observation,1),Z_gps(observation,2),'ko');  
    
    %Visual Odometry System
    if(observation > 1)
        R_t = H_i(i(observation):i(observation+1)-1,:);
        R = R_t(1:3,1:3);
        t = R_t(1:3,4);
        aux = rodrigues(R)*180/pi;
        %
        R_t(1,4) = R_t(1,4);
        R_t(3,4) = -1*R_t(3,4);
        %
        T = T*R_t;
        Re = T(1:3,1:3);
        auxe = rodrigues(Re)*180/pi;
        x = T*P;
    
        p0 = [x(1,1) x(3,1)];
        p1 = [x(1,2) x(3,2)];
        p2 = [x(1,3) x(3,3)];
        p3 = [x(1,4) x(3,4)];
        vectarrow(p0,p1,'r');hold on
        vectarrow(p2,p3,'b');hold on
        plot(x(1),x(3),'go');hold on
    
        uStore(1,observation) = t(1);
        uStore(2,observation) = t(3);
        uStore(3,observation) = aux(2);
        
        VOStore(1,observation) = x(1);
        VOStore(2,observation) = x(3);
        VOStore(3,observation) = 270-(auxe(2)+180);
%         Z_gps(observation,3)
%         270-(auxe(2)+180)
%         Z_gps(observation,3)-(270-(auxe(2)+180))
        VOErrorStore(1,observation) = (Z_gps(observation,1) - x(1));
        VOErrorStore(2,observation) = (Z_gps(observation,2) - x(3));
        VOErrorStore(3,observation) = (Z_gps(observation,3) - (270-(auxe(2)+180)));
    end
    
   
%     pause
end
quiver(Z_gps(length(i),1),Z_gps(length(i),2),cosd(Z_gps(length(i),3)),sind(Z_gps(length(i),3)),0.2,'k');hold on
plot(Z_gps(length(i),1),Z_gps(length(i),2),'ko');

hleg1 = legend('Odometria Visual','Ground Truth');
set(hleg1,'Location','NorthWest')
set(hleg1,'Interpreter','none')
% matlab2tikz( 'Visual-Odometry-RPro.tikz' );

%% VO Error - TODO: Create gt vs vo graph for error purpose
% figure(2);hold on;grid on;
% % errorbar(1:19,VOErrorStore(2,:),'b');
% subplot(3,1,1);plot(VOErrorStore(1,:));
% title('Odometry Error');ylabel('x');
% 
% subplot(3,1,2);plot(VOErrorStore(2,:));
% ylabel('y');
% 
% subplot(3,1,3);plot(VOErrorStore(3,:));hold on;
% ylabel('\theta');

%% Localization KF (GPS + VIO)
SNR = 30; % Signal-to-noise ratio dB
%Map Setup
figure(3);hold on;grid on;
axis(Explored_area); 
xlabel('x (m)');
ylabel('y (m)');

% Probabilistic covariace
P = [1 0 0; ...
     0 1 0; ...
     0 0 1];
 
Q = [0.04085 0 0; ...
     0 0.02099 0; ...
     0 0 0.05055];
 
R = [0.0156 0 0; ...
     0 0.0161 0; ...
     0 0 2.0]; %24.0334

u = uStore;
Xs = [Z_gps(1,1); Z_gps(1,2); Z_gps(1,3)];
for k = 1:N-1
    A = [1 0 0; ...
         0 1 0; ...
         0 0 1];
     
    B = [1 0 0; ...
         0 1 0; ...
         0 0 1];
     
    C = [1 0 0; ...
         0 1 0; ...
         0 0 1];
    % Prediction
    Xs = A*Xs + B*u(:,k);
    P = A*P*A' + B*Q*B';

    % Update
    H = [1 0 0; ...
         0 1 0; ...
         0 0 1];
    Z = [Z_gps(k,1); Z_gps(k,2); Z_gps(k,3)];
%     Z = [awgn(Z_gps(k+1,1),SNR,'measured'); ...
%          awgn(Z_gps(k+1,2),SNR,'measured'); ...
%          awgn(Z_gps(k+1,3),SNR,'measured')];
     
    Innov = Z - H*Xs;
    S = H*P*H'+ R;
    K = P*H'*inv(S);
    Xs = Xs + (K*Innov);
    P = (eye(3) - (K*H))*P;
    XStore(:,k) = Xs
    
    PStore(:,k) = sqrt(diag(P));
    XErrStore(1,k) = Z_gps(k,1) - Xs(1);
    XErrStore(2,k) = Z_gps(k,2) - Xs(2);
    XErrStore(3,k) = Z_gps(k,3) - Xs(3);
% %     pause
%     Innov = [NaN; NaN; NaN];
end
% A*[4.5229; 3.5558; 39.3616] + B*u(:,3)
plot(XStore(1,:),XStore(2,:),'r');
plot(Z_gps(:,1),Z_gps(:,2),'k');
plot(VOStore(1,:),VOStore(2,:),'b');
hleg2 = legend('Estimated','GPS(30dB)','Visual Odometry');
set(hleg2,'Location','NorthWest')
set(hleg2,'Interpreter','none')
% matlab2tikz( 'Visual-Odometry-EPro.tikz' );
% % plot(1.8100,0.6921,'*r');hold on;plot(2.6237,2.1451,'*r');            
%% Stats
% min_observation_index = 1;
% figure(4);
% subplot(3,1,1);plot(XErrStore(1,min_observation_index:N));hold on;
% plot(3*PStore(1,min_observation_index:N),'r');plot(-3*PStore(1,min_observation_index:N),'r');
% title('Erro e covariancia');ylabel('x');
% subplot(3,1,2);plot(XErrStore(2,min_observation_index:N));hold on;
% plot(3*PStore(2,min_observation_index:N),'r');plot(-3*PStore(2,min_observation_index:N),'r')
% ylabel('y');
% subplot(3,1,3);plot(XErrStore(2,min_observation_index:N));hold on;
% plot(3*PStore(3,min_observation_index:N),'r');plot(-3*PStore(3,min_observation_index:N),'r')
% ylabel('\theta');
% matlab2tikz( '3sigma-localization.tikz' );

% figure(4);
% % subplot(3,1,1);
% plot(XErrStore(1,:));hold on;
% plot(3*PStore(1,:),'r');plot(-3*PStore(1,:),'r');
% % title('Erro e covariancia');
% ylabel('x');
% % matlab2tikz( 'FVIO-X-res.tex' );
% figure(5);
% % subplot(3,1,2);
% plot(XErrStore(2,:));hold on;
% plot(3*PStore(2,:),'r');plot(-3*PStore(2,:),'r')
% ylabel('y');
% % matlab2tikz( 'FVIO-Y-res.tex' );
% figure(6);
% % subplot(3,1,3);
% plot(XErrStore(3,:));hold on;
% plot(3*PStore(3,:),'r');plot(-3*PStore(3,:),'r')
% ylabel('\theta');
% % matlab2tikz( 'FVIO-T-res.tex' );


%% Plot Observations pt2
%Map Setup
% figure(2);hold on;grid on;
% axis([0 2.3 0 1.75]); 
% xlabel('x (m)');
% ylabel('y (m)');
% 
% %Initial position and heading
% beta = (Z_gps(1,3)-90)*pi/180;
% Tgps = [cos(beta) 0 -sin(beta) Z_gps(1,1); ...
%         0 1 0 0; ...
%         sin(beta) 0 cos(beta) Z_gps(1,2); ...
%         0 0 0 1];
% P = [0 0.15 0 0 0; ...
%      0 0 0 0 0; ...
%      0 0 0 0.15 0; ...
%      1 1 1 1 1];
% T = Tgps;
% x = T*P;
% 
% %Vision Coordinate System
% p0 = [x(1,1) x(3,1)];
% p1 = [x(1,2) x(3,2)];
% p2 = [x(1,3) x(3,3)];
% p3 = [x(1,4) x(3,4)];
% vectarrow(p0,p1,'r');hold on
% vectarrow(p2,p3,'b');hold on
% plot(x(1),x(3),'ko');hold on

for observation = 1:length(i)-1
%     x = x_i(observation,1);
%     y = x_i(observation,2);
%     theta = x_i(observation,3);
%     
%     x_r = [X Y THETA];
%     DrawCar(x_r,0.127,2,'r');
%     
%     line([X,X + x],[Y,Y + y],'color','k');
%     X = X + x;
%     Y = Y + y;
%     THETA = THETA + theta;
    
%     line([Zr(1),Zr(1)],[Zr(2),Zr(2) + Z_i(observation,2)],'color','b');
%     Zr(2) = Zr(2) + Z_i(observation,2);
% %     DrawCar(Zr,'b');
%     Xr = [x_abs y_abs phi_abs];
% %     plot(x_abs,y_abs,'ko');
% %     DrawCar(Xr,'r')
%     line([x_abs,Zr(1)],[y_abs,Zr(2)],'color','r');
	
%     pause
end

% hleg1 = legend('Trajetoria Estimada','Encoder', 'Diferença');%ground truth
% set(hleg1,'Location','NorthWest')
% set(hleg1,'Interpreter','none')

%% Encoder's Covariance(R)
% Z_i = importdata('../../measurements/Sensors.csv');
% Z_i(:,2) = Z_i(:,2)*0.001;
% Z_i(:,3) = Z_i(:,3)*0.001;
% Z_i(:,4) = Z_i(:,4)*0.001;
% Z_i(:,5) = Z_i(:,5)*0.001;
% M = length(Z_i);
% %here goes std
% % Za_1 = mean(Z_i((3:length(Z_i)),2));
% % Za_2 = mean(Z_i((3:length(Z_i)),3));
% % Za_3 = mean(Z_i((3:length(Z_i)),4));
% % Za_4 = mean(Z_i((3:length(Z_i)),5));
% Za_1 = mean(Z_i(:,2));
% Za_2 = mean(Z_i(:,3));
% Za_3 = mean(Z_i(:,4));
% Za_4 = mean(Z_i(:,5));
% R_1 = 0;
% R_2 = 0;
% R_3 = 0;
% R_4 = 0;
% for k = 1:M
%     R_1 = R_1 + (Z_i(k,2) - Za_1)*(Z_i(k,2) - Za_1)';
%     R_2 = R_2 + (Z_i(k,3) - Za_2)*(Z_i(k,3) - Za_2)';
%     R_3 = R_3 + (Z_i(k,4) - Za_3)*(Z_i(k,4) - Za_3)';
%     R_4 = R_4 + (Z_i(k,5) - Za_4)*(Z_i(k,5) - Za_4)';
% %     pause
% end 
% R_1 = R_1/M;
% R_2 = R_2/M;
% R_3 = R_3/M;
% R_4 = R_4/M;
% R_covariance = diag([R_1 R_2 R_3 R_4]);
% sym(R_covariance)
