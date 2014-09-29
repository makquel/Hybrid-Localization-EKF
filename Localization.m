%% Localization(GPS + Encoder + Steering Sensor) Sensor Fusion using kalman filter

clear all
close all
clc

%% System parameters

L = 0.255; % Length of vehicle
W = 0.18; % Width of vehicle
dT = 0.252; % Sampling time

%% Load Wheels ands Steering Measurements
% Z_i = importdata('UNO/Sensors.csv');
% Z_i = importdata('DOS/Sensors.csv');
Z_i = importdata('TRES/Sensors.csv');
% Z_i = importdata('QUATRO/Sensors.csv');
Z_i(:,1) = AWrap(Z_i(:,1).*0.1); %awgn(10*ones(1,length(Z_i)),20,'measured'); % Omega <Noisy input>
Z_i(:,2) = Z_i(:,2)*0.001;% FL
Z_i(:,3) = Z_i(:,3)*0.001;% FR
Z_i(:,4) = Z_i(:,4)*0.001;% RL
Z_i(:,5) = Z_i(:,5)*0.001;% RR

%% Load GPS Measurements
% Z_gps = importdata('UNO/gps.csv');
% Z_gps = importdata('DOS/gps.csv'); % min = 8
Z_gps = importdata('TRES/gps.csv');
% Z_gps = importdata('QUATRO/gps.csv');

%% Average Speed
Speed = mean(Z_i(:,2)*0.001)/dT;

%% Simulation Setup
N = length(Z_gps);
SNR = 20; % Signal-to-noise ratio dB for Noisy GPS
min_observation_index = 2;
% Store Odometry results
XErrStore_o = NaN*zeros(2,N);
PStore_o = NaN*zeros(2,N);
XStore_o = NaN*zeros(2,N);
% Store Localization results
XErrStore_l = NaN*zeros(6,N);
PStore_l = NaN*zeros(6,N);
XStore_l = NaN*zeros(6,N);

%% Initial knowledge (prior @ k = 0)

% Standard deviation of controle input vector
SigmaV = 0.1; % 3cm/s std on speed
SigmaPhi = 4*pi/180; % steer inaccuracy
% Standard deviation of odometric variables measured
SigmaR1 = 0.005;
SigmaR2 = 0.005;
SigmaR3 = 0.005;
SigmaR4 = 0.005;
SigmaSt = 0.03;
% Covariance Matrices - Localization filter
P_l = 0.2*eye(6);
Q_l = diag([SigmaV^2 SigmaPhi^2]);
R_l = diag([0.1 0.1 0.1]);
% Covariance Matrices - Odometry filter
P_o = diag([0.1^2 0.1^2]);
Q_o = diag([0.1^2 0.1^2]);
R_o = diag([SigmaSt^2 SigmaR1^2 SigmaR2^2 SigmaR3^2 SigmaR4^2]);

% % Initial guess (prior @ k = 0)
steering = 10*pi/180;
delta = 0.1;
x_o = [delta; atan(L*steering/delta)];% \delta \omega mean(Z_i(3,2:5)a
x_l = [Z_gps(min_observation_index,1); Z_gps(min_observation_index,2); Z_gps(min_observation_index,3)*pi/180; 0; 0; 0];
xtrue = x_l;



%% Localization EKF
figure(1);hold on;grid on;axis([0 2.3 0 1.75]);
xlabel('x(m)');ylabel('y(m)');


for k = min_observation_index:N
    
    H_o = [-L*x_o(2)/(x_o(1)^2) L/x_o(1); ...
        1 -(W/2); ...
        1 (W/2); ...
        1 -(W/2); ...
        1 (W/2)];
  % Prediction
    x_o = x_o;% + awgn(x,SNR,'measured'); % Add white Gaussian noise 
    P_o = P_o + Q_o;
      
  % Update  
    Z_o = [tand(Z_i(k,1)); Z_i(k,4); Z_i(k,5); Z_i(k,2)*cosd(Z_i(k,1)); Z_i(k,3)*cosd(Z_i(k,1))]; 
    ze = H_o*x_o;
    ze(1) = L*x_o(2)/x_o(1);
    y_o = Z_o - ze;
    S_o = H_o*P_o*H_o'+ R_o;
    K_o = P_o*H_o'*inv(S_o);
    x_o = x_o + (K_o*y_o);
    P_o = (eye(2) - (K_o*H_o))*P_o;
    PStore_o(:,k) = sqrt(diag(P_o));
    XStore_o(:,k) = x_o;
    

    u = [XStore_o(1,k)/dT;-atan(L*XStore_o(2,k)./XStore_o(1,k))]; % WARNING: Angle correction 
    Q_l = Q_o;
    
    F = [1 0 -dT*u(1)*sin(x_l(3)) 0 0 0; ... 
         0 1 dT*u(1)*cos(x_l(3)) 0 0 0; ...
         0 0 1 0 0 0; ...
         0 0 0 0 0 0; ...
         0 0 0 0 0 0; ...
         0 0 0 0 0 0];   
    G = [dT*cos(x_l(3)) 0; ... 
         dT*sin(x_l(3)) 0; ...
         dT*tan(u(2))/L dT*u(1)*sec(u(2))^2; ...
         0 0; ...
         0 0;
         0 0];
    H_l = [1 0 0 0 1 0; ... % x + x_biased
           0 1 0 0 0 1; ... % y + y_biased
           0 0 1 1 0 0];    % theta + theta_biased

   % Prediction
    xx = AckermannModel(x_l(1:3),u,dT,L);
    x_l = [xx ; x_l(4:6)];
    P_l = F*P_l*F' + G*Q_l*G';
    
    xtrue = AckermannModel(xtrue,u,dT,L);  
    
   % Update  
    Z_l = [Z_gps(k,1); Z_gps(k,2); Z_gps(k,3)*pi/180];
    y_l = Z_l - H_l*x_l;
    S_l = H_l*P_l*H_l'+ R_l;
    K_l = P_l*H_l'*inv(S_l);
    x_l = x_l + (K_l*y_l);
    P_l = (eye(6) - (K_l*H_l))*P_l;
    PStore_l(:,k) = sqrt(diag(P_l));
    XStore_l(:,k) = x_l;
    
  
    XErrStore_l(1,k) = Z_gps(k,1)-x_l(1);
    XErrStore_l(2,k) = Z_gps(k,2)-x_l(2);
    XErrStore_l(3,k) = (Z_gps(k,1)*pi/180)-x_l(3);
  
    if(mod(k-1,4)==0)
        PlotEllipse(x_l,P_l,1);
        DrawCar(x_l,u(2),'r');
    end
    
    if(k > min_observation_index)
        line([Z_gps(k-1,1),Z_gps(k,1)],[Z_gps(k-1,2),Z_gps(k,2)],'color','k');
        line([XStore_l(1,k-1),XStore_l(1,k)],[XStore_l(2,k-1),XStore_l(2,k)],'color','r');
    end;
    pause(dT);

    
end
hleg1 = legend('Estado do Veiculo','Observacao do GPS','Covariancia');
set(hleg1,'Location','NorthWest')
set(hleg1,'Interpreter','none')

% matlab2tikz( 'localization-animation.tikz' );

%% Graphics: Odometric Filter
% figure(2);
% plot(XStore_o(1,min_observation_index:N),'k');hold on; grid on;
% plot(Z_i(min_observation_index:N,2),'r');plot(Z_i(min_observation_index:N,3),'g');plot(Z_i(min_observation_index:N,4),'y');plot(Z_i(min_observation_index:N,5),'b');
% ylabel('\Delta(k)');
% xlabel('k');
% legend('\Delta','FL','FR','RL','RR');
% 
% figure(3);
% plot(atan(L*XStore_o(2,min_observation_index:N)./XStore_o(1,min_observation_index:N))*180/pi,'k');hold on; grid on;
% plot(Z_i(min_observation_index:N,1),'r');
% ylabel('\omega(k)');
% xlabel('k');
% legend('\omega','\phi');

%% Graphics: Localization Filter
% figure(4);hold on;grid on;axis([0 2.3 0 1.75]);
% xlabel('x(m)');ylabel('y(m)');
% % for(k = min_observation_index:N)
% %     if(k > min_observation_index)
% %         line([XStore_l(1,k-1),XStore_l(1,k)],[XStore_l(2,k-1),XStore_l(2,k)],'color','r');
% %         line([Z_gps(k-1,1),Z_gps(k,1)],[Z_gps(k-1,2),Z_gps(k,2)],'color','k');
% %     end;
% %     if(k==N)
% %         DrawCar(XStore_l(1:3,k),0,'r');
% %     end; 
% % 
% % end;
% plot(XStore_l(1,min_observation_index:N),XStore_l(2,min_observation_index:N),'r');hold on;
% plot(Z_gps(min_observation_index:N,1),Z_gps(min_observation_index:N,2),'b');
% DrawCar(XStore_l(:,N),u(2),'r');
% hleg2 = legend('Posicao Estimada','Trajetoria real');
% set(hleg2,'Location','NorthWest')
% set(hleg2,'Interpreter','none')
% % matlab2tikz( 'Localization.tikz' );
%% Animation
fig_pos = [0 0 533 400]; % position and size of the figure window
fillplot_ax_pos = [30 30 500 370];        % position and size of fill plot

fig_col = [1 1 1];  % figure background color
text_col = [0 0 0]; % text color


movieflag = 1;
moviefilename = 'DR-localization.avi';
 
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
    u = [XStore_o(1,k)/dT;-atan(L*XStore_o(2,k)./XStore_o(1,k))];% WARNING: Angle correction
    X = XStore_l(:,k);
    P = [PStore_l(1,k) 0; 0 PStore_l(2,k)];
    hold on;   
    plot(Z_gps(min_observation_index:k,1), Z_gps(min_observation_index:k,2),'-.b');
    plot(XStore_l(1,min_observation_index:k), XStore_l(2,min_observation_index:k),'r');  
    PlotEllipse(X,P,1);
    DrawCar(XStore_l(:,k),u(2),'k');
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
%% Criterion
% figure(5);
% subplot(3,1,1);plot(XErrStore_l(1,min_observation_index:N));hold on;
% plot(3*PStore_l(1,min_observation_index:N),'r');plot(-3*PStore_l(1,min_observation_index:N),'r');
% title('Erro e covariancia');ylabel('x');
% subplot(3,1,2);plot(XErrStore_l(2,min_observation_index:N));hold on;
% plot(3*PStore_l(2,min_observation_index:N),'r');plot(-3*PStore_l(2,min_observation_index:N),'r')
% ylabel('y');
% subplot(3,1,3);plot(XErrStore_l(2,min_observation_index:N));hold on;
% plot(3*PStore_l(3,min_observation_index:N),'r');plot(-3*PStore_l(3,min_observation_index:N),'r')
% ylabel('\theta');
% % matlab2tikz( '3sigma-localization.tikz' );

% figure(5);
% % subplot(3,1,1);
% plot(XErrStore_l(1,:));hold on;
% plot(3*PStore_l(1,:),'r');plot(-3*PStore_l(1,:),'r');
% % title('Erro e covariancia');
% ylabel('x');
% matlab2tikz( 'FX-res.tex' );
% figure(6);
% % subplot(3,1,2);
% plot(XErrStore_l(2,:));hold on;
% plot(3*PStore_l(2,:),'r');plot(-3*PStore_l(2,:),'r')
% ylabel('y');
% matlab2tikz( 'FY-res.tex' );
% figure(7);
% % subplot(3,1,3);
% plot(XErrStore_l(3,:));hold on;
% plot(3*PStore_l(3,:),'r');plot(-3*PStore_l(3,:),'r')
% ylabel('\theta');
% matlab2tikz( 'FT-res.tex' );

