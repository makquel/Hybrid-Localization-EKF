clc
close all
clear all


%% GPS Simulator Plot

Z_gps = importdata('TRES/gps.csv');
x = Z_gps(:,1);
y = Z_gps(:,2);
theta = Z_gps(:,3);
SNR = 28;
X_gps = awgn(Z_gps(:,1),SNR,'measured');
Y_gps = awgn(Z_gps(:,2),SNR,'measured');
T_gps = awgn(Z_gps(:,3),SNR,'measured');
scale = 0.3;
n = 1; % minimum observationindex
N = length(Z_gps);
dT = 0.252; % Sampling time
%%

% figure(1);hold on;axis([-0 2.6 -0 1.75]);grid on;
% for k = n:N
%     quiver(Z_gps(k,1),Z_gps(k,2),cosd(Z_gps(k,3)),sind(Z_gps(k,3)),scale,'k');
% %     plot(Z_gps(k,1),Z_gps(k,2),'ko');
%     if(k > n)
%         line([Z_gps(k-1,1),Z_gps(k,1)],[Z_gps(k-1,2),Z_gps(k,2)],'color','r');%-.b
%      
%     end;
% %     pause(dT); % pause a bit to see animation   
% end
% % matlab2tikz( 'pre-GPS.tikz' );
% hleg = legend('Orientacao','Posicao');
% set(hleg,'Location','NorthWest');
% set(hleg,'Interpreter','none');


%% Animation
fig_pos = [0 0 533 400]; % position and size of the figure window
fillplot_ax_pos = [30 30 500 370];        % position and size of fill plot

fig_col = [1 1 1];  % figure background color
text_col = [0 0 0]; % text color
light_grey = [.4 .4 .4];
dark_grey = [.2 .2 .2];
nice_blue = [51/255 51/255 102/255];
light_red = [.6 .4 .4];

movieflag = 1;
moviefilename = 'gps.avi';
 
% only if our flag is set to 1, will we open up a movie object:
if movieflag == 1
    aviobj = avifile(moviefilename, 'fps', 2, 'compression', 'none');
end

fh= figure('color', fig_col, 'name', 'Tutorial animation movie','Position', fig_pos);
grid on;

tic;
fprintf('\nWe are entering the loop...\n');
% ***************  START THE BIG LOOP AFTER CREATING FIG. ***************
for k = n:N
 
%     titlestr = sprintf('GPS(%03d)', k); % our title will change
%     fprintf('Frame %d\n', k); % print loop counter
 
 
    % ***************  1) GPS  ***************
    fillplot_ax = axes;
    set(fillplot_ax, 'Units', 'pixels', 'Position', fillplot_ax_pos);
%     set(fillplot_ax, 'Units', 'pixels', 'Position', [50 50 1024-60 100]); 
    p1h(1) = plot(x,y, '-.b');
%     set(p1h(1), 'Color', light_red);
    hold on;   
    p1h(2)=plot(X_gps(1:k), Y_gps(1:k),'-r');
    p1h(3) = quiver(X_gps(1:k),Y_gps(1:k),cosd(T_gps(1:k)),sind(T_gps(1:k)),0.3,'k');
    hold off;
    axis([-0 2.6 -0 1.75]);
%     axis square;
%     th=title(titlestr);
    xlabel('x'); ylabel('y');
%     set(th, 'FontSize', 14); % make title font larger
 
    
     
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


% hold on
% h1=plot(t,x1);
% h2=plot(t,x2);
% h3=plot(t,x3);
% h4=plot(t,x4);
% hold off
% legend([h2 h4],'Second','Fourth')
%% 
% figure(3);
% 
quiver(X_gps,Y_gps,cosd(T_gps),sind(T_gps),0.5,'k');
hold on;
grid on;
plot(X_gps,Y_gps,'r');hold on;
plot(x,y,'-.b');hold on;
axis([0 2.3 0 1.75]);
hleg = legend('Orientacao','Posicao','Ground truth');
set(hleg,'Location','NorthWest');
set(hleg,'Interpreter','none');
% matlab2tikz( 'GPS.tikz' );


M = 5;
k = 9;
% k = 13;
X = ones(1,M);
Y = ones(1,M);
Theta = ones(1,M);
for i = 1:M;
    X(i) = X(i)*awgn(Z_gps(k,1),SNR,'measured'); 
    Y(i) = Y(i)*awgn(Z_gps(k,2),SNR,'measured');
    Theta(i) = Theta(i)*awgn(Z_gps(k,3),SNR,'measured');
end

dx = Z_gps(k,1);
dy = Z_gps(k,2);
dtheta = Z_gps(k,3);
Rx = 0;
Ry = 0;
Rtheta = 0;
for k = 1:M
    Rx = Rx + (X(k) - dx)*(X(k) - dx)';
    Ry = Ry + (Y(k) - dy)*(Y(k) - dy)';
    Rtheta = Rtheta + (Theta(k) - dtheta)*(Theta(k) - dtheta)';
end 

Rx = Rx/M;
Ry = Ry/M;
Rtheta = Rtheta/M;

R_covariance = diag([Rx Ry Rtheta]);

%% 
% figure(4)
% I = imread('0_igps_rect.jpg');
% imshow(I);
% hold on;
% 
% %# define points (in matrix coordinates)
% for obs = 1:length(Z_gps)-1
%     p1 = [Z_gps(obs,2)/0.003571429,Z_gps(obs,1)/0.003571429];
%     p2 = [Z_gps(obs+1,2)/0.003571429,Z_gps(obs+1,1)/0.003571429];
%     plot([p1(2),p2(2)],[p1(1),p2(1)],'Color','r','LineWidth',2,'LineStyle','-.'); %-.b
% end
% %imwrite(I,'gps.png','BitDepth',16);
% % 0.003571429

%% Test graph
% figure
% stream = RandStream('mrg32k3a','Seed',4);
% y1 = rand(stream,10,5);
% hb = bar(y1,'stacked'); colormap(summer); hold on
% y2 = rand(stream,10,1);
% hp = plot(1:10,y2,'marker','square','markersize',12,...
%        'markeredgecolor','y','markerfacecolor',[.6,0,.6],...
%        'linestyle','-','color','r','linewidth',2); hold off
% legend([hb,hp],'Carrots','Peas','Peppers','Green Beans',...
%         'Cucumbers','Eggplant','Location','SouthEastOutside')
