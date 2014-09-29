function DrawCar(Xr,omega, color)

    L = 0.255; % Length of vehicle
    W = 0.18; % Width of vehicle
    D = 0.04; % Wheel diameter 
    
    a = W/2;
    b = L/2;
    d = D/2;
    w = 0.12/2;
    h = 0.35/2;
%     omega = -45*pi/180;
    m = [cos(omega) -sin(omega); sin(omega) cos(omega)]*[0 d]';
    
%     P = [-a -b ; a -b ; a b; -a b]'; % Coordinates vector for a simple quadritaleral
    P = [a b ; -a b; ... %front axis
         0 b ; 0 -b; ... %rear axis
         a -b ; -a -b; ... % center axis
         a b; (a-m(1)) (b-m(2)); ...% FR
         -a b; (-a-m(1)) (b-m(2)); ... % FL -a (b+d); -a (b-d) 
         a (-b+d); a (-b-d); ... % RR
         -a (-b+d); -a (-b-d); ... % RL
         -w -h ; w -h ; w h; -w h; ...
         a b; (a+m(1)) (b+m(2));
         -a b; (-a+m(1)) (b+m(2))]'; 
     
    theta = Xr(3) - pi/2; % Theta angle for rotation

    P=[cos(theta) -sin(theta); ...
       sin(theta) cos(theta)]*P; % Rigid rotation
    
%     P(:,7:10) = [cos(omega) -sin(omega); ...
%                  sin(omega) cos(omega)]*P(:,7:10);

    P(1,:)=P(1,:) + Xr(1); % Shifting in X axis
    P(2,:)=P(2,:) + Xr(2); % Shifting in Y axis

%     A = zeros(4,4); % Connection nodes matrix
%     A(1,2) = 1; % n1 to n2; Node 1 connected to node 2
%     A(2,3) = 1; % n2 to n3
%     A(3,4) = 1; % n3 to n4
%     A(4,1) = 1; % n4 to n1
    A = zeros(22,22); % Connection nodes matrix
    A(1,2) = 1; % n1 to n2; Node 1 connected to node 2
    A(3,4) = 1; % n3 to n4
    A(5,6) = 1;
    A(7,8) = 1;
    A(9,10) = 1;
    A(11,12) = 1;
    A(13,14) = 1;
    A(15,16) = 1; % n1 to n2; Node 1 connected to node 2
    A(16,17) = 1; % n2 to n3
    A(17,18) = 1; % n3 to n4
    A(18,15) = 1; % n4 to n1
    A(19,20) = 1;
    A(21,22) = 1;
    [x y] = gplot(A,P');

    plot(x,y,color);

end
