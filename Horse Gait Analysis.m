clc;clear;close;
l1 = 5; l2 = 1.7985; l3 = 5.002; l4 = 2.211; l5 = 4.3707; l6 = 5.1081;
a1=0.93818; alpha1=pi/4; a3 =0.7149; alpha3=0.261;alpha6=pi/6;
theta1=8*pi/180;a6=1.5388;
N=200;
theta_3_4_5_and_6 = [pi/6,pi/3,pi/4,pi/2];
theta2 = linspace(0,4*pi,N);
theta3 = 0*theta2;
theta4 = 0*theta2;
theta5 = 0*theta2;
theta6 = 0*theta2;
%% Solving nonlinear equations, storing the solutions, and updating initial guess
for i=1:N
    xsol = fsolve(@(x)loop_closure_6bar_stephenson3(x,theta1,l1,l2,l3,l4,l5,l6,a1,alpha1,a3,alpha3,theta2(i)),theta_3_4_5_and_6);
    theta3(i) = xsol(1);
    theta4(i) = xsol(2);
    theta5(i) = xsol(3);
    theta6(i) = xsol(4);
    theta_3_4_5_and_6 = [ theta3(i), theta4(i), theta5(i), theta6(i)];
end

ext_xs=l4*cos(theta4)+a3*l3*cos(theta3-alpha3)+a6*l6*cos(theta6+alpha6);
ext_ys=l4*sin(theta4)+a3*l3*sin(theta3-alpha3)+a6*l6*sin(theta6+alpha6);
%% Plotting the movement of the mechanism    
for i=N:-1:1
    %UDF order-
    %x_coord_loop_1 = [0,l1*cos(theta1),l1*cos(theta1)+l2*cos(theta2(i)),l4*cos(theta4(i)),0];
    %y_coord_loop_1 = [0,l1*sin(theta1),l1*sin(theta1)+l2*sin(theta2(i)),l4*sin(theta4(i)),0];
    %x_coord_loop_2=[a1*l1*cos(theta1+alpha1),a1*l1*cos(theta1+alpha1)+l5*cos(theta5(i)),a1*l1*cos(theta1+alpha1)+l5*cos(theta5(i))+l6*cos(theta6(i)),l4*cos(theta4(i)),0];
    %y_coord_loop_2=[a1*l1*sin(theta1+alpha1),a1*l1*sin(theta1+alpha1)+l5*sin(theta5(i)),a1*l1*sin(theta1+alpha1)+l5*sin(theta5(i))+l6*sin(theta6(i)),l4*sin(theta4(i)),0];
    %extension_x=[l4*cos(theta4(i))+a3*l3*cos(theta3(i)-alpha3),l4*cos(theta4(i))+a3*l3*cos(theta3(i)-alpha3)+a6*l6*cos(theta6(i)+alpha6)];
    %extension_y=[l4*sin(theta4(i))+a3*l3*sin(theta3(i)-alpha3),l4*sin(theta4(i))+a3*l3*sin(theta3(i)-alpha3)+a6*l6*sin(theta6(i)+alpha6)];
    leg1=coords(i,l1,l2,l3,l4,l5,l6,a1,alpha1,a6,alpha6,a3,alpha3,theta1,theta2,theta3,theta4,theta5,theta6,N);
    leg2=coords(i+floor(N/4),l1,l2,l3,l4,l5,l6,a1,alpha1,a6,alpha6,a3,alpha3,theta1,theta2,theta3,theta4,theta5,theta6,N);
    leg3=coords(i+floor(N/4),l1,l2,l3,l4,l5,l6,a1,alpha1,a6,alpha6,a3,alpha3,theta1,theta2,theta3,theta4,theta5,theta6,N);
    leg4=coords(i,l1,l2,l3,l4,l5,l6,a1,alpha1,a6,alpha6,a3,alpha3,theta1,theta2,theta3,theta4,theta5,theta6,N);
    mod(i+floor(N/4),N+1)
    figure(1)
    plot(leg1(1,:),leg1(2,:), 'b-',LineWidth=2)
    hold on;
    plot(leg1(3,:),leg1(4,:), 'b-',LineWidth=2)
    plot(leg1(5,1:2),leg1(6,1:2),'b-',LineWidth=2)
    %Leg2
    plot(leg2(1,:)+10,leg2(2,:), 'b-',LineWidth=2)
    plot(leg2(3,:)+10,leg2(4,:), 'b-',LineWidth=2)
    plot(leg2(5,1:2)+10,leg2(6,1:2),'b-',LineWidth=2)
    %Leg3
    plot(leg3(1,:),leg2(2,:), 'b:',LineWidth=1.5)
    plot(leg3(3,:),leg2(4,:), 'b:',LineWidth=1.5)
    plot(leg3(5,1:2),leg2(6,1:2),'b:',LineWidth=1.5)
    %Leg4
    plot(leg4(1,:)+10,leg2(2,:), 'b:',LineWidth=1.5)
    plot(leg4(3,:)+10,leg2(4,:), 'b:',LineWidth=1.5)
    plot(leg4(5,1:2)+10,leg2(6,1:2),'b:',LineWidth=1.5)
    %Path traced
    plot(ext_xs(i:N),ext_ys(i:N),'r.',LineWidth=1)
    hold off;
    axis([-3*l1-l6 3*l1+l6 -3*l1-l6 3*l1+l6])
    axis("equal")
    grid on
    title("Horse Gait Mechanism",'FontSize',20,'FontName','Palatino Linotype')
    pause(0.05)
end
%% Loop-closure equations defined as a function for fsolve() to call
function F=loop_closure_6bar_stephenson3(x,theta1,l1,l2,l3,l4,l5,l6,a1,alpha1,a3,alpha3,theta2)
F(1) = l1*cos(theta1)+l2*cos(theta2)+l3*cos(x(1))-l4*cos(x(2));
F(2) = l1*sin(theta1)+l2*sin(theta2)+l3*sin(x(1))-l4*sin(x(2));
F(3) = a1*l1*cos(theta1+alpha1)+l5*cos(x(3))+l6*cos(x(4))-a3*l3*cos(x(1)-alpha3)-l4*cos(x(2));
F(4) = a1*l1*sin(theta1+alpha1)+l5*sin(x(3))+l6*sin(x(4))-a3*l3*sin(x(1)-alpha3)-l4*sin(x(2));
end
%% Coordinates generated an matrix for each leg with a function
function A = coords(i, l1, l2, l3, l4, l5, l6, a1, alpha1, a6, alpha6, a3, alpha3, theta1, theta2, theta3, theta4, theta5, theta6,N)
    if mod(i,N+1)==0
        i=1;
    end
    A(1, 1) = 0;
    A(1, 2) = l1 * cos(theta1);
    A(1, 3) = l1 * cos(theta1) + l2 * cos(theta2(mod(i,N+1)));
    A(1, 4) = l4 * cos(theta4(mod(i,N+1)));
    A(1, 5) = 0;

    A(2, 1) = 0;
    A(2, 2) = l1 * sin(theta1);
    A(2, 3) = l1 * sin(theta1) + l2 * sin(theta2(mod(i,N+1)));
    A(2, 4) = l4 * sin(theta4(mod(i,N+1)));
    A(2, 5) = 0;

    A(3, 1) = a1 * l1 * cos(theta1 + alpha1);
    A(3, 2) = a1 * l1 * cos(theta1 + alpha1) + l5 * cos(theta5(mod(i,N+1)));
    A(3, 3) = a1 * l1 * cos(theta1 + alpha1) + l5 * cos(theta5(mod(i,N+1))) + l6 * cos(theta6(mod(i,N+1)) );
    A(3, 4) = l4 * cos(theta4(mod(i,N+1)));

    A(4, 1) = a1 * l1 * sin(theta1 + alpha1);
    A(4, 2) = a1 * l1 * sin(theta1 + alpha1) + l5 * sin(theta5(mod(i,N+1)) );
    A(4, 3) = a1 * l1 * sin(theta1 + alpha1) + l5 * sin(theta5(mod(i,N+1)) ) + l6 * sin(theta6(mod(i,N+1)) );
    A(4, 4) = l4 * sin(theta4(mod(i,N+1)));

    A(5, 1) = l4 * cos(theta4(mod(i,N+1))) + a3 * l3 * cos(theta3(mod(i,N+1)) - alpha3 );
    A(5, 2) = l4 * cos(theta4(mod(i,N+1))) + a3 * l3 * cos(theta3(mod(i,N+1)) - alpha3 ) + a6 * l6 * cos(theta6(mod(i,N+1)) + alpha6);

    A(6, 1) = l4 * sin(theta4(mod(i,N+1))) + a3 * l3 * sin(theta3(mod(i,N+1)) - alpha3 );
    A(6, 2) = l4 * sin(theta4(mod(i,N+1))) + a3 * l3 * sin(theta3(mod(i,N+1)) - alpha3 ) + a6 * l6 * sin(theta6(mod(i,N+1)) + alpha6);
end