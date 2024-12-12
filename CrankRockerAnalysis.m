%% Initialization.
clear;clc;close all;
N=100;
theta_2=(20);
l1=30;l2=10;l3=22;l4=26;
alpha3=pi/4;a3=1.5;
%Initial guesses for theta3_4jj
theta3_4=[pi/2,pi/2];
theta_3=0*theta_2;
theta_4=0*theta_2;
%% Using fsolve for estimating theta_3 and theta_4.
for i=1:1:N
    xsol=fsolve(@(x)four_bar_equations(x,l1,l2,l3,l4,theta_2(i)),theta3_4);
    theta_3(i)=xsol(1);
    theta_4(i)=xsol(2);
    theta3_4=[xsol(1),xsol(2)];
    xsol(1)
    xsol(2)
end
coupler_xs=l1+l4*cos(theta_4)+a3*l3*cos(theta_3-alpha3);
coupler_ys=l4*sin(theta_4)+a3*l3*sin(theta_3-alpha3);

%% Plotting
for i=1:1:N
    x_coord=[0,l1,l1+l4*cos(theta_4(i)),l1+l4*cos(theta_4(i))+l3*cos(theta_3(i)),0];
    y_coord=[0,0,l4*sin(theta_4(i)),l4*sin(theta_4(i))+l3*sin(theta_3(i)),0];
    ter_x=[l2*cos(theta_2(i)),l1+l4*cos(theta_4(i)),l1+l4*cos(theta_4(i))+a3*l3*cos(theta_3(i)-alpha3),l2*cos(theta_2(i))];
    ter_y=[l2*sin(theta_2(i)),l4*sin(theta_4(i)),l4*sin(theta_4(i))+a3*l3*sin(theta_3(i)-alpha3),l2*sin(theta_2(i))];
    plot(x_coord,y_coord,'k-','LineWidth',2)
    hold on
    plot(x_coord,y_coord, 'bo','MarkerSize',5)
    %plot(ter_x,ter_y,'k-','LineWidth',2)
    %plot(coupler_xs(1:i),coupler_ys(1:i),'r-','LineWidth',1)
    hold off;
    axis([-l1-l4 2*l1+l4 -l1-l4 l1+l4])
    grid on
    title("Crank-rocker",'FontSize',30,'FontName','Palatino Linotype')
    pause(0.25)
end

%% User defined function for loop closure equations.
function F=four_bar_equations(x,l1,l2,l3,l4,theta2)
    F(1)=l1+l4*cos(x(2))+l3*cos(x(1))-l2*cos(theta2);
    F(2)=l4*sin(x(2))+l3*sin(x(1))-l2*sin(theta2);
end



