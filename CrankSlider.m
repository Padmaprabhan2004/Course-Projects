%% Initialization
clear;clc;close all;
N=50;
theta_2=linspace(45*pi/180,105*pi/180,N);
l2=11.57;l3=21.13;e=-5.48;
%Initial guesses for pred--[S,theta_3]
pred=[6,5*pi/6];
S=0*theta_2;
theta_3=0*theta_2;
b=0.5;

%% Using fsolve for convergence of S and theta_3
for i=1:1:N
    xsol=fsolve(@(x)crank_slider_eqn(x,l2,l3,e,theta_2(i)),pred);
    S(i)=xsol(1);
    theta_3(i)=xsol(2);
    pred=[xsol(1),xsol(2)];
end

%% Plotting
for i=1:1:N
    x_coord=[0,l2*cos(theta_2(i)),S(i)];
    y_coord=[0,l2*sin(theta_2(i)),-e];

    slider_x=[S(i)-b,S(i)+b,S(i)+b,S(i)-b,S(i)-b];
    slider_y=[-e+b,-e+b,-e-b,-e-b,-e+b];
    plot(x_coord,y_coord,'k-','LineWidth',2)
    hold on
    plot(slider_x,slider_y,'k-')
    plot(x_coord,y_coord,'bo')
    hold off;
    axis([-l2-l3 l2+l3 -l2-l3 l2+l3])
    grid on
    title("Crank-slider",'FontSize',30,'FontName','Palatino Linotype')
    pause(0.05)
end

%% UDF for crank slider loop closure equation.
function F=crank_slider_eqn(x,l2,l3,e,theta_2)
    F(1)=x(1)+l3*cos(x(2))-l2*cos(theta_2);
    F(2)=-e+l3*sin(x(2))-l2*sin(theta_2);
end