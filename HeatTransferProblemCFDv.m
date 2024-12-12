clc
clear
close all
%% Initialization
biot_nos=[0.2,1,5,10];
fo=10; %0.1 to 10
M=11;
dx=1/(M-1);
dt=1e-4;
t_step=9.9/dt;
figure(1)
for bi=biot_nos
    theta=ones(M,t_step+1); 
    theta(M,1)=1/(1+bi*dx);
    for t=2:1:t_step+1
        for x=2:1:M-1
            theta(x,t)=theta(x,t-1)+dt*(theta(x-1,t-1)-2*theta(x,t-1)+theta(x+1,t-1))/(dx^2);
        end
        theta(1,t)=theta(2,t);
        theta(M,t)=theta(M-1,t)/(1+bi*dx);
    end
    plot(1/bi,theta(6,39000)/theta(1,39000),'o','MarkerFaceColor','b');
    hold on;
    xscale("log")
    xlim([0.01,10])
    ylim([0,1])
end
xlabel('1/Bi');
ylabel('theta/theta_0');