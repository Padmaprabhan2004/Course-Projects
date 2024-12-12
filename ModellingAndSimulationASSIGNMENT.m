clc
clear all
close all
%% Initialization
t=linspace(0,10,1000000);
dt=t(2)-t(1);
s1=zeros(4,length(t)); %s1-state of mass 1 
s2=zeros(4,length(t)); %s2- state of mass 2

%% Known Values
m1=1;
m2=1;
k=1000;
c=5;
L=0.5;
g=[0,0]';

%% Initial Values
s1(1,1)=0;
s1(2,1)=0;
s1(3,1)=0;
s1(4,1)=1;

s2(1,1)=0.5;
s2(2,1)=0;
s2(3,1)=0;
s2(4,1)=-1;

%% Euler explicit integration
for i=1:1:length(t)-1
    s1(1:2,i+1)=s1(1:2,i)+s1(3:4,i)*dt;
    s2(1:2,i+1)=s2(1:2,i)+s2(3:4,i)*dt;
    
    w=s1(1:2,i)-s2(1:2,i);
    wd=s1(3:4,i)-s2(3:4,i);
    
    %spring forces on m1 and m2
    Fs1=-k*(norm(w,2)-L)*(w/norm(w,2));
    Fs2=-k*(norm(w,2)-L)*(-w/norm(w,2));
    % damping forces on m1 and m2
    Fd1=-c*wd;
    Fd2=c*wd;
    
    s1(3:4,i+1)=s1(3:4,i)+dt*(Fs1+Fd1+m1*g)/m1;
    s2(3:4,i+1)=s2(3:4,i)+dt*(Fs2+Fd2+m2*g)/m2;
end
%% Visualizing the system

figure(4)
for idx=1:10000:length(t)
    plot([s1(1,idx),s2(1,idx)],[s1(2,idx),s2(2,idx)],'k','LineWidth',2,'Marker','o');
    axis equal;
    title('Visualization')
    pause(0.01);
    hold off;
end
%% Plotting
figure(1)
subplot(2,1,1)
plot(t,s1(1,:),"Color",'k','LineWidth',2)
title('m1(x)')
xlabel('Time')
ylabel('X')
subplot(2,1,2)
plot(t,s1(2,:),"Color",'b','LineWidth',2)
title('m1(y)')
xlabel('Time')
ylabel('Y')

figure(2)
subplot(2,1,1)
plot(t,s2(1,:),"Color",'k','LineWidth',2)
title('m2(x)')
xlabel('Time')
ylabel('X')
subplot(2,1,2)
plot(t,s2(2,:),"Color",'b','LineWidth',2)
title('m2(y)')
xlabel('Time')
ylabel('Y')



%% Energy calculation
total_energies=zeros(1,length(t));
for i=1:1:length(t)
    ke=m1*(s1(3,i)^2+s1(4,i)^2)*0.5+m2*(s2(3,i)^2+s2(4,i)^2)*0.5;
    w=s1(1:2,i)-s2(1:2,i);
    pe=m1*9.81*s1(2,i)+m2*9.81*s2(2,i)+0.5*k*(norm(w,2)-L)^2;
    total_energies(1,i)=ke+pe;
end
figure(3)
plot(t,total_energies,'LineWidth',2)
title('Total Energy of the system')
xlabel('Time')
ylabel('Total Energy')




