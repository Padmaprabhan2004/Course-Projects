clc; close all;

%% given
L=1;
H=1;
U=100% velocity
Re=100;
nx=50;
ny=50;
delt=1e-6;
tol=1e-6;
delx=1/(nx-1);
dely=1/(ny-1);
nu=U*L/Re;

%% init
u=zeros(ny,nx);
v=zeros(ny,nx);
psi=zeros(ny,nx);
omega=zeros(ny,nx);
error=10;

%% ic
u=apply_x_velocity_BC(U,nx,ny,u);
v=apply_y_velocity_BC(U,nx,ny,v);
psi=apply_stream_function_BC(nx,ny,psi);
omega=apply_vorticity_BC(U,nx,ny,omega,psi,delx,dely);

timestep=0;
while error>tol
    timestep=timestep+1;
    omega_old=omega;
    for i=2:1:ny-1
        for j=2:nx-1
            %upwind scheme
            a = nu * ((omega(i+1, j) - 2 * omega(i, j) + omega(i-1, j)) / (dely^2) + ...
                      (omega(i, j+1) - 2 * omega(i, j) + omega(i, j-1)) / (delx^2));

            dwdx=0;
            dwdy=0;
            if(u(i,j)>0)
                dwdx=(omega(i,j)-omega(i,j-1))/delx;
            end
            if(u(i,j)<0)
                dwdx=(omega(i,j+1)-omega(i,j))/delx;
            end
            if(v(i,j)>0)
                dwdy=(omega(i,j)-omega(i-1,j))/dely;
            end
            if(v(i,j)<0)
                dwdy=(omega(i+1,j)-omega(i,j))/dely;
            end
            b=u(i,j)*dwdx;
            c=v(i,j)*dwdy;
            omega(i,j)=omega(i,j)+delt*(a-b-c);
        end
    end
    psi_error=10;
    while psi_error>tol
        psi_old=psi;
        for i=2:ny-1
            for j=2:nx-1
                psi(i,j)=0.25*(omega(i,j)*(dely^2)+psi(i,j-1)+psi(i,j+1)+psi(i-1,j)+psi(i+1,j));
            end
        end
        psi_err= max(max(abs(psi - psi_old)));
    end
    for i = 2:ny-1
        for j = 2:nx-1
            u(i, j) = (psi(i+1, j) - psi(i-1, j)) / (2 * dely);
            v(i, j) = (psi(i, j-1) - psi(i, j+1)) / (2 * delx);
        end
    end
    omega = apply_vorticity_BC(U, nx, ny, omega, psi, delx, dely);
    error = max(max(abs(omega - omega_old)));
end







%% Helper functions
function u = apply_x_velocity_BC(U, nx, ny, u)
    % Apply boundary conditions for horizontaol velocity
    u(1, :) = 0; % Bottom wall
    u(ny, :) = U; % Top wall
    u(:, 1) = 0; % Left wall
    u(:, nx) = 0; % Right wall
end

function v = apply_y_velocity_BC(U,nx, ny, v)
    % Apply boundary conditions for vertical velocity
    v(1, :) = 0; % Bottom wall
    v(ny, :) = 0; % Top wall
    v(:, 1) = 0; % Left wall
    v(:, nx) = 0; % Right wall
end

function psi = apply_stream_function_BC(nx, ny, psi)
    % Apply boundary conditions for stream function
    psi(1, :) = 0; % Bottom wall
    psi(ny, :) = 0; % Top wall
    psi(:, 1) = 0; % Left wall
    psi(:, nx) = 0; % Right wall
end

function omega = apply_vorticity_BC(U, nx, ny, omega, psi, delx, dely)
    omega(1, 1:nx)=-2*(psi(2, 1:nx)-psi(1,1:nx))/ dely^2; % Bottom wall
    omega(end, 1:nx)=-2*(psi(ny-1, 1:nx)-psi(ny,1:nx)) / dely^2 - 2 * U / dely; % Top wall
    omega(1:ny, 1)=-2*(psi(1:ny,2)-psi(1:ny,1))/ delx^2; % Left wall
    omega(1:ny, end)=2*(psi(1:ny,nx)-psi(1:ny, nx-1)) / delx^2; % Right wall
end

%Plotting--------------PLEASE CHANGE THIS ONLY AND TUNE THE ABOVE CODE
%ACCORDINGLY
x = 0:delx:L;
y = 0:dely:L;
V = sqrt((u.^2) + (v.^2));
figure;
pcolor(x,y,u);
daspect([1 1 1])
shading interp
colormap jet
map = colorbar;
map.Label.String = 'X-Velocity (m/s)';

figure;
pcolor(x,y,v);
daspect([1 1 1])
shading interp
colormap jet
map = colorbar;
map.Label.String = 'Y-Velocity (m/s)';

figure;
pcolor(x,y,omega_new);
daspect([1 1 1])
shading interp
colormap jet
map = colorbar;
map.Label.String = 'Vorticity (1/s)';

figure;
pcolor(x,y,psi);
daspect([1 1 1])
shading interp
colormap jet
map = colorbar;
map.Label.String = 'Stream Function (\psi)';

[X,Y] = meshgrid(x,y);
figure;
contourf(X,Y,psi,21);
daspect([1 1 1])
colormap jet
map = colorbar;
map.Label.String = 'Stream Function (\psi)';

%Validation
Data = readmatrix(".\data.xlsx");
y_ghia = Data(1:17,1);
u_ghia_Re400 = Data(1:17,2);
u_ghia_Re100 = Data(1:17,8);
x_ghia = Data(1:17,4);
v_ghia_Re400 = Data(1:17,5);
v_ghia_Re100 = Data(1:17,11);

if Re == 100
    figure;
    hold on
    grid on
    plot(u(:,0.5*(nx+1)),y,'k-','DisplayName', 'Present Solution','LineWidth',1.5,'MarkerSize',10);
    plot(u_ghia_Re100,y_ghia,'ro','DisplayName', 'Ghia Solution','LineWidth',1.5,'MarkerSize',10);
    title('Variation of x-Velocity along Vertical Line through Geometric Center for Re = 100')
    xlabel('u (m/s)');
    legend('Location', 'best');
    ylabel('y (m)');
    hold off

    figure;
    hold on
    grid on
    plot(x,v(0.5*(nx+1),:),'k-','DisplayName', 'Present Solution','LineWidth',1.5,'MarkerSize',10);
    plot(x_ghia,v_ghia_Re100,'ro','DisplayName', 'Ghia Solution','LineWidth',1.5,'MarkerSize',10);
    title('Variation of y-Velocity along Horizontal Line through Geometric Center for Re = 100')
    xlabel('x (m)');
    legend('Location', 'best');
    ylabel('v (m/s)');
    hold off
end

if Re == 400
    figure;
    hold on
    grid on
    plot(u(:,0.5*(nx+1)),y,'k-','DisplayName', 'Present Solution','LineWidth',1.5,'MarkerSize',10);
    plot(u_ghia_Re400,y_ghia,'ro','DisplayName', 'Ghia Solution','LineWidth',1.5,'MarkerSize',10);
    title('Variation of x-Velocity along Vertical Line through Geometric Center for Re = 400')
    xlabel('u (m/s)');
    legend('Location', 'best');
    ylabel('y (m)');
    hold off

    figure;
    hold on
    grid on
    plot(x,v(0.5*(nx+1),:),'k-','DisplayName', 'Present Solution','LineWidth',1.5,'MarkerSize',10);
    plot(x_ghia,v_ghia_Re400,'ro','DisplayName', 'Ghia Solution','LineWidth',1.5,'MarkerSize',10);
    title('Variation of y-Velocity along Horizontal Line through Geometric Center for Re = 400')
    xlabel('x (m)');
    legend('Location', 'best');
    ylabel('v (m/s)');
    hold off
end

%function definitions
function error = maximum(nx,ny,diff)
    error = diff(1,1);
    for i = 1:ny
        for j = 1:nx
            if abs(diff(i,j)) > error
               error = abs(diff(i,j));
            end
        end
    end
end

function u = x_velocity_BC(U,nx,ny,u)
    for j = 1:nx %Bottom wall
        u(1,j) = 0;
    end
    for i = 1:ny %Left wall
        u(i,1) = 0;
    end
    for i = 1:ny %Right wall
        u(i,nx) = 0;
    end
    for j = 1:nx %Top wall
        u(ny,j) = U;
    end
end

function v = y_velocity_BC(nx,ny,v)
    for j = 1:nx %Bottom wall
        v(1,j) = 0;
    end
    for i = 1:ny %Left wall
        v(i,1) = 0;
    end
    for i = 1:ny %Right wall
        v(i,nx) = 0;
    end
    for j = 1:nx %Top wall
        v(ny,j) = 0;
    end
end

function psi = stream_function_BC(nx,ny,psi)
    for j = 1:nx %Bottom wall
        psi(1,j) = 0;
    end
    for i = 1:ny %Left wall
        psi(i,1) = 0;
    end
    for i = 1:ny %Right wall
        psi(i,nx) = 0;
    end
    for j = 1:nx %Top wall
        psi(ny,j) = 0;
    end
end

function omega = vorticity_BC(U,nx,ny,omega,psi,delx,dely)
    for j = 1:nx %Bottom wall
        omega(1,j) = 2*(psi(1,j)-psi(2,j))/(dely^2);
    end
    for i = 1:ny %Left wall
        omega(i,1) = 2*(psi(i,1)-psi(i,2))/(delx^2);
    end
    for i = 1:ny %Right wall
        omega(i,nx) = 2*(psi(i,nx)-psi(i,nx-1))/(delx^2);
    end
    for j = 1:nx %Top wall
        omega(ny,j) = 2*(psi(ny,j)-psi(ny-1,j)-(dely*U))/(dely^2);
    end
end

function [omega_new,omega_error] = omega_solver(nx,ny,u,v,omega,omega_new,delx,dely,delt,nu)
    omega_diff = zeros(ny,nx);
    for i = 2:ny-1
        for j = 2:nx-1
            omega_old = omega_new(i,j);
            a = delt*nu/(delx^2);
            b = delt*nu/(dely^2);
            c = delt*u(i,j)/(2*delx);
            d = delt*v(i,j)/(2*dely);
            omega_new(i,j) = (((a+c)*omega_new(i,j-1)) + ((b+d)*omega_new(i-1,j)) + ((a-c)*omega_new(i,j+1)) + ((b-d)*omega_new(i+1,j)) + omega(i,j))/(1 + (2*(a + b)));
            omega_diff(i,j) = omega_new(i,j) - omega_old;
        end
    end
    omega_error = maximum(nx,ny,omega_diff);
end

function [psi,psi_error] = psi_solver(nx,ny,omega_new,psi,delx,dely)
    psi_diff = zeros(ny,nx);
    for i = 2:ny-1
        for j = 2:nx-1
            psi_old = psi(i,j);
            a = (psi(i,j-1)+psi(i,j+1))/(delx^2);
            b = (psi(i-1,j)+psi(i+1,j))/(dely^2);
            c = 2/(delx^2);
            d = 2/(dely^2);
            psi(i,j) = (omega_new(i,j)+a+b)/(c+d);
            psi_diff(i,j) = psi(i,j) - psi_old;
        end
    end
    psi_error = maximum(nx,ny,psi_diff);
end

function [u,v] = velocity_field(nx,ny,u,v,psi,delx,dely)
    for i = 2:ny-1
        for j = 2:nx-1
            u(i,j) = (psi(i+1,j) - psi(i-1,j))/(2*dely);
            v(i,j) = (psi(i,j-1) - psi(i,j+1))/(2*delx);
        end
    end
end