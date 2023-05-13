function [] = Postproc(COOR_v,COOR_p,u,v,p)
% This function is meant to plot the results obtained.
% Inputs:
%   - COOR_v: Matrix containing the coordinates of the nodes for the v-mesh
%   - COOR_p: Matrix containing the coordinates of the nodes for the p-mesh
%   - u: vector with the horizontal component of the velocity
%   - v: vector with the vertical component of the velocity
%   - p: vector with the pressure
    
    % Nodes coordinates
    x = COOR_v(:,1);
    y = COOR_v(:,2);
    xp = COOR_p(:,1);
    yp = COOR_p(:,2);
    
    
    figure();

    
    % u plot
    interpolant = scatteredInterpolant(x,y,u');
    % Grid
    [xx,yy] = meshgrid(unique(x),unique(y));  
    % Interpolate
    intensity_interp = interpolant(xx,yy);
    % Plot
    subplot(2,2,1)
    [~,h] = contourf(xx,yy,intensity_interp,300);
    set(h,'LineColor','none');
    colormap jet(300)
    shading interp
    axis equal
    colorbar
    title('$u~[m~s^{-1}]$','interpreter','latex')
    xlabel('$x~[m]$','interpreter','latex')
    ylabel('$y~[m]$','interpreter','latex')


    % v plot
    interpolant = scatteredInterpolant(x,y,v');
    % Grid
    [xx,yy] = meshgrid(unique(x),unique(y));
    % Interpolate
    intensity_interp = interpolant(xx,yy);
    % Plot
    subplot(2,2,2)
    [~,h] = contourf(xx,yy,intensity_interp,300);
    set(h,'LineColor','none');
    colormap jet(300)
    shading interp
    axis equal
    colorbar
    title('$v~[m~s^{-1}]$','interpreter','latex')
    xlabel('$x~[m]$','interpreter','latex')
    ylabel('$y~[m]$','interpreter','latex')

    % Pressure Plot
    interpolant = scatteredInterpolant(xp,yp,p);
    % Grid
    [xx,yy] = meshgrid(unique(xp),unique(yp)); 
    % Interpolate
    intensity_interp = interpolant(xx,yy);
    % Plot
    subplot(2,2,3)
    [~,h] = contourf(xx,yy,intensity_interp,200);
    set(h,'LineColor','none');
    colormap jet(200)
    axis equal
    shading interp
    colorbar
    title('$p~[Pa]$','interpreter','latex')
    xlabel('$x~[m]$','interpreter','latex')
    ylabel('$y~[m]$','interpreter','latex')

    % Streamlines Plot
    subplot(2,2,4)
    interpolantu = scatteredInterpolant(x,y,u');
    interpolantv = scatteredInterpolant(x,y,v');
    intensity_interpu = interpolantu(xx,yy);
    intensity_interpv = interpolantv(xx,yy);
    k = 4;
    quiver(x,y,u',v',2);
    starty = [linspace(0,1,k),linspace(0,0.1,0),linspace(0,0.1,0)];
    startx = [0.5*ones(k,1);linspace(0,0.1,0)';linspace(0.9,1,0)'];
% startx = 1*rand(k,1); 
% starty = 1*rand(k,1); [0.5,100000]
    streamline(xx,yy,intensity_interpu,intensity_interpv,startx,starty,[0.5,100000])
    axis equal
    axis([0 1 0 1])
    title('Streamlines and Velocity Field','interpreter','latex')
    xlabel('$x~[m]$','interpreter','latex')
    ylabel('$y~[m]$','interpreter','latex')
end