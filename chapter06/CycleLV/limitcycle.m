% plot limit cycle of sfun

close all;
clear all; 

    sfunction=@sfun;
    alpha = 1; 
    % save the sizes of x1/x2 for use later    
    x1 = [0.0:.2:3.5];
    x2 = [0.0:.2:3.5];
    lenx1 = length(x1);
    lenx2 = length(x2);

    % data variables
    uv = [0 1];
    outx = zeros(lenx1, lenx2); % holds the x positions
    outy = zeros(lenx1, lenx2); % holds the y positions
    outu = zeros(lenx1, lenx2); % holds the x component of the slope
    outv = zeros(lenx1, lenx2); % holds the y component of the slope

    fig = figure(1);

    % points
    x  = 0;
    y  = 0;
    len = .25*(x1(2)-x1(1)); % length of marks
    dx = 0;
    dy = 0;
    a  = 0; % angle
    
    % trajectory variables
    pathcolor = [1 0 0];
    t  = 0;
    dt = .001; % time step to simulate
    simlen = 8; % length of the simulation in seconds

    % create the axis and setup the callbacks
    axis([x1(1)-len, x1(lenx1)+len, x2(1)-len, x2(lenx2)+len]);

    
    % calculates all the slope field vectors for later use
    for i = 1:lenx1
        for j = 1:lenx2
            % store the position of the slope vector
            outx(i,j) = x1(i);
            outy(i,j) = x2(j);

            % store the slope vector
            uv = sfunction(0, [x1(i) x2(j)], alpha);
            outu(i,j) = uv(1);
            outv(i,j) = uv(2);
        end
    end
    
    mode = 1;
    % Loops to plot each tick mark
    if (mode == 1)
        for i=1:lenx1
            for j=1:lenx2
                x = outx(i,j);
                y = outy(i,j);
                a = atan2(outv(i,j),outu(i,j));

                dx = len*cos(a);
                dy = len*sin(a);

                line([x-dx,x+dx], [y-dy,y+dy], 'HitTest','Off',...
                    'Color',[.5 .5 .5])
            end
        end
    else
        quiv=quiver(outx,outy,outu,outv);
        set(quiv,'HitTest','Off');
    end
 
    
    
        % initial point P1
        x = 0.2;
        y = 1.0;
        t = 0;     
        H = alpha*x+y-log(x^alpha * y)
        while ( t < simlen )
        uv  = sfunction(t,[x,y],alpha);
        dx = uv(1)*dt;
        dy = uv(2)*dt;

   %     line([x,x+dx], [y,y+dy], 'HitTest','Off','Color',pathcolor);
        line([x,x+dx], [y,y+dy],'Color','blue','LineStyle','-','Linewidth',1)
     
        x = x + dx;
        y = y + dy;
        t = t + dt;      
        end
        
        % arrow 
        ha = annotation('arrow');  % store the arrow information in ha
        ha.Parent = gca;           % associate the arrow the the current axes
        ha.X = [x x+dx];          % the location in data units
        ha.Y = [y y+dy];   

        ha.Color = 'blue';
        ha.HeadStyle ='vback3';
        ha.LineWidth  = 1;          % make the arrow bolder for the picture
        ha.HeadWidth  = 10;
        ha.HeadLength = 10;
        
        text(x+1000*dx,y+200*dy,'H=2.8')
 
         % initial point P2
        x = 0.5;
        y = 1.0;
        t = 0;   
        H = alpha*x+y-log(x^alpha * y)
        while ( t < simlen )
        uv  = sfunction(t,[x,y],alpha);
        dx = uv(1)*dt;
        dy = uv(2)*dt;

   %     line([x,x+dx], [y,y+dy], 'HitTest','Off','Color',pathcolor);
        line([x,x+dx], [y,y+dy],'Color','blue','LineStyle','-','Linewidth',1)
   
        x = x + dx;
        y = y + dy;
        t = t + dt;      
        end
        
                % arrow 
        ha = annotation('arrow');  % store the arrow information in ha
        ha.Parent = gca;           % associate the arrow the the current axes
        ha.X = [x x+dx];          % the location in data units
        ha.Y = [y y+dy];   

        ha.Color = 'blue';
        ha.HeadStyle ='vback3';
        ha.LineWidth  = 1;          % make the arrow bolder for the picture
        ha.HeadWidth  = 10;
        ha.HeadLength = 10;
        
              text(x+20*dx,y-1000*dy,'H=2.2')
 
% xL = xlim;
% yL = ylim;
% line([0 0], yL);  %y-axis
% line(xL, [0 0]);  %x-axis
           
    print('-depsc2', 'exLVline01.eps'); 
    print('-dpdf', 'exLVline01.pdf','-b0');



