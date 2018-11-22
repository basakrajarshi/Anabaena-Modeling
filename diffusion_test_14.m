function y_out = diffusion_test()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% May 2018, Orit Peleg, orit.peleg@colorado.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
dbstop if error;
set(0,'Defaultlinelinewidth',3.5, 'DefaultlineMarkerSize',12,...
    'DefaultTextFontSize',5, 'DefaultAxesFontSize',18);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The agent at position (0,0) produces a chemical signal that spreads out following a 
%diffusion equation. 

%The output of the code is a movie showing the concentration of the 
%chemical signal as a function of time in the arena. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define   
A_1=0.2;
A_2=0.4;%initial chemical concentration per active agent 

D_1=0.5;
D_2=0.7;% diffusion coefficient
delta_t_1 = 0.1;
delta_t_2 = 0.1;% time integration constant (time discretization)
delta_x_1 = 0.1;
delta_x_2 = 0.1% size of a cell in the simulation (space discretization)
B_1=0;
B_2=0;% ignore for now

SAVE_MOVIE = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        
t_array_1=0:delta_t_1:0.5;
t_array_2=0:delta_t_2:0.5;
min_x_1 = -5; max_x_1=5;
min_x_2 = -5; max_x_2=5;

source_and_sink_propagation_c1(A_1,B_1,D_1,t_array_1,delta_x_1,min_x_1,max_x_1,SAVE_MOVIE);
y_out=0; 

source_and_sink_propagation_c2(A_2,B_2,D_2,t_array_2,delta_x_2,min_x_2,max_x_2,SAVE_MOVIE);
y_out=0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function source_and_sink_propagation_c1(A,B,D,t_array,delta_x,min_x,max_x,SAVE_MOVIE)

source_x1 = -1.0;
source_y1 = 0.0;
sink_x2 = 1.0;
sink_y2 = 0.0;

[x,y] = meshgrid(min_x:delta_x:max_x,min_x:delta_x:max_x);

if SAVE_MOVIE==1
    fig_densities=figure();
    loops = length(2:length(t_array));
    clear F
    F(loops) = struct('cdata',[],'colormap',[]);
end

for t_i=2:length(t_array)
    t=t_array(t_i);
    %c = 0.1 * ones(size(x));
    c = 0.05 * ones(size(x));

    curr_t = t;
    curr_c_1 =(A/(curr_t^0.5))*exp(-(((x-source_x1).^2)+((y-source_y1).^2))./(4*D*curr_t))+B;
    curr_c_2 =(A/(curr_t^0.5))*exp(-(((x-sink_x2).^2)+((y-sink_y2).^2))./(4*D*curr_t))+B;
    c = c + curr_c_1 - curr_c_2;
    %c = c + curr_c_1;
    
    if (t_i>1) && SAVE_MOVIE==1
        figure(fig_densities);
        surf(x,y,c,'EdgeColor','none'); view(2); title(['t=',num2str(t)]); xlim([min_x max_x]); ylim([min_x max_x]);
        daspect([1 1 1]); colorbar; caxis([0 A]); hold on;
        
        pause(0.001);
        F(t_i-1) = getframe(fig_densities);
        pause(0.001);
        
        pause(0.1);
    end
    
end

if SAVE_MOVIE==1
    v = VideoWriter(['A',num2str(A),'B',num2str(A),'D',num2str(A),'test.mp4'],'MPEG-4');
    open(v);
    writeVideo(v,F);
    close(v);
end

end

function source_and_sink_propagation_c2(A,B,D,t_array,delta_x,min_x,max_x,SAVE_MOVIE)

source_x1 = 1.0;
source_y1 = 0.0;
sink_x2 = -1.0;
sink_y2 = 0.0;

[x,y] = meshgrid(min_x:delta_x:max_x,min_x:delta_x:max_x);

if SAVE_MOVIE==1
    fig_densities=figure();
    loops = length(2:length(t_array));
    clear F
    F(loops) = struct('cdata',[],'colormap',[]);
end

for t_i=2:length(t_array)
    t=t_array(t_i);
    %c = 0.1 * ones(size(x));
    c = 0.05 * ones(size(x));

    curr_t = t;
    curr_c_1 =(A/(curr_t^0.5))*exp(-(((x-source_x1).^2)+((y-source_y1).^2))./(4*D*curr_t))+B;
    curr_c_2 =(A/(curr_t^0.5))*exp(-(((x-sink_x2).^2)+((y-sink_y2).^2))./(4*D*curr_t))+B;
    c = c + curr_c_1 - curr_c_2;
    %c = c + curr_c_1;
    
    if (t_i>1) && SAVE_MOVIE==1
        figure(fig_densities);
        surf(x,y,c,'EdgeColor','none'); view(2); title(['t=',num2str(t)]); xlim([min_x max_x]); ylim([min_x max_x]);
        daspect([1 1 1]); colorbar; caxis([0 A]); hold on;
        
        pause(0.001);
        F(t_i-1) = getframe(fig_densities);
        pause(0.001);
        
        pause(0.1);
    end
    
end

if SAVE_MOVIE==1
    v = VideoWriter(['A',num2str(A),'B',num2str(A),'D',num2str(A),'test.mp4'],'MPEG-4');
    open(v);
    writeVideo(v,F);
    close(v);
end

end

