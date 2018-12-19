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
A_1=0.1;
A_2=0.1;
A_3=0.1;%initial chemical concentration per active agent 

D_1=0.9;
D_2=0.7;
D_3=0.5;% diffusion coefficient

delta_t = 0.1;
%delta_t_2 = 0.1;
%delta_t_3 = 0.1;% time integration constant (time discretization)

delta_x = 0.1;
%delta_x_2 = 0.1;
%delta_x_3 = 0.1;% size of a cell in the simulation (space discretization)

B_1=0;
B_2=0;
B_3=0;% ignore for now

SAVE_MOVIE = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        
t_array=0:delta_t:1.0;
%t_array_2=0:delta_t_2:1.5;
%t_array_3=0:delta_t_3:1.5;

min_x = -5; max_x=5;
%min_x_2 = -5; max_x_2=5;
%min_x_3 = -5; max_x_3=5;

source_and_sink_propagation_c(A_1,A_2,A_3,B_1,B_2,B_3,D_1,D_2,D_3,t_array,delta_x,min_x,max_x,SAVE_MOVIE);
y_out=0; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function source_and_sink_propagation_c(A_1,A_2,A_3,B_1,B_2,B_3,D_1,D_2,D_3,t_array,delta_x,min_x,max_x,SAVE_MOVIE)

% First simulation
%source_x1_1 = -4.0;
%source_y1_1 = 0.0;
%source_x2_1 = -2.0;
%source_y2_1 = 0.0
%source_x3_1 = 0.0
%source_y3_1 = 0.0
%source_x4_1 = 2.0
%source_y4_1 = 0.0
%sink_x5_1 = 4.0;
%sink_y5_1 = 0.0;

source_x_1 = [4]
source_y_1 = [0]
sink_x_1 = [-4 -2 0 2]
sink_y_1 = [0 0 0 0]


% Second simulation
%source_x1_2 = -4.0;
%source_y1_2 = 0.0;
%source_x2_2 = -2.0
%source_y2_2 = 0.0
%source_x3_2 = 0.0
%source_y3_2 = 0.0
%sink_x4_2 = 2.0
%sink_y4_2 = 0.0
%source_x5_2 = 4.0;
%source_y5_2 = 0.0;

source_x_2 = [2]
source_y_2 = [0]
sink_x_2 = [-4 -2 0 4]
sink_y_2 = [0 0 0 0]


% Third simulation
%source_x1_3 = -4.0;
%source_y1_3 = 0.0;
%source_x2_3 = -2.0
%source_y2_3 = 0.0
%sink_x3_3 = 0.0
%sink_y3_3 = 0.0
%source_x4_3 = 2.0
%source_y4_3 = 0.0
%source_x5_3 = 4.0;
%source_y5_3 = 0.0;

source_x_3 = [0]
source_y_3 = [0]
sink_x_3 = [-4 -2 2 4]
sink_y_3 = [0 0 0 0]

[x,y] = meshgrid(min_x:delta_x:max_x,min_x:delta_x:max_x);

if SAVE_MOVIE==1
    %fig_densities_1=figure();
    %loops_1 = length(2:length(t_array_1));
    %clear F_1
    %F_1(loops_1) = struct('cdata',[],'colormap',[]);
    
    fig_densities=figure();
    loops = length(2:length(t_array));
    clear F
    F(loops) = struct('cdata',[],'colormap',[]);
end

for t_i=2:length(t_array)
    t=t_array(t_i);
    
    %c = 0.1 * ones(size(x));
    c_1 = (A_1/2) * ones(size(x));
    c_2 = (A_2/2) * ones(size(x));
    c_3 = (A_3/2) * ones(size(x));
    
    for j=1:length(sink_x_1)
        curr_t = t;
        curr_c_1_1 =(A_1/(curr_t^0.5))*exp(-(((x-sink_x_1(j)).^2)+((y-sink_y_1(j)).^2))./(4*D_1*curr_t))+B_1;
        %curr_c_2_1 =(A_1/(curr_t^0.5))*exp(-(((x-sink_x_1(1)).^2)+((y-sink_y_1(1)).^2))./(4*D_1*curr_t))+B_1;
        %curr_c_2_1 =(A_1/(curr_t^0.5))*exp(-(((x-sink_x(1)).^2)+((y-sink_y(1)).^2))./(4*D*curr_t))+B;
        %c = c + curr_c_1 + curr_c_2;
        %c = c + curr_c_1_1 - curr_c_2_1;
        c_1 = c_1 - curr_c_1_1
    end
    
    curr_c_2_1 =(A_1/(curr_t^0.5))*exp(-(((x-source_x_1(1)).^2)+((y-source_y_1(1)).^2))./(4*D_1*curr_t))+B_1;
    c_1 = c_1 + curr_c_2_1;
    
    for k=1:length(sink_x_2)
        curr_t = t;
        curr_c_1_2 =(A_2/(curr_t^0.5))*exp(-(((x-sink_x_2(k)).^2)+((y-sink_y_2(k)).^2))./(4*D_2*curr_t))+B_2;
        %curr_c_2_2 =(A_2/(curr_t^0.5))*exp(-(((x-sink_x_2(1)).^2)+((y-sink_y_2(1)).^2))./(4*D_2*curr_t))+B_2;
        %curr_c_2_1 =(A_1/(curr_t^0.5))*exp(-(((x-sink_x(1)).^2)+((y-sink_y(1)).^2))./(4*D*curr_t))+B;
        %c = c + curr_c_1 + curr_c_2;
        %c = c + curr_c_1_1 - curr_c_2_1;
        c_2 = c_2 - curr_c_1_2
    end
    
    curr_c_2_2 =(A_2/(curr_t^0.5))*exp(-(((x-source_x_2(1)).^2)+((y-source_y_2(1)).^2))./(4*D_2*curr_t))+B_2;
    c_2 = c_2 + curr_c_2_2;
    
    for l=1:length(sink_x_3)
        curr_t = t;
        curr_c_1_3 =(A_3/(curr_t^0.5))*exp(-(((x-sink_x_3(l)).^2)+((y-sink_y_3(l)).^2))./(4*D_3*curr_t))+B_3;
        %curr_c_2_3 =(A_3/(curr_t^0.5))*exp(-(((x-sink_x_3(1)).^2)+((y-sink_y_3(1)).^2))./(4*D_3*curr_t))+B_3;
        %curr_c_2_1 =(A_1/(curr_t^0.5))*exp(-(((x-sink_x(1)).^2)+((y-sink_y(1)).^2))./(4*D*curr_t))+B;
        %c = c + curr_c_1 + curr_c_2;
        %c = c + curr_c_1_1 - curr_c_2_1;
        c_3 = c_3 - curr_c_1_3
    end
    
    curr_c_2_3 =(A_3/(curr_t^0.5))*exp(-(((x-source_x_3(1)).^2)+((y-source_y_3(1)).^2))./(4*D_3*curr_t))+B_3;
    c_3 = c_3 + curr_c_2_3;

    
    
    if (t_i>1) && SAVE_MOVIE==1
        %figure(fig_densities_1);
        subplot(1,3,1);
        surf(x,y,c_1,'EdgeColor','none'); view(2); title(['t=',num2str(t)]); xlim([min_x max_x]); ylim([min_x max_x]);
        daspect([1 1 1]); colorbar; caxis([0 A_1]); hold on;
        
        pause(0.001);
        F(t_i-1) = getframe(fig_densities);
        pause(0.001);
        
        pause(0.1);
        
        
        
        %figure(fig_densities_2);
        subplot(1,3,2);
        surf(x,y,c_2,'EdgeColor','none'); view(2); title(['t=',num2str(t)]); xlim([min_x max_x]); ylim([min_x max_x]);
        daspect([1 1 1]); colorbar; caxis([0 A_2]); hold on;
        
        pause(0.001);
        F(t_i-1) = getframe(fig_densities);
        pause(0.001);
        
        pause(0.1);
        
        
        
        %figure(fig_densities_3);
        subplot(1,3,3);
        surf(x,y,c_3,'EdgeColor','none'); view(2); title(['t=',num2str(t)]); xlim([min_x max_x]); ylim([min_x max_x]);
        daspect([1 1 1]); colorbar; caxis([0 A_3]); hold on;
        
        pause(0.001);
        F(t_i-1) = getframe(fig_densities);
        pause(0.001);
        
        pause(0.1);
    end
    
end

if SAVE_MOVIE==1
    v_1 = VideoWriter(['A',num2str(A_1),'B',num2str(A_1),'D',num2str(A_1),'test.mp4'],'MPEG-4');
    open(v_1);
    writeVideo(v_1,F);
    close(v_1);
    
    v_2 = VideoWriter(['A',num2str(A_2),'B',num2str(A_2),'D',num2str(A_2),'test.mp4'],'MPEG-4');
    open(v_2);
    writeVideo(v_2,F);
    close(v_2);
    
    v_3 = VideoWriter(['A',num2str(A_3),'B',num2str(A_3),'D',num2str(A_3),'test.mp4'],'MPEG-4');
    open(v_3);
    writeVideo(v_3,F);
    close(v_3);
end

end