function y_out_1 = diffusion_test()
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
A_1=0.1; %initial chemical concentration per active agent for 1 
A_2=0.1; %initial chemical concentration per active agent for 2

D_1=0.5; % diffusion coefficient for 1
D_2=0.5; % diffusion coefficient for 2
delta_t_1 = 0.1; % time integration constant (time discretization) for 1
delta_t_2 = 0.1; % time integration constant (time discretization) for 2
delta_x_1 = 0.1; % size of a cell in the simulation (space discretization)
delta_x_2 = 0.1; % size of a cell in the simulation (space discretization)
B_1=0; % ignore for now
B_2=0; % ignore for now

SAVE_MOVIE_1 = 1;
SAVE_MOVIE_2 = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        
t_array_1=0:delta_t_1:0.3; 
t_array_2=0:delta_t_2:0.3;
min_x_1 = -5; max_x_1=5;
min_x_2 = -5; max_x_2=5;

source_propagation(A_1,A_2,B_1,B_2,D_1,D_2,t_array_1,t_array_2,delta_x_1,delta_x_2,min_x_1,min_x_2,max_x_1,max_x_2,SAVE_MOVIE_1,SAVE_MOVIE_2);
y_out_1=0; 


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function source_propagation(A_1,A_2,B_1,B_2,D_1,D_2,t_array_1,t_array_2,delta_x_1,delta_x_2,min_x_1,min_x_2,max_x_1,max_x_2,SAVE_MOVIE_1,SAVE_MOVIE_2)

source_x_1 = 1;
source_x_2 = 1;
source_y_1 = 1;
source_y_2 = 1;

[x1,y1] = meshgrid(min_x_1:delta_x_1:max_x_1,min_x_1:delta_x_1:max_x_1);
[x2,y2] = meshgrid(min_x_2:delta_x_2:max_x_2,min_x_2:delta_x_2:max_x_2);

if SAVE_MOVIE_1==1 && SAVE_MOVIE_2==1
    fig_densities_1=figure();
    fig_densities_2=figure();
    loops_1 = length(2:length(t_array_1));
    loops_2 = length(2:length(t_array_2));
    clear F1
    clear F2
    F1(loops_1) = struct('cdata',[],'colormap',[]);
    F2(loops_2) = struct('cdata',[],'colormap',[]);
end

for t_i_1=2:length(t_array_1)
    t_1=t_array_1(t_i_1);
    c_1 = zeros(size(x1));

    curr_t_1 = t_1;
    curr_c_1 =(A_1/(curr_t_1^0.5))*exp(-(((x1-source_x_1).^2)+((y1-source_y_1).^2))./(4*D_1*curr_t_1))+B_1;
    c_1 = c_1 + curr_c_1;

    
    if (t_i_1>1) && SAVE_MOVIE_1==1  
        figure(fig_densities_1);
        surf(x1,y1,c_1,'EdgeColor','none'); view(2); title(['t=',num2str(t_1)]); xlim([min_x_1 max_x_1]); ylim([min_x_1 max_x_1]);
        daspect([1 1 1]); colorbar; caxis([0 A_1]); hold on;
        
        pause(0.001);
        F_1(t_i_1) = getframe(fig_densities_1);
        pause(0.001);
        
        pause(0.1);
        
    end
    
end

for t_i_2=2:length(t_array_2)
    t_2=t_array_2(t_i_2);
    c_2 = zeros(size(x2));

    curr_t_2 = t_2;
    curr_c_2 =(A_2/(curr_t_2^0.5))*exp(-(((x2-source_x_2).^2)+((y2-source_y_2).^2))./(4*D_2*curr_t_2))+B_2;
    c_2 = c_2 + curr_c_2;
    
    if (t_i_2>1) && SAVE_MOVIE_2==1
        figure(fig_densities_2);
        surf(x2,y2,c_2,'EdgeColor','none'); view(2); title(['t=',num2str(t_2)]); xlim([min_x_2 max_x_2]); ylim([min_x_2 max_x_2]);
        daspect([1 1 1]); colorbar; caxis([0 A_2]); hold on;
        
        pause(0.001);
        F_2(t_i_2) = getframe(fig_densities_2);
        pause(0.001);
        
        pause(0.1);
    end
end



if SAVE_MOVIE_1==1 && SAVE_MOVIE_2==1
    v_1 = VideoWriter(['A',num2str(A_1),'B',num2str(A_1),'D',num2str(A_1),'test.mp4'],'MPEG-4');
    v_2 = VideoWriter(['A',num2str(A_2),'B',num2str(A_2),'D',num2str(A_2),'test.mp4'],'MPEG-4');
    open(v_1);
    open(v_2);
    writeVideo(v_1,F_1);
    writeVideo(v_2,F_2);
    close(v_1);
    close(v_2);
end

end


