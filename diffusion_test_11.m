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

D_1=1.0; % diffusion coefficient for 1
D_2=1.0; % diffusion coefficient for 2
delta_t_1 = 0.1; % time integration constant (time discretization) for 1
delta_t_2 = 0.1; % time integration constant (time discretization) for 2
delta_x_1 = 0.1; % size of a cell in the simulation (space discretization)
delta_x_2 = 0.1; % size of a cell in the simulation (space discretization)
B_1=0; % ignore for now
B_2=0; % ignore for now

SAVE_MOVIE_1 = 1;
SAVE_MOVIE_2 = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        
t_array_1=0:delta_t_1:1.0; 
t_array_2=0:delta_t_2:1.0;
min_x_1 = -5; max_x_1=5;
min_x_2 = -5; max_x_2=5;

source_propagation(A_1,B_1,D_1,t_array_1,delta_x_1,min_x_1,max_x_1,SAVE_MOVIE_1);
y_out_1=0; 

source_propagation(A_2,B_2,D_2,t_array_2,delta_x_2,min_x_2,max_x_2,SAVE_MOVIE_2);
y_out_2=0; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function source_propagation(A,B,D,t_array,delta_x,min_x,max_x,SAVE_MOVIE)

source_x = 1;
source_y = 1;

[x,y] = meshgrid(min_x:delta_x:max_x,min_x:delta_x:max_x);

if SAVE_MOVIE==1
    fig_densities=figure();
    loops = length(2:length(t_array));
    clear F
    F(loops) = struct('cdata',[],'colormap',[]);
end

for t_i=2:length(t_array)
    t=t_array(t_i);
    c = zeros(size(x));

    curr_t = t;
    curr_c =(A/(curr_t^0.5))*exp(-(((x-source_x).^2)+((y-source_y).^2))./(4*D*curr_t))+B;
    c = c + curr_c;

    
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


