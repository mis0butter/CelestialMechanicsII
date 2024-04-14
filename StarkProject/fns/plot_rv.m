function plot_rv( rv_hist, plot_title ) 

    if ~exist('plot_title', 'var') 
        plot_title = 'Cartesian Trajectory'; 
    end 

%     figure() ; 
        plot3( rv_hist(:,1), rv_hist(:,2), rv_hist(:,3) ) ; 
        hold on ; grid on ; 
        scatter3( rv_hist(:,1), rv_hist(:,2), rv_hist(:,3) ) ; 
        view(30, 30) 
        xlabel('x (LU)') ; ylabel('y (LU)') ; zlabel('z (LU)') ; 
        title(plot_title) ; 

end 



