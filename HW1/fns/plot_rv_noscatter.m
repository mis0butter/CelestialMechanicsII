function plot_rv_noscatter( rv_hist, plot_title ) 

    if ~exist('plot_title', 'var') 
        plot_title = 'Cartesian Trajectory'; 
    end 

%     figure() ; 
        line_rv = plot3( rv_hist(:,1), rv_hist(:,2), rv_hist(:,3) ) ; 
        color   = get(line_rv, 'color') ; 
        set(line_rv, 'color', [ color 0.4 ]) ; 
        
        hold on ; grid on ; 
%         scatter3( rv_hist(:,1), rv_hist(:,2), rv_hist(:,3) ) ; 
        view(30, 30) 
        xlabel('x (LU)') ; ylabel('y (LU)') ; zlabel('z (LU)') ; 
        title(plot_title) ; 

end 



