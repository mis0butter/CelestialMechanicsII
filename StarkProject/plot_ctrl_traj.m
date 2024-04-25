function f = plot_ctrl_traj( X_hist, acc_hist, alphaSun, order ) 

acc_hist_norm = [ ] ; 
for i = 1 : size(acc_hist, 1) 
    acc_hist_norm = [ acc_hist_norm ; norm(acc_hist(i,:)) ] ; 
end 

f = figure() ; 
    subplot(3,1,1) ; 
        semilogy( acc_hist_norm, 'b' ) ; hold on ; grid on ; 
        scatter( [ 1 : length(acc_hist) ], acc_hist_norm, 'b' ) ; 
        xlabel( {'Segment'; ''} ) ; 
        ylabel( 'Acc mag (LU/TU^2)' ) ; 
        title( sprintf( 'alpha = %.2g, order = %d', alphaSun, order ) ) ;         
    subplot(3,1,2:3) ; 
        scatter3( 0,0,0, '+k' ) ; hold on ; grid on ; 
        plot3( X_hist(:,1), X_hist(:,2), X_hist(:,3), 'b' ) ; 
        scatter3( X_hist(:,1), X_hist(:,2), X_hist(:,3), 'b' )
        scatter3( X_hist(1,1), X_hist(1,2), X_hist(1,3), '^g', 'filled' ) ; 
        scatter3( X_hist(end,1), X_hist(end,2), X_hist(end,3), 'diamondr', 'filled' ) ; 
        view(0,90)
        xlabel('x') ; ylabel('y') ; zlabel('z') ; 
        title( sprintf( 'trajectory' ) ) ; 
        
end 

