function f = plot_dH_traj( X_hist, mu, accuracy_hist, alphaSun, order ) 

oe_hist = [ ] ; E_hist = [ ] ; 
for i = 1 : size(X_hist, 1) 
    
    rv = X_hist( i, 1:6 ) ; 
    oe = rv2oe( rv, mu ) ;  
    
    e  = oe(2) ; 
    nu = oe(6) ; 
    E  = nu2E_fn( nu, e ) ; 
    
    oe_hist = [ oe_hist ; oe' ] ; 
    E_hist  = [ E_hist ; E ] ; 
    
end 

% make E_hist propagate forward 
for i = 1 : length(E_hist) - 1 
    
    if E_hist(i) > E_hist(i+1) 
        E_hist(i+1) = pi + ( pi - E_hist(i+1) ) ; 
    end 
    
end 

f = figure() ; 
    subplot(3,1,1) ; 
        semilogy( E_hist, accuracy_hist(:,1), 'b' ) ; hold on ; grid on ; 
        scatter( E_hist, accuracy_hist(:,1), 'b' ) ; 
        xlabel( {'Eccentric anomaly E'; ''} ) ; 
        ylabel( 'delta H' ) ; 
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


