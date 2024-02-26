function eps = eps_fn( u, uprime, mu ) 

    r   = dot( u, u ) ; 
    eps = 2 * dot( uprime, uprime ) / r - mu / r ; 

end 