function Lu = Lu_fn(u) 

    Lu = [  u(1)    -u(2)   -u(3)   u(4)  ; 
            u(2)    u(1)    -u(4)   -u(3) ; 
            u(3)    u(4)    u(1)    u(2)  ; 
            u(4)    -u(3)   u(2)    -u(1) ] ; 

end 