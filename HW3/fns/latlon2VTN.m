function VTN = latlon2VTN( Tmag, lat, lon ) 

    V = Tmag * cos( lon ) * cos( lat ) ; 
    T = Tmag * sin( lon ) * cos( lat ) ; 
    N = Tmag * sin( lat ) ; 
    
    VTN = [ V T N ]' ; 

end 