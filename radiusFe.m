function rFe = radiusFe(mFe, mFeO, ep)
rhoFe = 7874.0; 

rFe = ( 0.75 * mFe / (pi * rhoFe) ).^(1/3);

end