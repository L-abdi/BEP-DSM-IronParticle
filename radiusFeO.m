function rFeO = radiusFeO(mFe, mFeO, ep)
rhoFeO  = 5745.0;
rhoFe = 7874.0; 

rFeO = ( 0.75 * mFeO / (pi * rhoFeO) ).^(1/3) + ( 0.75 * mFe / (pi * rhoFe) ).^(1/3);

return
