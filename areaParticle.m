function Ap = areaParticle(mFe,mFeO,ep)

rp = radiusFeO(mFe,mFeO,ep);

Ap =  4.0 * pi *rp.^2.0;

return