function mdot_R = kineticRate(mFe,mFeO,ep)
rhoFeO  = 5745.0;
gamma_FeO_O2 = 4.490418390689651;
k0FeO = 2.669800000000000e-04;
TaFeO = 20319;
rFe = radiusFe(mFe,mFeO,ep);
rp = radiusFeO(mFe,mFeO,ep);
Tp = temperatureParticle(mFe,mFeO,ep);
Ap = areaParticle(mFe,mFeO,ep);

XFeO = rp - rFe;

dXdt = (rp-XFeO)./...
        (rp .* XFeO) .* k0FeO * exp(-TaFeO ./ Tp);
    
mdot_R = gamma_FeO_O2 * rhoFeO * Ap .* dXdt;

return