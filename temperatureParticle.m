%% Calculate iron particle's temperature
function Tp = temperatureParticle(mFe, mFeO, ep)

Tp = 300.0;
Tpn = 999999.9;
itr = 1;
diff = 1.0e8;
while ( abs(diff) > 1.0e-1 && itr < 40)
    f = ep - energyParticle(mFe, mFeO, Tp);
    dfdT = -cpParticle(mFe, mFeO, Tp);
    Tpn = Tp - f / dfdT;
    diff = Tp - Tpn;
    Tp = Tpn;
    itr = itr + 1;
end

if (itr >= 40)
    fprintf('N-R Root finding wrong Tp = %g\n',Tp)
end

return