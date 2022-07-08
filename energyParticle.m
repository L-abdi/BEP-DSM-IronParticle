%% Calculate iron particle's internal energy
function ep = energyParticle(mFe, mFeO, Tp)
MW_Fe = 55.845e-3;                                      % Pure iron molar weight, [kg/mol]
MW_FeO = 71.844e-3;                                     % Iron-oxide molar weight, [kg/mol]

ep = mFe / MW_Fe * hFe(Tp) + mFeO / MW_FeO * hFeO(Tp);
return
