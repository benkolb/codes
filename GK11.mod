// Dynare model file to calculate the GK model
// Original code created by Peter Karadi, July 2010
// Modified by Benedikt Kolb, Feb. 2014
// changes include:
// - correction of equation 6, "Value of banks' net wealth"
// - analytic calculation of steady state values in steady_state_model (also, no more initval)
// - analytic calculation of parameters calibrated on steady-state values as "hashed parameters" in model block (b del lam chi omg); therefore, deep parameters can now be estimated
// - deleting descriptive variables not required for simulation/estimation (Welf, VMPK, W, Keff)
// - changing variable and parameter names
// The plotted IRFs replicate the black lines in Figure 2 of the paper (all in % deviations from SS).

parameters 
bet sig hab 
varphi 
zet tet alf gov acI eps gam gmP kpP kpY 
rho_a rho_i rho_k rho_g rho_p 
sig_a sig_i sig_k sig_g sig_p sig_n 
kap tau
// steady states 
LSS GSS FSS ISS DELSS SPRSS
;
 
var 
Y Ym K L I C G Q PI PIstar RHO LAM Rk R N Ne Nn NU ETA PHI Z X Pm U MU D F H Rn DEL In
// exog. shock processes
EXA % technology (TFP)
EXK % capital efficiency
EXG % government spending
// for plotting
SPR % premium / spread Rk-R
;
 
varexo % shock to ...
e_a    % ... technology (TFP) 
e_k    % ... capital efficiency
e_g    % ... government spending
e_n    % ... net worth
e_i    % ... investment technology
;
 
bet = 0.99;
sig = 1;
hab = 0.815;
varphi = 0.276;
zet = 7.2;
tet = 0.97155955;
alf = 0.33;
gov = 0.2;     % SS government spending over GDP
acI = 1.728;
eps = 4.167;
gam = 0.779;
gmP = 0.241;
kpP = 1.5;
kpY = -0.125;
rho_i = 0;
rho_k = 0.66;
sig_k = 0.05;
rho_a = 0.95;
sig_a = 0.01;
rho_g = 0.95;
sig_g = 0.01;
sig_n = 0.01;
sig_i = 0.01;
rho_p = 0.66;
sig_p = 0.072;
kap = 10;
tau = 0.001;
GSS = 0.16975710;
LSS = 0.33333333;
FSS = 4;          % leverage
ISS = 0.14153927;
DELSS = 0.025; 
SPRSS = 0.01/4;   % spread Rk-R
% omg = 0.00222778;
% lam = 0.38149499;
% chi = 3.41080850;
% b   = 0.03760101;
% del = 0.02041451;

model;
#KSS   = (alf*(eps-1)/eps/(1/bet+SPRSS-1+DELSS))^(1/(1-alf))*LSS;
#YSS   = KSS^alf*LSS^(1-alf);
#b     = alf*(eps-1)/eps*YSS/KSS;
#del   = DELSS - b/(1+zet);
#omg   = (1-tet*(SPRSS*FSS + 1/bet))*KSS/FSS/KSS;
#ZSS   = (1-omg*FSS)/tet;
#ETASS = (1-tet)/(1-tet*bet*ZSS);
#NUSS  = (1-tet)*bet*SPRSS/(1-bet*tet*ZSS);
#lam   = NUSS + ETASS/FSS;
#RHOSS = (1-bet*hab)*((1-hab)*(YSS-DELSS*KSS-gov*YSS))^(-sig);
#chi   = RHOSS*(eps-1)/eps*(1-alf)*YSS/LSS^(1+varphi);


//Households
//1. Marginal utility of consumption
exp(RHO) =  (exp(C)-hab*exp(C(-1)))^(-sig) - bet*hab*(exp(C(+1))-hab*exp(C))^(-sig);

//2. Euler equation
bet*exp(R)*exp(LAM(+1)) = 1;

//3. Stochastic discount rate
exp(LAM) = exp(RHO)/exp(RHO(-1));

//4. Labor market equilibrium
chi*exp(L)^varphi =  exp(RHO)*exp(Pm)*(1-alf)*exp(Y)/exp(L);

//Financial Intermediaries
//5. Value of banks' capital
exp(NU) = (1-tet)*bet*exp(LAM(+1))*(exp(Rk(+1))-exp(R))+bet*exp(LAM(+1))*tet*exp(X(+1))*exp(NU(+1));

//6. Value of banks' net wealth
exp(ETA) = bet*exp(LAM(+1))*( (1-tet)*exp(R) + tet*exp(Z(+1))*exp(ETA(+1)) );
//NOTE: The corresponding equation in Peter Karadi's code is:
// exp(ETA) =  (1-tet) + bet*exp(LAM(+1))*tet*exp(Z(+1))*exp(ETA(+1));

//7. Optimal leverage
exp(PHI) = exp(ETA)/(lam-exp(NU));

//8. Growth rate of banks' capital
exp(Z) = (exp(Rk)-exp(R(-1)))*exp(PHI(-1)) + exp(R(-1));

//9. Growth rate of banks' net wealth
exp(X) = exp(PHI)/exp(PHI(-1))*exp(Z);

//Aggregate capital, net worth
//10. Aggregate capital
exp(Q)*exp(K) = exp(PHI)*exp(N);

//11. Banks' net worth
exp(N) = exp(Ne) + exp(Nn);

//12. Existing banks' net worth accumulation
exp(Ne) = tet*exp(Z)*exp(N(-1))*exp(-sig_n*e_n);

//13. New banks' net worth
exp(Nn) = omg*exp(Q)*exp(EXK)*exp(K(-1));

//Final goods producer
//14. Return to capital
exp(Rk) =  (exp(Pm)*alf*exp(Ym)/exp(K(-1))+exp(EXK)*(exp(Q)-exp(DEL)))/exp(Q(-1));

//15. Production function
exp(Ym) = exp(EXA)*(exp(EXK)*exp(U)*exp(K(-1)))^alf*exp(L)^(1-alf);

//Capital Goods Producer
//16. Optimal investment decision
exp(Q) = 1 + acI/2*((In+ISS)/(In(-1)+ISS)-1)^2+acI*((In+ISS)/(In(-1)+ISS)-1)*(In+ISS)/(In(-1)+ISS)-bet*exp(LAM(+1))*acI*((In(+1)+ISS)/(In+ISS)-1)*((In(+1)+ISS)/(In+ISS))^2;

//17. Depreciation rate
exp(DEL) = del + b/(1+zet)*exp(U)^(1+zet);

//18. Optimal capacity utilization rate
exp(Pm)*alf*exp(Ym)/exp(U) = b*exp(U)^zet*exp(EXK)*exp(K(-1));

//19. Net investment
In = exp(I) - exp(DEL)*exp(EXK)*exp(K(-1));

//20. Capital accumulation equation
exp(K) = exp(EXK)*exp(K(-1)) + In; 

//21. Government consumption
exp(G) = GSS*exp(EXG);

//Equilibrium
//22. Aggregate resource constraint
exp(Y) = exp(C)+exp(G)+exp(I)+acI/2*((In+ISS)/(In(-1)+ISS)-1)^2*(In+ISS);

//23. Wholesale, retail output
exp(Ym) = exp(Y)*exp(D);

//24. Price dispersion
exp(D) = gam*exp(D(-1))*exp(PI(-1))^(-kap*eps)*exp(PI)^eps+(1-gam)*((1-gam*exp(PI(-1))^(kap*(1-gam))*exp(PI)^(gam-1))/(1-gam))^(-eps/(1-gam));

//25. Markup
exp(MU) = 1/exp(Pm);

//26. Optimal price choice I
exp(F) =  exp(Y)*exp(Pm)+bet*gam*exp(LAM(+1))*exp(PI(+1))^eps*(exp(PI))^(-eps*kap)*exp(F(+1));

//27. Optimal price choice II
exp(H) = exp(Y)+bet*gam*exp(LAM(+1))*exp(PI(+1))^(eps-1)*exp(PI)^(kap*(1-eps))*exp(H(+1));

//28. Optimal price choice III
exp(PIstar) = eps/(eps-1)*exp(F)/exp(H)*exp(PI);

//29. Price index
(exp(PI))^(1-eps) = gam*exp(PI(-1))^(kap*(1-eps))+(1-gam)*(exp(PIstar))^(1-eps);

//30. Fisher equation
exp(Rn) = exp(R)*exp(PI(+1));

//31. Interest rate rule
exp(Rn)      =   exp(Rn(-1))^rho_i*((1/bet)*exp(PI)^kpP*(exp(MU)/(eps/(eps-1)))^(kpY))^(1-rho_i)*exp(sig_i*e_i);

//Shocks
//32. TFP shock
EXA = rho_a*EXA(-1) - sig_a*e_a;

//33. Capital quality shock
EXK = rho_k*EXK(-1) - sig_k*e_k;

//34. Government consumption shock
EXG = rho_g*EXG(-1) - sig_g*e_g;

//For plotting:
//35. Premium
exp(SPR) = exp(Rk(+1))/exp(R);
end;


steady_state_model; % all not mentioned are zero
L   = log(LSS);
PHI = log(FSS);
R   = log(1/bet);
Pm  = log((eps-1)/eps);
MU  = log(1/exp(Pm));
DEL = log(DELSS);
Rk  = log(exp(R)+SPRSS);
Rn  = R;
K   = log((alf*exp(Pm)/(exp(Rk)-1+DELSS))^(1/(1-alf))*LSS); 
Y   = log(exp(K)^alf*exp(L)^(1-alf));
Ym  = Y;
G   = log(gov*exp(Y));
I   = log(DELSS*exp(K));
C   = log(exp(Y)-exp(I)-exp(G));
RHO = log((1-bet*hab)*((1-hab)*exp(C))^(-sig));
N   = log(exp(K)/exp(PHI));
Z   = log(SPRSS*FSS + 1/bet);
X   = Z;
Ne  = log(tet*exp(Z)*exp(N));
Nn  = log(exp(N)-exp(Ne));
NU  = log(((1-tet)*bet*SPRSS)/(1-bet*tet*exp(X)));
ETA = log((1-tet)/(1-tet*bet*exp(Z)));
F   = log(exp(Y)*exp(Pm)/(1-bet*gam));
H   = log(exp(Y)/(1-bet*gam));
SPR = log(exp(Rk)/exp(R));
end;


// SIMULATION.

resid;
check;
steady;
model_diagnostics;
options_.pruning        = 1;
options_.relative_irf   = 1;
options_.nocorr         = 1;
options_.nomoments      = 1;
options_.noprint        = 1;
options_.nograph        = 0;
options_.simul_seed     = 42;



shocks;
var e_a = 0;
var e_k = 1; 
var e_g = 0; 
var e_n = 0; 
var e_i = 0;
end;

stoch_simul(order=1, periods=2000, irf=40) EXK R SPR Y C I K L Q N PI Rn;

