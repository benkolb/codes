// Dynare model file to calculate the GK model
// Original code created by Peter Karadi, July 2010
// Modified by Benedikt Kolb, Feb. 2014
// changes include:
// - correction of equation 6, "Value of banks' net wealth"
// - analytic calculation of steady state values in steady_state_model (also, no more initval)
// - analytic calculation of parameters calibrated on steady-state values as "hashed parameters" in model block (utc del lam chi omg); therefore, deep parameters can now be estimated
// - deleting descriptive variables not required for simulation/estimation (Welf, VMPK, W, Keff)
// - changing variable and parameter names
// The plotted IRFs replicate the black lines in Figure 2 of the paper (all in % deviations from SS).

parameters  % explanations below
bet sig hab fri zet tet alf gov acI eps gam kpP kpY kap
// exogenous-process parameters
rho_a rho_i rho_k rho_g rho_p 
sig_a sig_i sig_k sig_g sig_p sig_n 
// steady states 
LSS GSS FSS ISS DELSS SPRSS
;
 
var 
Y 		% output
Ym 		% intermediate output
K 		% capital
L 		% labour supply
I 		% investment
C 		% household consumption
G 		% government spending
Qk 		% price of capital
PI 		% inflation
PIo		% optimal inflation
RHO 	% stochastic discount factor
LAM		% growth rate of stoch. discount factor
Rk 		% return on capital	
R 		% return on safe bond	
N 		% total bank net worth
Ne 		% net worth of existing banks
Nn 		% net worth of entering banks
NU 		% discounted marginal value of capital for bank
ETA 	% discounted marginal value of net worth for bank
PHI		% bank leverage (assets/net worth) 
Z 		% growth rate of bank capital
X 		% growth rate of bank net worth
Pm 		% relative price intermediate goods
U 		% capital utilisation
MU 		% markup
D 		% price dispersion
F 		% recursive optimal inflation term I
H 		% recursive optimal inflation term II
Rn 		% nominal rate on safe bonds
DEL 	% gross depreciation
In  	% net investment
// exog. shock processes
EXA 	% technology (TFP)
EXK 	% capital efficiency
EXG 	% government spending
// for plotting
SPR 	% premium / spread Rk-R
;
 
varexo % shock to ...
e_a    % ... technology (TFP) 
e_k    % ... capital efficiency
e_g    % ... government spending
e_n    % ... net worth
e_i    % ... investment technology
;
 
bet = 0.99;    	% discount factor
sig = 1;       	% intertemp. elasticity of substitution
hab = 0.815;   	% habit parameter
fri = 0.276;   	% Frisch labour elasticity
zet = 7.2;		% utilisation elasticity of depreciation
tet = 0.97155955; % survival probability bankers [why not 0.5^(1/40)?]
alf = 0.33;     % capital share in production
gov = 0.2;      % SS government spending over GDP
acI = 1.728;    % adjustment cost parameter investment
eps = 4.167;    % elasticity of substitution goods
gam = 0.779;    % Calvo probability
kpP = 1.5;      % Taylor rule weight on inflation
kpY = -0.125;   % Taylor rule weight on markup
kap = 10;
// persistence (rho) and standard deviation (sig) of ...
rho_i = 0;   	% ... investment-specific shock
sig_i = 0.01;
rho_k = 0.66;	% ... capital efficiency shock
sig_k = 0.05;
rho_a = 0.95;	% ... TFP shock
sig_a = 0.01;
rho_g = 0.95;  	% ... government spending shock
sig_g = 0.01;
rho_p = 0.66;	% ... intertemp. preference shock
sig_p = 0.072;
sig_n = 0.01;	% ... net worth shock (rho_n=0)

// steady state...
GSS = 0.16975710;	% ... government spending
LSS = 0.33333333; 	% ... labour supply
FSS = 4;          	% ... leverage
ISS = 0.14153927; 	% ... investment
DELSS = 0.025;    	% ... net depreciation
SPRSS = 0.01/4;   	% ... spread Rk-R
% omg = 0.00222778; % ... share of net worth to new banks
% lam = 0.38149499; % ... share of divertable assets
% chi = 3.41080850; % ... labour disutility parameter
% utc = 0.03760101; % ... utilisation cost parameter
% del = 0.02041451; % ... gross depreciation rate

model;
#KSS   = (alf*(eps-1)/eps/(1/bet+SPRSS-1+DELSS))^(1/(1-alf))*LSS;
#YSS   = KSS^alf*LSS^(1-alf);
#utc     = alf*(eps-1)/eps*YSS/KSS;
#del   = DELSS - utc/(1+zet);
#omg   = (1-tet*(SPRSS*FSS + 1/bet))*KSS/FSS/KSS;
#ZSS   = (1-omg*FSS)/tet;
#ETASS = (1-tet)/(1-tet*bet*ZSS);
#NUSS  = (1-tet)*bet*SPRSS/(1-bet*tet*ZSS);
#lam   = NUSS + ETASS/FSS;
#RHOSS = (1-bet*hab)*((1-hab)*(YSS-DELSS*KSS-gov*YSS))^(-sig);
#chi   = RHOSS*(eps-1)/eps*(1-alf)*YSS/LSS^(1+fri);


//Households
//1. Marginal utility of consumption
exp(RHO) =  (exp(C)-hab*exp(C(-1)))^(-sig) - bet*hab*(exp(C(+1))-hab*exp(C))^(-sig);

//2. Euler equation
bet*exp(R)*exp(LAM(+1)) = 1;

//3. Stochastic discount rate
exp(LAM) = exp(RHO)/exp(RHO(-1));

//4. Labor market equilibrium
chi*exp(L)^fri =  exp(RHO)*exp(Pm)*(1-alf)*exp(Y)/exp(L);


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
exp(Qk)*exp(K) = exp(PHI)*exp(N);

//11. Banks' net worth
exp(N) = exp(Ne) + exp(Nn);

//12. Existing banks' net worth accumulation
exp(Ne) = tet*exp(Z)*exp(N(-1))*exp(-sig_n*e_n);

//13. New banks' net worth
exp(Nn) = omg*exp(Qk)*exp(EXK)*exp(K(-1));

//Final goods producer
//14. Return to capital
exp(Rk) =  (exp(Pm)*alf*exp(Ym)/exp(K(-1))+exp(EXK)*(exp(Qk)-exp(DEL)))/exp(Qk(-1));

//15. Production function
exp(Ym) = exp(EXA)*(exp(EXK)*exp(U)*exp(K(-1)))^alf*exp(L)^(1-alf);

//Capital Goods Producer
//16. Optimal investment decision
exp(Qk) = 1 + acI/2*((In+ISS)/(In(-1)+ISS)-1)^2+acI*((In+ISS)/(In(-1)+ISS)-1)*(In+ISS)/(In(-1)+ISS)-bet*exp(LAM(+1))*acI*((In(+1)+ISS)/(In+ISS)-1)*((In(+1)+ISS)/(In+ISS))^2;

//17. Depreciation rate
exp(DEL) = del + utc/(1+zet)*exp(U)^(1+zet);

//18. Optimal capacity utilization rate
exp(Pm)*alf*exp(Ym)/exp(U) = utc*exp(U)^zet*exp(EXK)*exp(K(-1));

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
exp(PIo) = eps/(eps-1)*exp(F)/exp(H)*exp(PI);

//29. Price index
(exp(PI))^(1-eps) = gam*exp(PI(-1))^(kap*(1-eps))+(1-gam)*(exp(PIo))^(1-eps);

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

stoch_simul(order=1, periods=2000, irf=40) EXK R SPR Y C I K L Qk N PI Rn;

