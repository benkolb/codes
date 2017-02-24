// Dynare model file to calculate the GK model
// Original code created by Peter Karadi, July 2010
// Modified by Benedikt Kolb, Feb. 2017
// changes include:
// - correction of equation 6, "Value of banks' net wealth"
// - no more capital utilisation, net investment
// - substitute for some variables (MU, LAMBDA, Ym, X, Ne, Nn)

parameters  % explanations below
bet sig hab fri zet tet alf gov acI eps gam kpP kpY kap del
// exogenous-process parameters
rho_a rho_i rho_k rho_g rho_p 
sig_a sig_i sig_k sig_g sig_p sig_n 
// steady states 
LSS FSS SPRSS
;
 
var 
Y 		% output
K 		% capital
L 		% labour supply
I 		% investment
C 		% household consumption
G 		% government spending
Q 		% price of capital
PIP 	% inflation
PIS		% optimal inflation
RHO 	% stochastic discount factor
RK 		% return on capital	
R 		% return on safe bond	
N 		% total bank net worth
NU 		% discounted marginal value of capital for bank
ETA 	% discounted marginal value of net worth for bank
PHI		% bank leverage (assets/net worth) 
Z 		% growth rate of bank capital
PM 		% relative price intermediate goods
D 		% price dispersion
F 		% recursive optimal inflation term I
H 		% recursive optimal inflation term II
RN 		% nominal rate on safe bonds
// exog. shock processes
EXA 	% technology (TFP)
EXK 	% capital efficiency
EXG 	% government spending
// for plotting
SPR 	% premium / spread RK-R
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
tet = 0.97155955; % survival probability bankers, see GK11
alf = 0.33;     % capital share in production
gov = 0.2;      % SS government spending over GDP
acI = 1.728;    % adjustment cost parameter investment
eps = 4.167;    % elasticity of substitution goods
gam = 0.779;    % Calvo probability
kpP = 1.5;      % Taylor rule weight on inflation
kpY = -0.125;   % Taylor rule weight on markup
kap = 10;
del = 0.025;  	% depreciation rate

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
LSS = 1/3; 	% ... labour supply
FSS = 4;          	% ... leverage
SPRSS = 0.01/4;   	% ... spread RK-R

model;
#KSS   = (alf*(eps-1)/eps/(1/bet+SPRSS-1+del))^(1/(1-alf))*LSS;
#YSS   = KSS^alf*LSS^(1-alf);
#omg   = (1-tet*(SPRSS*FSS + 1/bet))*KSS/FSS/KSS;
#ZSS   = (1-omg*FSS)/tet;
#ETASS = (1-tet)/(1-tet*bet*ZSS);
#NUSS  = (1-tet)*bet*SPRSS/(1-bet*tet*ZSS);
#lam   = NUSS + ETASS/FSS;
#RHOSS = (1-bet*hab)*((1-hab)*(YSS*(1-gov)-del*KSS))^(-sig);
#chi   = RHOSS*(eps-1)/eps*(1-alf)*YSS/LSS^(1+fri);


//HOUSEHOLDS
//1. Marginal utility of consumption
exp(RHO) =  (exp(C)-hab*exp(C(-1)))^(-sig) - bet*hab*(exp(C(+1))-hab*exp(C))^(-sig);

//2. Euler equation
bet*exp(R)*exp(RHO(+1))/exp(RHO) = 1;

//3. Labor market equilibrium
chi*exp(L)^fri =  exp(RHO)*exp(PM)*(1-alf)*exp(Y)/exp(L);


//BANKS
//4. Value of banks' capital
exp(NU) = (1-tet)*bet*exp(RHO(+1))/exp(RHO)*(exp(RK(+1))-exp(R))+bet*exp(RHO(+1))/exp(RHO)*tet*exp(PHI(+1))/exp(PHI)*exp(Z(+1))*exp(NU(+1));

//5. Value of banks' net wealth
exp(ETA) = bet*exp(RHO(+1))/exp(RHO)*( (1-tet)*exp(R) + tet*exp(Z(+1))*exp(ETA(+1)) );
//NOTE: The corresponding equation in Peter Karadi's code is:
// exp(ETA) =  (1-tet) + bet*exp(RHO(+1))/exp(RHO)*tet*exp(Z(+1))*exp(ETA(+1));

//6. Optimal leverage
exp(PHI) = exp(ETA)/(lam-exp(NU));

//7. Growth rate of banks' capital
exp(Z) = (exp(RK)-exp(R(-1)))*exp(PHI(-1)) + exp(R(-1));

//8. Aggregate capital
exp(Q)*exp(K) = exp(PHI)*exp(N);

//9. Banks' net worth
exp(N) = tet*exp(Z)*exp(N(-1))*exp(-sig_n*e_n) + omg*exp(Q)*exp(EXK)*exp(K(-1));


//GOOD AND CAPITAL PRODUCTION
//10. Return to capital
exp(RK) =  (exp(PM)*alf*exp(Y)*exp(D)/exp(K(-1))+exp(EXK)*(1-del)*exp(Q))/exp(Q(-1)); % different to GK11!

//11. Production function
exp(Y) = exp(EXA)*(exp(EXK)*exp(K(-1)))^alf*exp(L)^(1-alf)/exp(D);

//12. Optimal investment decision
exp(Q) = 1 + acI/2*(exp(I)/exp(I(-1))-1)^2 + acI*(exp(I)/exp(I(-1))-1)*exp(I)/exp(I(-1)) - bet*exp(RHO(+1))/exp(RHO)*acI*(exp(I(+1))/exp(I)-1)*(exp(I(+1))/exp(I))^2;

//13. Capital accumulation equation
exp(K) = (1-del)*exp(EXK)*exp(K(-1)) + exp(I); 


//RETAILER
//14. Price dispersion
exp(D) = gam*exp(D(-1))*exp(PIP(-1))^(-kap*eps)*exp(PIP)^eps+(1-gam)*((1-gam*exp(PIP(-1))^(kap*(1-gam))*exp(PIP)^(gam-1))/(1-gam))^(-eps/(1-gam));

//15. Optimal price choice I
exp(F) =  exp(Y)*exp(PM)+bet*gam*exp(RHO(+1))/exp(RHO)*exp(PIP(+1))^eps*(exp(PIP))^(-eps*kap)*exp(F(+1));

//16. Optimal price choice II
exp(H) = exp(Y)+bet*gam*exp(RHO(+1))/exp(RHO)*exp(PIP(+1))^(eps-1)*exp(PIP)^(kap*(1-eps))*exp(H(+1));

//17. Optimal price choice III
exp(PIS) = eps/(eps-1)*exp(F)/exp(H)*exp(PIP);

//18. Price index
(exp(PIP))^(1-eps) = gam*exp(PIP(-1))^(kap*(1-eps))+(1-gam)*(exp(PIS))^(1-eps);


//EQUILIBRIUM
//19. Government consumption
exp(G) = gov*YSS*exp(EXG);

//20. Aggregate resource constraint
exp(Y) = exp(C)+exp(G)+exp(I)+acI/2*(exp(I)/exp(I(-1))-1)^2*exp(I);

//21. Fisher equation
exp(RN) = exp(R)*exp(PIP(+1));

//22. Interest rate rule
exp(RN)      =   exp(RN(-1))^rho_i*((1/bet)*exp(PIP)^kpP*((eps-1)/eps/exp(PM))^(kpY))^(1-rho_i)*exp(sig_i*e_i);


//ADDITIONAL VARIABLES
//23. Premium
exp(SPR) = exp(RK(+1))/exp(R);


//SHOCKS
//24. TFP shock
EXA = rho_a*EXA(-1) - sig_a*e_a;

//25. Capital quality shock
EXK = rho_k*EXK(-1) - sig_k*e_k;

//26. Government consumption shock
EXG = rho_g*EXG(-1) - sig_g*e_g;

end;


steady_state_model; % all not mentioned are zero
L   = log(LSS);
PHI = log(FSS);
R   = log(1/bet);
PM  = log((eps-1)/eps);
RK  = log(exp(R)+SPRSS);
RN  = R;
K   = log((alf*exp(PM)/(exp(RK)-1+del))^(1/(1-alf))*LSS); 
Y   = log(exp(K)^alf*exp(L)^(1-alf));
G   = log(gov*exp(Y));
I   = log(del*exp(K));
C   = log(exp(Y)-exp(I)-exp(G));
RHO = log((1-bet*hab)*((1-hab)*exp(C))^(-sig));
N   = log(exp(K)/exp(PHI));
Z   = log(SPRSS*FSS + 1/bet);
NU  = log(((1-tet)*bet*SPRSS)/(1-bet*tet*exp(Z)));
ETA = log((1-tet)/(1-tet*bet*exp(Z)));
F   = log(exp(Y)*exp(PM)/(1-bet*gam));
H   = log(exp(Y)/(1-bet*gam));
SPR = log(exp(RK)/exp(R));
end;


// SIMULATION.

model_diagnostics;
resid;
check;
steady;

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

stoch_simul(order=1, periods=2000, irf=40) EXK SPR Q Y PIP RN I K N;

