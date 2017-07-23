// Endogenous variables
var at c ics p ps r rb b bh rp ut v y ;

//at ut shocks
//c consumption
//p inflation
//ics auxiliary variable X2t
//r nominal interest rate
//ps inflation*
//r nominal interest rate
//v measure of price dispersion 
//y output


varexo eps_a eps_u ;

parameters SIGMA BETA CHI TETA GYSS RP B MU EPSILON PHI TAYLOR_P TAYLOR_Y TAYLOR_R INFSS RZLB RHO_A RHO_U;



//////////////////
model;

#XSS = EPSILON/(EPSILON-1);
#pss = INFSS;
#psss =  ((1-TETA*pss^(EPSILON-1))/(1-TETA))^(1/(1-EPSILON)) ;
#vss = (1-TETA)*psss^(-EPSILON)/(1-TETA*pss^EPSILON);
#rss = pss/BETA;
#icsss = psss/(1-TETA*BETA*pss^(EPSILON-1))/(1-GYSS);
#PSI = icsss*(1-TETA*BETA*pss^EPSILON)/(XSS*vss^PHI) ;
#YSS = (icsss*(1-TETA*BETA*pss^EPSILON)/(XSS*PSI*vss^PHI))^(1/(1+PHI)) ;



c^(-SIGMA)=CHI/bh+BETA*rb*((c(1)^(-SIGMA))/(p(1)));// qe euler equation4.5

rb=r+rp;// decompose the bond return into two conponents // yield decomposition4.6

b=bh; // bond market clearing equilibrium4.7

rp-RP=MU*(b-B); // risk premium channel of bond purchase4.8

RP=rp*B;// during nomal times


ics - ( (EPSILON/(EPSILON-1))*PSI*v^PHI*((y/at)^(1+PHI)) + TETA*BETA*ut*((p(+1)^EPSILON*ics(+1))) ); // A.6

ics/ps - ( y/c+TETA*BETA*ut*((p(+1)^(EPSILON-1)/ps(+1)*ics(+1))) ); //A.7

r - (TAYLOR_R*r(-1)+(1-TAYLOR_R)*(INFSS/BETA*((p/INFSS)^TAYLOR_P)*(y/YSS)^TAYLOR_Y)); //taylor rule A.8

1 - (TETA*p^(EPSILON-1)+(1-TETA)*ps^(1-EPSILON));    //(A.11)

v - (TETA*p^EPSILON*v(-1)+(1-TETA)*ps^(-EPSILON));   //(A.12)

y - ( c+GYSS*y )  ;//(A.13+A.10)


//shock evolve process

log(ut) = RHO_U*log(ut(-1))+eps_u ; //demand shock process

log(at) = RHO_A*log(at(-1))+eps_a ; 

end;



shocks;
  var eps_a; stderr 0.002;
  var eps_u; stderr 0.002;
end;


steady;


stoch_simul(order=1,nocorr,nomoments,irf=0,print);