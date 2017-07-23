
// Endogenous variables
var at c ics p ps r rb b bh rp ut v y ;

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

RP=rp*B;





ics - ( XSS*PSI*v^PHI*((y/at)^(1+PHI)) + TETA*BETA*ut*((p(+1)^EPSILON*ics(+1))) ); //same



ics/ps - ( y/c+TETA*BETA*ut*((p(+1)^(EPSILON-1)/ps(+1)*ics(+1))) );



r - RZLB ;//the ZLB

1 - (TETA*p^(EPSILON-1)+(1-TETA)*ps^(1-EPSILON));

v - (TETA*p^EPSILON*v(-1)+(1-TETA)*ps^(-EPSILON));

y - ( c+GYSS*y )  ;

log(ut) = RHO_U*log(ut(-1))+eps_u ;

log(at) = RHO_A*log(at(-1))+eps_a ;

end;

