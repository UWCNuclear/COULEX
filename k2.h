
const double pie = 3.141592654;
//const double u = 931.494;
double frac(double h1,double h2)
 {
 double h3;
  h3 = h1/h2;
 return (h3);
 }
double CSC(double csc0)
{
 double CSC1;
CSC1=frac(1,sin(csc0));
return CSC1;
}
double si_detconfig(int u, double ang)
{
 double u1;
 if (u ==1)
  { 
   u1 = ang;
  } 
 else if (u ==2)
 {
 u1 = pie - ang;
 }
return(u1);
}
 
double si_swp(int u, double u2)
 {
 double u1;
 if (u==1)
  {
   u1 = u2;
 }
else if (u==2)
 {
 u1 = -u2;
 }
return(u1);
}
double d2tar(double u)
 {
 double u1;
 u1 = pow(10,3)*u;
 return (u1);
}
double red_en( double m1, double m2,double  eb,double  e_ex){
double eps;
eps = eb -( 1+ frac(m1,m2))*e_ex;
return (eps);
}
double lab_ang( double r, double dtr){
double la;
la=atan(r/dtr);
return (la);
}
double ang_degrees(double angle) {
           double deg;
           deg = (180/pie)*angle;
           return (deg);
}
double cm_ang( double m1,double  m2,double  eb,double lag, double e_ex){
double c1,c2,c3;
c1=red_en(m1,m2,eb,e_ex);
c2=frac(m1,m2)*sqrt(frac(eb,c1));
c3 = lag + asin(c2*lag);
return (c3); 
}

double rec_ang(double cma){
double s1;
s1=pie- cma;
return(s1);
}
double tau( double m1, double m2,double  eb,double  e_ex)
 {
 double g1,g2,g3;
 g1= frac(m1,m2);
 g2 = sqrt(1-frac(e_ex,eb)*(1+g1));
 g3 = frac(g1,g2);
 return(g3);
}

double elabf( double m1,double  m2,double  eb,double lag, double e_ex)
 {
 double c1,c2,c3,c4;
 c1 = tau(m1,m2,eb,e_ex);
 c2 = (1+pow(c1,2)+2*c1*cos(lag))*red_en(m1,m2,eb,e_ex);
 c3 = frac(m2,m1+m2);
 c4 = pow(c3,2)*c2;
 return (c4);
 }

double tt(double p,double p_m){
        double dx_tar;
        dx_tar =(p/p_m);
        return (dx_tar);
}

double adx(double dxxx,double thtl){
       double ap_th;
       ap_th = 0.5*dxxx + (dxxx/(2*cos(thtl)));
       return (ap_th);
}
double beta_b(double eee , double mmm){
	double beta_g=sqrt(frac(2*eee,mmm));
	return(beta_g);
}
double omeg(int c, int cc){
//c= lambda, cc= mu
 double ccc;
 ccc= c-cc;
 return(ccc);
}
double omeg_conj(int d, int dd){
//d=lambda,dd=mu
double ddd;
 ddd=-d-dd;
return(ddd);
}
double orb_int_1(int f,int ff){
        double fff;
          fff=frac(2*pie,tgamma(frac(omeg(f,ff)+1,2)));
         return(fff);
}
double orb_int_2(double g,double gg){
 double ggg;
 ggg=exp(-gg*(CSC(frac(g,2))+frac(pie,2)));
 //ggg=exp(-gg*CSC(frac(gg,2));
  return(ggg);
}
double orb_int_3(double h, int hh,int hhh){
double hhhh;
hhhh=pow(h,frac(omeg(hh,hhh)-1,2));
return (hhhh);
}
double orb_int_4(double i, int ii, int iii){
  double iiii;
  double iiiii;
  iiii=frac(omeg_conj(ii,iii)-1,2);
  iiiii=pow(2*CSC(frac(i,2)),iiii);
 return(iiiii);
}
double orb_int(int k, int kk, double kkk, double kkkk)
//lambda,mu,cm_ang,adia
{
 double k5;
 k5 = orb_int_1(k,kk)*orb_int_2(kkk,kkkk)*orb_int_3(kkkk,k,kk)*orb_int_4(kkk,k,kk);
 return (k5);
 }
double cond_0(double l){
double ll;
ll = sqrt(pow(CSC(l),2)-1);
return (ll);
}
double cond_1(double lll){
double l4;
l4=frac(1,cbrt(lll));
return(l4);
}
double orb_int_5(int m, int mm){
double mmm;
mmm= pow(-1,omeg(m,mm)/2);
return(mmm);
}
double orb_int_6(double n){
double nn;
nn = frac(1,CSC(n/2));
return(nn);
}
double orb_int_7(int o, double o1, double o2){
//lambda,cm_ang,adia
double o3,o4;
o3=frac(o2,2*CSC(o1/2));
o4=pow(o3,(o-1)/2);
return(o4);
}

double orb_int_8(int p,int p1){
//lambda,mu
double p2;
p2=tgamma(frac(omeg(p1,p)+1,2));
return(p2);
}
double orb_int_9(double q){
// q=adia;
double q1;
q1 = exp(-1*(pie*q)/2);
return q1;
}
double hyp_geom(double r,int r1, double r2)
{
//r=0.5*(omeg(mu,lambda)+1),r1=1-lambda,r2=2*adia*CSC(cm_ang/2); 
 double r3,r4,r5; 
 r3 = frac(exp(r2)*pow(r2,r-r1),tgamma(r));
 r4 = frac(pow(-1*r2,-1*r),tgamma(r1-r));
 r5 = tgamma(r1)*(r3+r4);
 return (r5);
}
double whit_asym(double s, int s1, double s2){
 //s=0.5*(omeg(mu,lambda)+1),s1=1-lambda,s2=2*adia*CSC(cm_ang/2); 
 double whit_fac1,whit_fac2,whit_func;
 whit_fac1=frac(tgamma(1-s1),tgamma(s+(s1)));
 whit_fac2=frac(tgamma(s1-1),tgamma(s))*pow(s2,1-s1);
 whit_func=whit_fac1*hyp_geom(s,s1,s2)+whit_fac2*hyp_geom(s+1-s1,2-s1,s2);
 return(whit_func);
 }
 double whit_fac(int t,double t1, double t2){
  //t=lambda,t1=cm_ang,t2=adia;
  double wf1,wf2,wf3;
   wf1 = 2*t2*CSC(t1/2);
   wf2 = exp(-1*(wf1/2));
   wf3 = wf2*pow(wf1,0.5*(1-t));
  return(wf3);
 }
double orb_int_asym(int u, int u1, double u2, double u3)
{
//u=lambda,u1=mu,u2=cm_ang,u3=adia;
 double obs;
 obs = orb_int_5(u,u1)*orb_int_6(u)*orb_int_7(u,u2,u3)*orb_int_8(u,u1)*orb_int_9(u3)*whit_fac(u,u1,u2)*whit_asym(0.5*(omeg(u,u1)+1),1-u,2*u3*CSC(u2/2));
  return(obs);
}
int fact(int n) { 
   if ((n==0)||(n==1))
      return 1; 
   else
      return n*fact(n-1);
}
int dfact(int n)
 
{
       int i;double res=1.0;
 
       for(i=n;i>=1;i-=2)
       {
       res *=i;
       }
 
       return res;
}
double sph_harm_900(int o, int oo){
 double ooo;
 ooo=sqrt(frac(2*o+1,4*pie))*frac(sqrt(fact(o-oo)*fact(o+oo)),dfact(o-oo)*(o+oo))*pow(-1,frac(oo+o,2));
 return (ooo);
}
double k_func(int p, int pp, double ppp, double pppp){
double p5;
p5 = sqrt(pie)*frac(dfact(2*p-1),fact(p-1))*sph_harm_900(p,pp)*orb_int(p,pp,ppp,pppp);
return(p5);
}
