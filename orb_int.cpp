//Kinematics and orbital integral Calculator
#include <iostream>
#include <cmath> 
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include "k2.h"
#include "orb_info_pr325.h"
#define nn 24
using namespace std;
double nrad;
double nu_r(double mass) {
	nrad = 1.25 * pow(mass,0.3333);
return (nrad);
}
const double hbarc = 1.973*pow(10,-13);
const double n_a = 6.23 * pow(10,23);
const double pi = 3.141592654 ;
const double fm = pow(10,-15);
const double e =  1.6*pow(10,-19);
const double barn = pow(10,-28);
const double nA = pow(10,-9);
const double hour = 3600;
const double alpha = 1.44;
double ecoulb( int zproj, int ztar, int mproj, int mtar)
{
 double v1;
 v1 = 1.44*(zproj*ztar)/(nu_r(mproj) + nu_r(mtar));
 return (v1);
 }
double esafe( int zproj, int ztar, int mproj, int mtar, double sd)
 {
  double u1;
  u1 = 1.44*((mproj + mtar)/(mtar))*((zproj*ztar)/(nu_r(mproj) + nu_r(mtar) + sd));
  return (u1);
 }

double doca(double x1, double x2,double x3,double x4, double x5)
 {
 double d1; 
  d1 = ((alpha*x1*x2)/(x3))*(1 + (x4/x5));
  return (d1);
}
double ang_degg( double x1) 
 {
  double x2; 
  x2 = x1*(180/pi);
  return (x2);
 }
double oparam(double xx)
 {
  double op1;
  op1 = (pi-xx)/2;
  return op1;
 }
double adia_param(double som1, double som2)
{
 double xxx;
 xxx = som2-som1;
 return xxx;
 } 

int main()
{ 
int lam,mew;
double  u, m_t, m_p, z_t, z_p, beamenergy_lab, v_lab, mu, a, b, beamenergy_cm,d_ft,d,m_rat,lamdabar;	
double  data_array[24][2]; 
double adia_f[24];
double  r_i,r_o,theta_labi,theta_labo,theta_cmi,theta_cmo,theta_mid;
double sigma_ruti,sigma_ruto,sigma_rut,dx_t,rho_a,eta,r_t,r_b,de_srim,elas_p,adx_t;
double s_thetacm;
//double a2;
double Theta_labav[24];
double aptdx[nn];
double de_loss[nn];
double ebeam[nn];
double beta_r[nn];
double orbital_int[nn];
double orbital_int_a[nn];
double cond1[nn];
double cond2[nn];
double oi[1000];
double oi1[1000];
double oi2[1000];
double oi3[1000];
double oi4[1000];
double oi5[1000];
double oi6[1000];
double oi7[1000];
double oi8[1000];
double oi9[1000];
double gtest,gtest2;
double safe_eng;
/*double n23[16384];
double e23[16384];
double efr[16384];*/
/*int i1=50;
int i2 = 70;*/
fstream out("OI_TEST_e112");
u = 931.494;			// amu in MeV/c^2
const double alpha = 1.44;
m_t = 208;
m_p = 18;
z_t = 82;
z_p = 8;
beamenergy_lab	= 53;
d_ft = 25;	//1.44 MeV.fm
/*cout << "Enter target mass in amu\n";
cin >> m_t;
cout << "Enter projectile mass in amu\n";
cin >> m_p;
cout << "Enter Z value of target\n";
cin >> z_t;
cout << "Enter Z value of beam\n";
cin >> z_p;
cout << "Enter beam Energy in MeV\n";
cin >> beamenergy_lab;
cout << "Enter distance from target\n";
cin >> d_ft ;*/
/*cout << "Enter beam current\n";
cin >> i_b;
cout << "Enter efficiency\n";
cin >> eff;*/
/*cout << "Enter target thickness\n";
cin >> dx_t;
cout << "Enter srim de/dx @ the beam energy entered\n";
 cin >> de_srim;
 cout << "Enter elastic peak energy in MeV\n";
 cin >> elas_p;*/
/*cout << " Enter lambda\n";
 cin >> lam;
cout << " Enter mu\n";
 cin >> mew;*/
v_lab = sqrt((2 * beamenergy_lab) / (u * m_p));
mu = u * m_t * m_p / ( m_t + m_p);
a=doca(z_p,z_t,alpha,mu,v_lab,v_lab);
//a = (z_p * z_t * alpha) / (mu * v_lab * v_lab ); 
b = 2*a;
m_rat = m_p/m_t;
beamenergy_cm = 0.5 * mu * v_lab * v_lab;
d = pow(10,3)*d_ft;
lamdabar = hbarc/sqrt(2*beamenergy_lab*m_p*u);
eta = a*fm/lamdabar;
r_t = nu_r(m_t);
r_b =nu_r(m_p);
rho_a = (dx_t * n_a)/m_t;
adx_t = tt(dx_t,rho_a);
gtest=0.349065850;
safe_eng = esafe(z_p,z_t,m_p,m_t,6.5);
cout << a << '\t' << eta << '\t'<< safe_eng<<              endl;

//i_ips = (etai_b*nA)/e;
        ifstream in_file("RUTDATA.dat", ios::binary);
        //Check if the file is open
        if(!in_file.is_open()){
                cout << "File not opened..." << endl;
                return 1;
        }
        //read the file (two columns) to the array
        for(int j=0; !in_file.eof(); j++){
                in_file >> data_array[j][0];
                in_file >> data_array[j][1];
                
        }

        //Display the array
 for(int j=0; j<nn; j++){
        r_i = data_array[j][0];
        r_o = data_array[j][1];
	theta_labi = pi-atan(r_i/d);	
	theta_labo = pi-atan(r_o/d);
        Theta_labav[j] = (theta_labi + theta_labo)/2;
	theta_cmi = theta_labi + asin(m_rat*sin(theta_labi));
	theta_cmo = theta_labo + asin(m_rat*sin(theta_labo));
        theta_mid = (theta_cmi + theta_cmo)/2;
	sigma_ruti = (2*pow(a,2)*pi)/(cos(theta_cmi)-1);
	sigma_ruto = (2*pow(a,2)*pi)/(cos(theta_cmo)-1);
	sigma_rut = (sigma_ruto-sigma_ruti)*(pow(fm,2)/barn);
	//n_c = (rho_a * eff *sigma_rut * i_ips *barn);
        aptdx[j]= adx(adx_t,Theta_labav[j]);
       de_loss[j]=aptdx[j]*de_srim;
       ebeam[j]=elabf(m_p,m_t,beamenergy_lab,theta_mid,elas_p) -de_loss[j];
       beta_r[j]=beta_b(ebeam[j],m_p);
       //adia_f[j]=adia_param(eta,doca(z_p,z_t,alpha,mu,beta_r[j],beta_r[j])*fm/lamdabar);
       adia_f[j]=doca(z_p,z_t,alpha,mu,v_lab,beta_r[j])*fm/lamdabar;
        cout <<  adia_f[j] <<endl;
       orbital_int[j]=orb_int(lam,mew,0.0,adia_f[j])
       s_thetacm=0.72*((z_p*z_t)/(beamenergy_lab))*(1+m_rat)*(1+(1/sin(0.5*theta_mid)))-(r_t+r_b);
	//cout << (theta_labi)*(180/pi)<< '\t' << (theta_labo)*(180/pi) <<'\t'<< s_thetacm << endl;
         //cout << ang_degg(theta_mid)<< '\t'<< orbital_int[j]<<endl;
          cond1[j]=cond_0(theta_mid);
          cond2[j]=cond_1(adia_f[j]);
          orbital_int_a[j]=orb_int_asym(lam,mew,0.0,adia_f[j]);
         //cout<<ang_degg(theta_mid)<< '\t'<< cond1[j]<<'\t'<<cond2[j]<<'\t'<< orbital_int[j]<<'\t'<<whit_asym(0.5*(omeg(mew,lam)+1), 1-lam,2*adia_f[j]*CSC(theta_mid/2))<<endl;
           //cout<<adia_f[j]<< '\t'<< orbital_int[j]<<endl;
        //cout  << (theta_mid)*(180/(3.14592654))<<'\t' << s_thetacm<<endl;
        }
        
  /*for(int jj = 0;jj<1000;jj++){
     gtest2 = 0.0018*jj;
     oi[jj]=orb_int(lam,mew,3.141592654-8*gtest,gtest2);
     oi1[jj]=orb_int(2,0,3.141592654-7*gtest,gtest2);
     oi2[jj]=orb_int(2,-2,3.141592654-7*gtest,gtest2);
     oi3[jj]=orb_int(2,2,3.141592654-7*gtest,gtest2);
     oi4[jj]=orb_int(3,0,3.141592654-7*gtest,gtest2);
     oi5[jj]=orb_int(3,2,3.141592654-7*gtest,gtest2);
     oi6[jj]=orb_int(3,-2,3.141592654-7*gtest,gtest2);
     oi7[jj]=2*oi1[jj]*oi4[jj]+3*oi3[jj]*oi5[jj]+3*oi2[jj]*oi6[jj];
     oi8[jj]=fabs(oi1[jj]*oi1[jj])+1.5*fabs(oi3[jj]*oi3[jj])+1.5*fabs(oi2[jj]*oi2[jj]);
     oi9[jj]=frac(oi7[jj],oi8[jj]);
   out << gtest2 <<'\t'<< k_func(lam,mew,pie-5*gtest,gtest2)<<endl;
   }*/
//cout<<sph_harm_900(lam,mew,) <<endl;

//gtest=tgamma(1-lam)/tgamma(0.5*(omeg(mew,lam)+1)+lam);
//gtest=tgamma(1-lam);
//cout<<gtest<<endl;
/*gtest=tgamma(0.5);
cout<< gtest<<endl;*/
/*-cout<<"valuearray Theta_labav[1:24]"<<endl;
cout<<Theta_labav[0]<<" "<<Theta_labav[1]<<" "<<Theta_labav[2]<<" "<<Theta_labav[3]<<" "<<Theta_labav[4]<<" "<<Theta_labav[5]<<endl;
cout<<Theta_labav[6]<<" "<<Theta_labav[7]<<" "<<Theta_labav[8]<<" "<<Theta_labav[9]<<" "<<Theta_labav[10]<<" "<<Theta_labav[11]<<endl;
cout<<Theta_labav[12]<<" "<<Theta_labav[13]<<" "<< Theta_labav[14]<<" "<<Theta_labav[15]<<" "<<Theta_labav[16]<<" "<<Theta_labav[17]<<endl;
cout<<Theta_labav[18]<<" "<<Theta_labav[19]<<" "<<Theta_labav[20]<<" "<<Theta_labav[21]<<" "<<Theta_labav[22]<<" "<< Theta_labav[23]<<endl;*/        
/*srand((unsigned)time(0)); 
for (int h =0; h < 16384;h++)
{
 e23[h] = i1 + int(1.0*(i2-i1+1)*rand()/(RAND_MAX+1.0));
 efr[h] = doca(z_p,z_t,e23[h],m_p,m_t);
 n23[h]= (rho_a * eff * i_ips *barn)*((2*pow(efr[h],2)*pi)/(cos(114*pi/180)-1));
 out1<< e23[h] << " " << -1*n23[h]<<endl;
 }*/
/*cout << "v_lab is\n" << v_lab << "c" << endl;
cout << "Reduced mass of two body system is\n" << mu << "MeV/c^2" << endl;
cout << "Distance of closest approach is\n" << a << "fm" << endl;
cout << "impact parameter is\n" << b << "fm " << endl;
cout << "sommerfled parameter is\n" << eta << endl;
cout << "Beam Energy in cm frame\n" << beamenergy_cm << "MeV" << endl;*/
return 0;
}
