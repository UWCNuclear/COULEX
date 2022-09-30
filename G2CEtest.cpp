/*Code to determine correlated error using GOSIA2 using method outlined in GOSIA 2012 manual ( T. Czosnyka, D. Cline, C. Y. Wu, A. B. Hayes, \textit{Gosia user manual for simulation and analysis of Coulomb excitation experiments}, 109-114 (2012). )*/
#include <iostream> // IO
#include <sstream>  //string stream
#include <fstream> // ifstream and ofstream
#include <cstdio>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <string> // use strings and some IO with them
#include <cstdlib>
using namespace std;
using std::istringstream;
#define NSTEPS 19
int main (int argc, char **argv )
{
if (argc!=2)
{
cout<< " Usage : " << argv[0] << "input_filename" << endl;
return 1;
}
istringstream filenamein(argv[1]);
string filename, outfilename, targetfilename, targetoutfilename;
string best, errorfile;
filenamein >> filename;
cout<< " Calculate the correlated errors using : " << filename<< endl;
ifstream in,in2,inb;
ofstream out;
ifstream inPO, inTO ; // look in output files for chisq
ofstream outchi; // write-out the chisq values
in.open(filename.c_str(),ios::in); // now look at that Gosia input file
int no;
string line;
in >> no;
cout<< endl;
cout<< " Number of file : " << no << endl;

while (!in.eof())
{
getline(in, line);
if (line.find ("22,3,1")!= string::npos) in>>outfilename;
if (line.find ("26,3,1")!= string::npos) in>>targetfilename;
if (line.find ("12,3,1")!= string::npos) in>>best; 
} 

in.close();
in2.open (targetfilename.c_str(),ios::in);
while (!in2.eof())
{
getline(in2,line);
if (line.find("22,3,1")!=string::npos) in2>>targetoutfilename;
}
in2.close();
errorfile=filename;
size_t sz=errorfile.size();
errorfile.replace(sz-4,4,".error.dat");
system((string("rm ")+errorfile).c_str());
cout << "Previuous chisq_total grid table  " << errorfile.c_str()<< "    was cleaned." <<endl;
cout << "Output file: " << outfilename<< endl;
cout << "Input file target exc. :" << targetfilename<< endl;
cout << "Output file target exc. : " << targetoutfilename<< endl;
cout << "Current chisq_total grid table written to : " << errorfile<< endl;
cout << endl;
double m20start;
double step =0.000526316;
double step2=0.000526316;
double m22start=-0.290239 ;
//int NSTEPS=19;
inb.open(best.c_str(),ios::in);
inb>>m20start;
inb.close();
cout << " Starting transitional matrix element: " << m20start <<endl;
cout << " Step size : " << step <<  endl;
int s1 ,s2;
double m20,m22 ,ll ,ul;
double m20i, m22j;
for (int i=0; i <NSTEPS; i++) // TME calculation
{
m20i=m20start-(NSTEPS-i)*step;
for ( int j =0; j <NSTEPS; j++) // DME calculation (can be modified accordingly)
{
m22j=m22start+step2*(NSTEPS-j);
cout<< j <<'\t'<< m22j<<'\t'<< m20i<<endl;
in.open(filename.c_str(), ios::in); // make adjustments Gosia input file
out.open("temp.inp", ios::out);
while (!in.eof( ))
{
 getline(in,line);
  if (line.find("ME")!=string::npos)
   {
    out << line << endl;
    getline(in,line);
    out << line << endl;
    in >> s1 >> s2 >> m20 >> ll >> ul;
    out<< " " << s1 << " " << s2 << " " << m20i << " " << ul << " " << ul << endl;
in >> s1 >> s2 >> m22 >> ll >> ul;
out << " " << s1 << " " << s2 << " " << m22j << " " << ul << " " << ul << endl;
getline(in ,line );
getline(in , line); 
}
if (line.find("NTAP")!=string::npos)
{
out << " 3 ! NTAP " << endl;
getline(in,line);
}
if (line.find("OP,MINI")!=string::npos)
 {
 getline(in,line);
 break;
 }
 else
 {
  out<<line<<endl;
 }
}
out.close();
in.close( );

system("cat temp.inp INTI.txt > temp2.inp ");
system((string("mv temp2.inp ")+filename).c_str( ));
system((string("gosia2<")+filename+string(" 2>&1 1>/dev/null")).c_str()); // here calculate the corrected yields , keeping the M. E constant
in.open( filename.c_str(), ios::in); // now manipulate Gosia inputfile
out.open("temp.inp",ios::out);
while (!in.eof())
{
getline (in,line);
if (line.find ( "NTAP" )!=string::npos)
{
out << " 4 ! NTAP " << endl;
break;
if (line.find("OP,INTI")!= string::npos)
{
 getline(in, line);
break;
}
else
 {
  out << line << endl;
 }
}
system("cat temp.inp MINI.txt > temp2.inp ");
system((string("mv temp2.inp ")+filename).c_str());
system((string("gosia2 < " )+filename+string(" 2>&1 1>/dev/null")).c_str()); // here doing the minimization, keeping the M.E . constant
out.close();
in.close();
char chiP[20], chiT[20];
inPO.open(outfilename.c_str(), ios::in);
inTO.open(targetoutfilename.c_str(),ios::in);
outchi.open(errorfile.c_str(), ios::app);
while (!inPO.eof())
{
 getline(inPO,line);
 if (line.find("CHISQ=")!=string::npos)
  {
   line.copy(chiP,12,17);
   chiP[12]='\0';
 break;
  }
}
while (!inTO.eof())
{
  getline(inTO, line);
   if (line.find("CHISQ=")!=string::npos)
   {
     chiP[12]='\0';
     break;
     }
    }
while (!inTO.eof())
{
line.copy(chiT,12,17);
 chiT[12]= '\0';
 break;
 }
inPO.close();
inTO.close();
outchi << m22j << '\t' << m20i << '\t' <<'\t'<< chiP << '\t'<< chiT << endl;
outchi.close();
cout << "  calculation status : NSTEPS " << i * NSTEPS+j+1 << " out of " << NSTEPS*NSTEPS << " was done . " << endl;
} // end o f m22 l o o p
outchi.open( errorfile.c_str(), ios::app);
outchi.close();
cout << " After step " << (i +1)*NSTEPS << " the error data looks the following: " << endl;
system((string( "cat  ")+errorfile ).c_str());
} // end o f m20 l o o p
cout << " The correlated error file  " << errorfile.c_str( ) << " was sucessfully calculated . " << endl;
}
}
