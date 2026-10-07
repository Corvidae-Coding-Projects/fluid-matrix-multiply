// D3Q19 BGK lattice Boltzmann solver, full volumetric 3D velocity and pressure.
// No hydraulic resistors or prescribed pipe profiles are used by this solver.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
#include <omp.h>
using namespace std;
constexpr int Q=19;
const int c[Q][3]={{0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,-1,0},{1,-1,0},{-1,1,0},{1,0,1},{-1,0,-1},{1,0,-1},{-1,0,1},{0,1,1},{0,-1,-1},{0,1,-1},{0,-1,1}};
const int opp[Q]={0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17};
const double w[Q]={1./3,1./18,1./18,1./18,1./18,1./18,1./18,1./36,1./36,1./36,1./36,1./36,1./36,1./36,1./36,1./36,1./36,1./36,1./36};
struct Node { int x,y,z,port=-1,inside=-1; array<int,Q> nb; };
struct State {double rho,ux,uy,uz;};
State macro(const double* f){State s{0,0,0,0}; for(int q=0;q<Q;q++){s.rho+=f[q];s.ux+=c[q][0]*f[q];s.uy+=c[q][1]*f[q];s.uz+=c[q][2]*f[q];}s.ux/=s.rho;s.uy/=s.rho;s.uz/=s.rho;return s;}
double eq(int q,const State&s){double cu=c[q][0]*s.ux+c[q][1]*s.uy+c[q][2]*s.uz;return w[q]*s.rho*(1+3*cu+4.5*cu*cu-1.5*(s.ux*s.ux+s.uy*s.uy+s.uz*s.uz));}
int main(int argc,char**argv){
 if(argc<4){cerr<<"Usage: solver scale case output-prefix [max-steps] [tolerance] [pressure-Pa] [coarse-field.bin]\nCases: e0 e1 e2 e3 mix signed duct zero\n";return 2;}
 int scale=stoi(argv[1]); string cas=argv[2],prefix=argv[3]; int maxSteps=argc>4?stoi(argv[4]):30000*scale*scale; double tol=argc>5?stod(argv[5]):1e-7,pressure=argc>6?stod(argv[6]):1.;
 if(scale<1||scale>4||pressure<0||pressure>20){cerr<<"Invalid scale or pressure\n";return 2;}
 bool duct=cas=="duct"; const int nx=(duct?64:88)*scale,ny=(duct?12:88)*scale,nz=(duct?12:48)*scale;
 const double tau=.8,nu=(tau-.5)/3,dx=1e-5/scale,dt=nu*dx*dx/1e-6,velUnit=dx/dt,pUnit=1000*velUnit*velUnit,flowUnit=dx*dx*velUnit*1e12; // flow in nL/s
 array<double,4> inputs{0,0,0,0};
 if(cas.size()==2&&cas[0]=='e'&&cas[1]>='0'&&cas[1]<='3') inputs[cas[1]-'0']=1;
 else if(cas=="mix") inputs={.8,.4,.2,.6};
 else if(cas=="signed") inputs={.8,.4,-.2,.6};
 else if(cas=="custom"){ifstream file("inputs.txt");for(double &v:inputs)if(!(file>>v)||!isfinite(v)||abs(v)>1){cerr<<"inputs.txt requires four values in [-1,1]\n";return 2;}}
 else if(duct) inputs[0]=1;
 else if(cas!="zero"){cerr<<"Unknown case\n";return 2;}
 const int centers[4]={16,34,52,70}; int widths[4][4]={{8,6,4,6},{4,8,6,4},{6,4,8,6},{4,6,4,8}};
 ifstream geometry("widths.txt");if(geometry)for(auto &row:widths)for(int &v:row)if(!(geometry>>v)||v<0||v>10||v%2){cerr<<"widths.txt requires 16 even widths in [0,10], in 10 micrometre units\n";return 2;}
 auto fluid=[&](int x,int y,int z){
   if(x<0||x>=nx||y<0||y>=ny||z<0||z>=nz)return false;
   double X=(x+.5)/scale,Y=(y+.5)/scale,Z=(z+.5)/scale;
   if(duct)return Y>=2&&Y<10&&Z>=2&&Z<10;
   for(int i=0;i<4;i++)if(abs(X-centers[i])<6&&Y<80&&abs(Z-36)<6)return true;
   for(int j=0;j<4;j++)if(X>=8&&abs(Y-centers[j])<6&&abs(Z-10)<6)return true;
   for(int j=0;j<4;j++)for(int i=0;i<4;i++)if(abs(X-centers[i])<widths[j][i]/2.&&abs(Y-centers[j])<widths[j][i]/2.&&Z>=10&&Z<36)return true;
   return false;
 };
 vector<int> grid(nx*ny*nz,-1); vector<Node> nodes;
 auto flat=[&](int x,int y,int z){return (z*ny+y)*nx+x;};
 for(int z=0;z<nz;z++)for(int y=0;y<ny;y++)for(int x=0;x<nx;x++)if(fluid(x,y,z)){grid[flat(x,y,z)]=nodes.size();Node n;n.x=x;n.y=y;n.z=z;nodes.push_back(n);}
 auto index=[&](int x,int y,int z){return fluid(x,y,z)?grid[flat(x,y,z)]:-1;};
 vector<int> ports;
 for(int k=0;k<(int)nodes.size();k++){auto &n=nodes[k];for(int q=0;q<Q;q++)n.nb[q]=index(n.x-c[q][0],n.y-c[q][1],n.z-c[q][2]);
   if(duct){if(n.x==0){n.port=0;n.inside=index(1,n.y,n.z);}if(n.x==nx-1){n.port=4;n.inside=index(nx-2,n.y,n.z);}}
   else {if(n.y==0){for(int i=0;i<4;i++)if(abs((n.x+.5)/scale-centers[i])<6)n.port=i;n.inside=index(n.x,1,n.z);}
     if(n.x==nx-1){for(int j=0;j<4;j++)if(abs((n.y+.5)/scale-centers[j])<6)n.port=4+j;n.inside=index(nx-2,n.y,n.z);}}
   if(n.port>=0){if(n.inside<0){cerr<<"Bad boundary\n";return 3;}ports.push_back(k);}
 }
 const int N=nodes.size();vector<double> f(N*Q),next(N*Q),lastU(N*3,0.);
 for(int k=0;k<N;k++)for(int q=0;q<Q;q++)f[k*Q+q]=w[q];
 // Optional coarse-field warm start. It only initializes the full 3D solve;
 // every fine fluid voxel is subsequently evolved by the same LBM equations.
 if(argc>7){
   ifstream seed(argv[7],ios::binary);if(!seed){cerr<<"Cannot open warm start\n";return 2;}
   int cx=duct?64:88,cy=duct?12:88,cz=duct?12:48;vector<array<float,4>> init(cx*cy*cz);vector<bool> valid(cx*cy*cz,false);float a[7];
   while(seed.read((char*)a,sizeof(a))){int x=int(a[0]/10),y=int(a[1]/10),z=int(a[2]/10),idx=(z*cy+y)*cx+x;if(idx>=0&&idx<(int)init.size()){init[idx]={a[3],a[4],a[5],a[6]};valid[idx]=true;}}
   for(int k=0;k<N;k++){auto n=nodes[k];int idx=((n.z/scale)*cy+n.y/scale)*cx+n.x/scale;if(!valid[idx]){cerr<<"Warm start geometry mismatch\n";return 3;}auto a=init[idx];State s{1+3*a[0]/pUnit,a[1]/(1000*velUnit),a[2]/(1000*velUnit),a[3]/(1000*velUnit)};for(int q=0;q<Q;q++)f[k*Q+q]=eq(q,s);}
 }
 // Interior link planes, away from imposed reservoir nodes. These sum exact
 // streamed mass flux, including all diagonal lattice links crossing a plane.
 struct Link{int a,b,q;}; array<vector<Link>,8> links;
 for(int k=0;k<N;k++){auto&n=nodes[k];for(int q=1;q<Q;q++){
   int axis=duct?0:1,plane=4*scale,port=-1;
   if((axis==0?n.x:n.y)==plane&&c[q][axis]==1){if(duct)port=0;else if(n.z>28*scale)for(int i=0;i<4;i++)if(abs((n.x+.5)/scale-centers[i])<6)port=i;}
   if(n.x==nx-4*scale-1&&c[q][0]==1){if(duct)port=4;else if(n.z<18*scale)for(int j=0;j<4;j++)if(abs((n.y+.5)/scale-centers[j])<6)port=4+j;}
   if(port>=0){int dest=index(n.x+c[q][0],n.y+c[q][1],n.z+c[q][2]);if(dest>=0)links[port].push_back({k,dest,q});}
 }}
 auto fluxes=[&](){array<double,8> flux{};for(int p=0;p<8;p++)for(auto l:links[p])flux[p]+=f[l.a*Q+l.q]-f[l.b*Q+opp[l.q]];return flux;};
 ofstream hist(prefix+"-history.csv");hist<<"step,time_s,velocity_relative_change,max_mach,mass_balance_relative,q0_nLs,q1_nLs,q2_nLs,q3_nLs\n";hist<<setprecision(12);
 int checkEvery=200*scale*scale,stable=0,steps=0;double residual=1,maxMach=0,balance=1,rhoMin=1,rhoMax=1;bool converged=false;
 auto start=chrono::steady_clock::now();cerr<<cas<<" scale="<<scale<<" fluid cells="<<N<<" threads="<<omp_get_max_threads()<<"\n";
 for(int t=1;t<=maxSteps;t++){
   #pragma omp parallel for schedule(static)
   for(int k=0;k<N;k++){
     double a[Q];const auto &n=nodes[k];for(int q=0;q<Q;q++)a[q]=n.nb[q]>=0?f[n.nb[q]*Q+q]:f[k*Q+opp[q]];
     State s=macro(a);for(int q=0;q<Q;q++)next[k*Q+q]=a[q]-(a[q]-eq(q,s))/tau;
   }
   #pragma omp parallel for schedule(static)
   for(int pk=0;pk<(int)ports.size();pk++){
     int k=ports[pk],j=nodes[k].inside;State s=macro(&next[j*Q]),target=s;target.rho=1+3*pressure*(nodes[k].port<4?inputs[nodes[k].port]:0)/pUnit;
     for(int q=0;q<Q;q++)next[k*Q+q]=eq(q,target)+next[j*Q+q]-eq(q,s);
   }
   f.swap(next);steps=t;
   if(t%checkEvery==0||t==maxSteps){
     double change=0,norm=0,maxU=0;rhoMin=1;rhoMax=1;
     #pragma omp parallel for reduction(+:change,norm) reduction(max:maxU,rhoMax) reduction(min:rhoMin) schedule(static)
     for(int k=0;k<N;k++){State s=macro(&f[k*Q]);double u[3]={s.ux,s.uy,s.uz};double usq=0;for(int d=0;d<3;d++){change+=(u[d]-lastU[k*3+d])*(u[d]-lastU[k*3+d]);norm+=u[d]*u[d];lastU[k*3+d]=u[d];usq+=u[d]*u[d];}maxU=max(maxU,sqrt(usq));rhoMin=min(rhoMin,s.rho);rhoMax=max(rhoMax,s.rho);}
     residual=maxU<1e-12?0:sqrt(change/max(norm,1e-40));maxMach=maxU*sqrt(3.);auto flux=fluxes();double in=0,out=0,through=0;for(int i=0;i<4;i++){in+=flux[i];out+=flux[4+i];through+=abs(flux[i])+abs(flux[4+i]);}balance=abs(in-out)/max(through/2,1e-30);
     hist<<t<<","<<t*dt<<","<<residual<<","<<maxMach<<","<<balance;for(int j=0;j<4;j++)hist<<","<<flux[4+j]*flowUnit;hist<<"\n";hist.flush();
     if(!isfinite(residual)||maxMach>.2||rhoMin<.8){cerr<<"Unstable run\n";return 4;}
     if(t%(5*checkEvery)==0)cerr<<cas<<" step="<<t<<" residual="<<residual<<" mass="<<balance<<" Ma="<<maxMach<<"\n";
     if(residual<tol)stable++;else stable=0;if(stable>=3){converged=true;break;}
   }
 }
 double seconds=chrono::duration<double>(chrono::steady_clock::now()-start).count();auto flux=fluxes();
 ofstream field(prefix+"-field.bin",ios::binary);
 for(int k=0;k<N;k++){auto n=nodes[k];State s=macro(&f[k*Q]);float a[7]={(float)((n.x+.5)*dx*1e6),(float)((n.y+.5)*dx*1e6),(float)((n.z+.5)*dx*1e6),(float)((s.rho-1)*pUnit/3),(float)(s.ux*velUnit*1000),(float)(s.uy*velUnit*1000),(float)(s.uz*velUnit*1000)};field.write((char*)a,sizeof(a));}
 ofstream js(prefix+".json");js<<setprecision(12)<<"{\"case\":\""<<cas<<"\",\"scale\":"<<scale<<",\"fluid_cells\":"<<N<<",\"shape\":["<<nx<<","<<ny<<","<<nz<<"],\"dx_m\":"<<dx<<",\"dt_s\":"<<dt<<",\"tau\":"<<tau<<",\"viscosity_Pas\":0.001,\"density_kgm3\":1000,\"pressure_Pa\":"<<pressure<<",\"inputs\":[";for(int i=0;i<4;i++)js<<(i?",":"")<<inputs[i];js<<"],\"steps\":"<<steps<<",\"time_s\":"<<steps*dt<<",\"wall_seconds\":"<<seconds<<",\"converged\":"<<(converged?"true":"false")<<",\"velocity_relative_change\":"<<residual<<",\"check_interval_steps\":"<<checkEvery<<",\"max_mach\":"<<maxMach<<",\"rho_min\":"<<rhoMin<<",\"rho_max\":"<<rhoMax<<",\"mass_balance_relative\":"<<balance<<",\"input_mass_equivalent_nLs\":[";for(int i=0;i<4;i++)js<<(i?",":"")<<flux[i]*flowUnit;js<<"],\"output_mass_equivalent_nLs\":[";for(int i=0;i<4;i++)js<<(i?",":"")<<flux[4+i]*flowUnit;js<<"]}\n";
 cerr<<"DONE "<<cas<<" scale="<<scale<<" steps="<<steps<<" converged="<<converged<<" seconds="<<seconds<<"\n";
 return converged?0:5;
}
