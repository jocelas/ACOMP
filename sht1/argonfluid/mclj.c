/**********************************************************************
 *
 * File: mclj.c
 *
 * NVT-Monte Carlo of Lennard-Jonesium
 *
 * 12-May-2007 (MN)
 * 19-Apr-2012
 * 25-Feb-2016 (added pdb output + input file dump, MS)
 * 10-Mar-2020 (added extented xyz output, CL)
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#define TRUE 1
#define FALSE 0
#define min(a,b) ((a)<(b)?(a):(b))

void usage(char** argv, int help);

unsigned short seed[3];
int n,nt,ntjob,ntprint,ntskip;
int nacc,ndr;
double disp,dr,rho,t;
double alnmax,c,c2,c2m1,drm1,r2max,rc,rc2,su0,sw0;
double *x,*y,*z;
double *un,*umatb;
double **umat;
long long *ag;
double accr,au,au2,aw;


/**********************************************************************/

void getcf(char **argv) {

/**********************************************************************/

  FILE *fpi;
  int i;

  alnmax=log(0.1*DBL_MAX);

  /* Read startup/checkpoint file */

  if((fpi=fopen("mclj_in.dat","r"))==NULL) usage(argv,0);

  fread(&n,sizeof(int),1,fpi);
  fread(&nt,sizeof(int),1,fpi);
  fread(&ntjob,sizeof(int),1,fpi);
  fread(&ntprint,sizeof(int),1,fpi);
  fread(&ntskip,sizeof(int),1,fpi);
  fread(seed,sizeof(unsigned short),3,fpi);

  fread(&disp,sizeof(double),1,fpi);
  fread(&dr,sizeof(double),1,fpi);
  fread(&rho,sizeof(double),1,fpi);
  fread(&t,sizeof(double),1,fpi);

  /* Box parameters */

  c=cbrt(n/rho);
  c2=0.5*c;
  c2m1=1.0/c2;
  r2max=c2*c2;
  drm1=1.0/dr;
  ndr=c2*drm1;

  /* Allocate arrays */

  ag=(long long*)malloc(ndr*sizeof(long long));
  umat=(double**)malloc(n*sizeof(double*));
  umatb=(double*)malloc(n*n*sizeof(double));
  un=(double*)malloc(n*sizeof(double));
  x=(double*)malloc(n*sizeof(double));
  y=(double*)malloc(n*sizeof(double));
  z=(double*)malloc(n*sizeof(double));

  /* Positions & accumulated averages */

  fread(x,sizeof(double),n,fpi);
  fread(y,sizeof(double),n,fpi);
  fread(z,sizeof(double),n,fpi);

  if(nt>0) {

    fread(&accr,sizeof(double),1,fpi);
    fread(&au,sizeof(double),1,fpi);
    fread(&au2,sizeof(double),1,fpi);
    fread(&aw,sizeof(double),1,fpi);

    fread(ag,sizeof(long long),ndr,fpi);

  } else {

    accr=0.0;
    au=0.0;
    au2=0.0;
    aw=0.0;

    for(i=0;i<=ndr-1;i++) {
      ag[i]=(long long)0;
    }

  }

  fclose(fpi);

  /* Cutoff & tail corrections */

  rc=c2;
  rc2=rc*rc;
  su0=2.0*M_PI*rho*n*(4.0/(9.0*pow(rc,9))-4.0/(3.0*pow(rc,3)));
  sw0=2.0*M_PI*rho*n*(24.0/(3.0*pow(rc,3))-48.0/(9.0*pow(rc,9)));

  /* Initialize RNG & acceptance counter */

  seed48(seed);
  nacc=0;

  return;

}


/**********************************************************************/

int dumpcf(char ** argv) {

/**********************************************************************/

  FILE *fpi;
  int i;

  /* Read startup/checkpoint file */

  getcf(argv);

  printf("n = %d\n",n);
  printf("nt = %d\n",nt);
  printf("ntjob = %d\n",ntjob);
  printf("ntprint = %d\n",ntprint);
  printf("ntskip = %d\n",ntskip);
  printf("seed = %ld\n",(size_t)seed);

  printf("disp = %f\n",disp);
  printf("dr   = %f\n",dr  );
  printf("rho  = %f\n",rho );
  printf("t    = %f\n",t   );

  for(i=0;i<n;i++) printf("%f %f %f\n",x[i],y[i],z[i]);

  return 0;
}


/**********************************************************************/

void uinit() {

/**********************************************************************/

  int i,j;
  double dx,dy,dz,r2,rm6;

  /* Initialize pair interaction matrix */

  for(i=0;i<=n-1;i++) {
    umat[i]=umatb+i*n;
  }

  for(i=0;i<=n-2;i++) {
    umat[i][i]=0.0;
    for(j=i;j<=n-1;j++) {
      dx=x[j]-x[i];
      dx=dx-(int)(dx*c2m1)*c;
      dy=y[j]-y[i];
      dy=dy-(int)(dy*c2m1)*c;
      dz=z[j]-z[i];
      dz=dz-(int)(dz*c2m1)*c;
      r2=dx*dx+dy*dy+dz*dz;
      if(r2<rc2) {
        rm6=1.0/(r2*r2*r2);
        umat[i][j]=(4.0*rm6-4.0)*rm6;
      } else {
        umat[i][j]=0.0;
      }
      umat[j][i]=umat[i][j];
    }
  }
  umat[n-1][n-1]=0.0;

  return;

}

/**********************************************************************/

void move() {

/**********************************************************************/

  int accept;
  int i,j;
  double dx,dy,dz,r2,rm6,su,xin,yin,zin;

  /* Random particle */

  i=min((int)(drand48()*n),n-1);

  /* Trial move */

  xin=x[i]+disp*(drand48()-0.5);
  xin=xin-(int)(xin*c2m1)*c;
  yin=y[i]+disp*(drand48()-0.5);
  yin=yin-(int)(yin*c2m1)*c;
  zin=z[i]+disp*(drand48()-0.5);
  zin=zin-(int)(zin*c2m1)*c;

  /* Energy difference */

  su=0.0;
  for(j=0;j<=n-1;j++) {
    if(j==i) {
      un[j]=0.0;
    } else {
      dx=x[j]-xin;
      dx=dx-(int)(dx*c2m1)*c;
      dy=y[j]-yin;
      dy=dy-(int)(dy*c2m1)*c;
      dz=z[j]-zin;
      dz=dz-(int)(dz*c2m1)*c;
      r2=dx*dx+dy*dy+dz*dz;
      if(r2<rc2) {
        rm6=1.0/(r2*r2*r2);
        un[j]=(4.0*rm6-4.0)*rm6;
      } else {
        un[j]=0.0;
      }
      su=su+(un[j]-umat[i][j]);
    }
  }

  /* Acceptance test */

  if(su<=0.0) {
    accept=TRUE;
  } else {
    if(su/t<alnmax) {
      if(drand48()<=exp(-su/t)) {
        accept=TRUE;
      } else {
        accept=FALSE;
      }
    } else {
      accept=FALSE;
    }
  }

  /* Accept move & update pair interaction matrix */

  if(accept) {
    nacc=nacc+1;
    x[i]=xin;
    y[i]=yin;
    z[i]=zin;
    for(j=0;j<=n-1;j++) {
      umat[i][j]=un[j];
      umat[j][i]=un[j];
    }
  }

  return;

}

/**********************************************************************/

void writePDB() {

/**********************************************************************/
    FILE *fpo;
    int i;
        if(nt==1) {
          fpo=fopen("traj.pdb","w");
    } else { 
          fpo=fopen("traj.pdb","a");
    }
    fprintf(fpo,"REMARK    GENERATED BY A COMPUTER\n");
    fprintf(fpo,"TITLE     LAYERS t= %f\n",t);
    fprintf(fpo,"REMARK    THIS IS A SIMULATION BOX\n");
    fprintf(fpo,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%10s%3d\n",
              c,c,c,
              90.,90.,90.,"P  1 ",1);
    fprintf(fpo,"MODEL         %d\n",nt);

    for(i=0;i<n;i++) { 
        fprintf(fpo,"%-6s%5d %4s%1s%3s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
                    "ATOM",
                    i+1,
                           "LJ",
                    " ",
                           "LJ", 
                    " ",
                    (i+1)%10000, // pdb accepts at most 4-digit residue numbers
                    " ",
                   x[i],
                   y[i],
                   z[i],
               0.0,
                           0.0,
                    " ",
                    "0");
    }    
        fprintf(fpo,"TER\nENDMDL\n");
    fclose(fpo);
} 

/**********************************************************************/

void writeXYZ() {

/**********************************************************************/

  FILE *fpo;
  int i;
  
  if(nt==ntprint) {
    fpo=fopen("traj.xyz", "w");
  } else { 
    fpo=fopen("traj.xyz", "a");
  }
  
  fprintf(fpo, "%d\n", n);
  
  /* extented XYZ format for use with e. g. Ovito */
  fprintf(fpo, "Lattice=\"%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f\" ", c, c, c);
  fprintf(fpo, "Properties=species:S:1:pos:R:3 ");
  fprintf(fpo, "Time=%d\n", nt);

  for(i=0;i<n;i++) { 
    fprintf(fpo, "Ar %15.8e %15.8e %15.8e\n", x[i], y[i], z[i]);
  }
  
  fclose(fpo);

}

/**********************************************************************/

void means() {

/**********************************************************************/

  int i,j,k;
  double dx,dy,dz,p,r2,rm6,su,sw;

  /* Potential energy, virial & g(r) */

  su=su0;
  sw=sw0;
  for(i=0;i<=n-2;i++) {
    for(j=i+1;j<=n-1;j++) {
      dx=x[j]-x[i];
      dx=dx-(int)(dx*c2m1)*c;
      dy=y[j]-y[i];
      dy=dy-(int)(dy*c2m1)*c;
      dz=z[j]-z[i];
      dz=dz-(int)(dz*c2m1)*c;
      r2=dx*dx+dy*dy+dz*dz;
      if(r2<rc2) {
        rm6=1.0/(r2*r2*r2);
        su=su+umat[i][j];          /* Potential energy */
        sw=sw+(24.0-48.0*rm6)*rm6; /* Virial */
      }
      if(r2<r2max) {
        k=sqrt(r2)*drm1;
        if(k<=ndr-1) {
          ag[k]=ag[k]+(long long)1;
        }
      }
    }
  }
  su=su/n;
  sw=sw/n;

  /* Print control variables and write trajectory file */

  if((ntprint>0)&&(nt%ntprint==0)) {
    p=rho*(t-sw/3.0);
    printf(" %10d %12.5e %12.5e %12.5e\n",\
       nt,(double)nacc/(n*ntskip),su,p);
    writeXYZ();
  }

  /* Accumulate averages */

  accr=accr+(double)nacc/(n*ntskip);
  au=au+su;
  au2=au2+su*su;
  aw=aw+sw;

  /* Clear acceptance counter */

  nacc=0;

  return;

}

/**********************************************************************/

void putcf() {

/**********************************************************************/

  FILE *fpo;
  unsigned short *lastseed;

  /* RNG seed */

  lastseed=seed48(seed);
  seed[0]=lastseed[0];
  seed[1]=lastseed[1];
  seed[2]=lastseed[2];

  /* Write checkpoint file */

  fpo=fopen("mclj_out.dat","w");

  fwrite(&n,sizeof(int),1,fpo);
  fwrite(&nt,sizeof(int),1,fpo);
  fwrite(&ntjob,sizeof(int),1,fpo);
  fwrite(&ntprint,sizeof(int),1,fpo);
  fwrite(&ntskip,sizeof(int),1,fpo);
  fwrite(seed,sizeof(unsigned short),3,fpo);

  fwrite(&disp,sizeof(double),1,fpo);
  fwrite(&dr,sizeof(double),1,fpo);
  fwrite(&rho,sizeof(double),1,fpo);
  fwrite(&t,sizeof(double),1,fpo);

  fwrite(x,sizeof(double),n,fpo);
  fwrite(y,sizeof(double),n,fpo);
  fwrite(z,sizeof(double),n,fpo);

  fwrite(&accr,sizeof(double),1,fpo);
  fwrite(&au,sizeof(double),1,fpo);
  fwrite(&au2,sizeof(double),1,fpo);
  fwrite(&aw,sizeof(double),1,fpo);

  fwrite(ag,sizeof(long long),ndr,fpo);

  fclose(fpo);

  /* Deallocate arrays */

  free(ag);
  free(umat);
  free(umatb);
  free(un);
  free(x);
  free(y);
  free(z);

  return;

}

/**********************************************************************/

void usage(char** argv, int help) {

/**********************************************************************/

   if(!help) { 
      exit(printf("Usage: %s [dump|help] \n",argv[0]));
   }else { 
    printf("If no options are given, the program runs the simulation\n");
    printf("according to the binary input file 'mclj_in.dat' as \n");
    printf("generated by one of the associated utilities (lmclj,gmclj).\n\n");
    printf("If the 'dump' option is specified, the binary input file is\n");
    printf("written to stdout. The 'help' option prints this information\n");
    exit(0);
   }
}

/**********************************************************************/

int main(int argc, char** argv) {

/**********************************************************************/

  if(argc==2) {
    if(!strcmp(argv[1],"dump"))  
        exit(dumpcf(argv));
    if(!strcmp(argv[1],"help")) 
            usage(argv,1);
        usage(argv,0);
  } 

  int i,j;

  /* Read startup/checkpoint file & initialize */

  getcf(argv);
  uinit();

  /* Do ntskip*ntjob passes */

  for(i=1;i<=ntjob;i++) {
    nt=nt+1;
    for(j=1;j<=ntskip*n;j++) {
      move();
    }
    means();
  }

  /* Write checkpoint file */
  
  putcf();

  return 0;

}

/**********************************************************************/
