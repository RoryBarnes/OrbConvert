/**********************************************************
 * orbconvert.c -- Convert between orbital elements and   *
 * Cartesian coordinates                                  *
 *                                                        *
 * Rory Barnes                                            *
 *                                                        *
 * Thu Oct  5 13:42:27 PDT 2017                           *
 *                                                        *
 **********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#define dot(a,b)        (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])


// Switch to Gauss' Constant?
/* Gravitation Constant; From Luzum et al., 2011; value recommended by IAU NSFA in Prsa et al. 2016 */
// Units are SI
#define BIGG       6.67428e-11       
#define MSUN       1.988416e30       // Solar Mass; Prsa et al. 2016
#define MEARTH     5.972186e24       // Earth Mass; Prsa et al. 2016
#define AUM        1.49598e11
#define PI         3.1415926535
#define RADDEG     0.017453292519444445 // Degrees per radian

// Option formatting
#define LINE       128
#define OPTLEN     24
#define MAXLINES   128
#define MAXARRAY   24

// Orbital Element Systems
#define KEPLERIAN  0
#define COMETARY   1
#define ASTEROIDAL 2
#define CARTESIAN  3

// Origin Options
#define BARY       10
#define HELIO      11
// No Jacobi for now

typedef struct {
  double dMass,dPeriDist,dEpoch,dTPeri;
  double dSemi;
  double dEcc;
  double dIncl;
  double dLasc;
  double dArgPeri;
  double dMeanAnom;
  double *x,*v;
  char name[16];
} BODY;

typedef struct {
  int iInputCoord;     // Input Coordinate System
  int iInputOrigin;    // Input Origin
  int iOutputCoord;    // Ouptut Coordinate System
  int iOutputOrigin;   // Output Origin
} OPTIONS;

/************* Options functions **********************/

void LineExit(int iLine) {
  fprintf(stderr,"\tLine %d\n",iLine);
  exit(1);
}

char *sLower(char cString[]) {
  int iPos;
  for (iPos=0;cString[iPos];iPos++)
  cString[iPos] = tolower(cString[iPos]);
    
  return cString;
}

/* Is the first non-white space a #? I so, return 1 */
int CheckComment(char cLine[],int iLen) {
  int iPos;

  for (iPos=0;iPos<iLen;iPos++) {
    if (!isspace(cLine[iPos])) {
      if (cLine[iPos] == 35) // # is ASCII code 35
      return 1;
    } else
      return 0;
  }
  return 0;
}

void GetLine(char cFile[MAXLINES][LINE],char cOption[],char cLine[],int *iLine) {
  int iLen,bDone=0,iLineTmp=0,i;
  char cWord[OPTLEN],cTmp[LINE];
  
  iLen=strlen(cOption);

  memset(cLine,'\0',LINE);
  memset(cTmp,'\0',LINE);
  memset(cWord,'\0',OPTLEN);

  i=0;
  while(!bDone) {
    strcpy(cTmp,cFile[i]);
    if (!CheckComment(cTmp,iLen)) {
      sscanf(cTmp,"%s",cWord);
      // XXX Add check for comments embedded in the option here
      if (memcmp(cWord,cOption,iLen+1) == 0) {
        /* Parameter Found! */
        if (bDone) {
	  fprintf(stderr,"Multiple occurences of parameter %s found.\n",cOption);
          fprintf(stderr,"\tlines: %d and %d\n",(*iLine+1),iLineTmp+1);
          exit(1);
        }
        strcpy(cLine,cTmp);
        *iLine=iLineTmp;
        bDone=1;
      }
    }
    iLineTmp++;
    i++;
    memset(cTmp,'\0',LINE);
    memset(cWord,'\0',OPTLEN);
  }
}

void GetNextValidLine(char cFile[MAXLINES][LINE],int iStart,char cLine[],int *iLine) {
  FILE *fp;
  int iPos,iLineTmp,ok=1;

  *iLine=0;

  for (iLineTmp=0;iLineTmp<iStart;iLineTmp++) {
    fgets(cLine,LINE,fp);
    (*iLine)++;
  }

  /* If EOF, return */
  if (fgets(cLine,LINE,fp) == NULL) {
    sprintf(cLine,"null");
    return;
  }

  /* Now check for blank line, comment (# = 35), continue ($ = 36)
     or blank line (line feed = 10). */

  for (iPos=0;iPos<LINE;iPos++) {
    if (cLine[iPos] == 36 || cLine[iPos] == 35 || cLine[iPos] == 10) {
      /* First character is a $, # or \n: continue */
      GetNextValidLine(cFile,iStart+1,cLine,iLine);
      fclose(fp);
      return;
    }
    if (!isspace(cLine[iPos])) {
      /* Found next valid line */
      fclose(fp);
      return;
    }
  }
  /* If made it here, line was blank */
  GetNextValidLine(cFile,iStart+1,cLine,iLine);
  fclose(fp);
}

/* Where is the first non-white-space character in a line? */

int GetPos(char cLine[]) {
  int iPos;

  for (iPos=0;iPos<strlen(cLine);iPos++)
    if (!isspace(cLine[iPos]))
      return iPos;

  /* Shouldn't be possible to get here */
  return 0;
}

/* Separate a line into words. cInput is an array of strings, each
   containing one word. This routine also checks if the final word has
   a trailing $, if so, it is an array that continues to the next
   line. */

void GetWords(char cLine[],char cInput[MAXARRAY][OPTLEN],int *iNumWords,int *bContinue) {
  int iPos,iPosStart,iWord;
  char cTmp[OPTLEN];

  //iPos0=GetPos(cLine);
  iWord=0;
  /* Use GetPos to avoid white space */
  for (iPos=GetPos(cLine);iPos<strlen(cLine);iPos++) {
  //for (iPos=GetPos(cLine);iPos<strlen(cLine)-GetPos(cLine);iPos++) {
    /* MEM: Something is wrong here, but it is intermittent. Sometimes a call
       here produces a memory error with valgrind. On 12/14/16 a run without the
       next print statements produced an error, but with them in, the error
       disappeared. After commenting out again, the problem was still gone. */
    /* DEBUG
    printf("%s\n",cLine);
    printf("%d %d %d\n",(int)strlen(cLine),GetPos(cLine),iPos);
    fflush(stdout);
    */
    iPosStart=0;
    while (!isspace(cLine[iPos])) {
      if (cLine[iPos] != 35) { // 35 is ASCII code for #
        /* Fill word in */
        cInput[iWord][iPosStart] = cLine[iPos];
        iPosStart++;
        iPos++;
      } else {
        /* If at the start of the word, we must decrement iWord
          so that when it is incremented at the end of the loop
          the correct number of words is returned. If at the end
          of a word, everything should be fine. */
        if (iPosStart==0) iWord--;

        iPos=strlen(cLine);
        break;
      }
    }
    /* Now advance to next word */
    while (isspace(cLine[iPos]))
      iPos++;

    iPos--;
    iWord++;
  }
  /* Is the last character a $? If so, remove it and adjust iNumWords */
  if (cInput[iWord-1][strlen(cInput[iWord-1])-1] == 36) {
    *bContinue=1;
    if (strlen(cInput[iWord-1]) == 1)
      *iNumWords = iWord-1;
    else
      *iNumWords = iWord;
    cInput[iWord-1][strlen(cInput[iWord-1])-1] = '\0';
  } else {
    *bContinue=0;
    *iNumWords=iWord;
  }
}

void AddOptionDouble(char cFile[MAXLINES][LINE],char cOption[],double *dInput,int *iLine) {
  char cTmp[OPTLEN],cLine[LINE];

  GetLine(cFile,cOption,cLine,iLine);
  if(*iLine >= 0)
      sscanf(cLine,"%s %lf",cTmp, dInput);
}

void AddOptionInt(char cFile[MAXLINES][LINE],char cOption[],int *iInput,int *iLine) {
  char cTmp[OPTLEN],cLine[LINE];

  GetLine(cFile,cOption,cLine,iLine);
  if(*iLine >= 0)
      sscanf(cLine,"%s %d",cTmp,iInput);
}

void AddOptionBool(char cFile[MAXLINES][LINE],char cOption[],int *iInput,int *iLine) {

  AddOptionInt(cFile,cOption,iInput,iLine);
  if (*iInput == 0 || *iInput == 1)
    return;
  else {
    fprintf(stderr,"ERROR: %s must be either 0 or 1.\n",cOption);
    LineExit(*iLine);
  }
}

void AddOptionString(char cFile[MAXLINES][LINE],char cOption[],char cInput[],int *iLine) {
  char cTmp[OPTLEN],cLine[LINE];

  GetLine(cFile,cOption,cLine,iLine);
  sscanf(cLine,"%s %s",cTmp,cInput);
}

/********** Read Functions *************/

void ReadSemi(char cFile[MAXLINES][LINE],BODY *body) {
  int lTmp=-1;
  double dTmp;

  AddOptionDouble(cFile,"dSemi",&dTmp,&lTmp);
  if (lTmp >= 0) {
    if (dTmp <= 0) {
      fprintf(stderr,"ERROR: Mass must be > 0.\n");
      LineExit(lTmp);
    }
    body->dSemi = dTmp*AUM;
  } else {
    printf("ERROR: Semi-major Axis is not set.\n");
  }
}

void ReadEcc(char cFile[MAXLINES][LINE],BODY *body) {
  int lTmp=-1;
  double dTmp;

  AddOptionDouble(cFile,"dEcc",&dTmp,&lTmp);
  if (lTmp >= 0) {
    if (dTmp <= 0) {
      fprintf(stderr,"ERROR: Eccentricity must be > 0.\n");
      LineExit(lTmp);
    }
    body->dEcc = dTmp;
  } else {
    printf("ERROR: Semi-major Axis is not set.\n");
  }
}

void ReadIncl(char cFile[MAXLINES][LINE],BODY *body) {
  int lTmp=-1;
  double dTmp;

  AddOptionDouble(cFile,"dIncl",&dTmp,&lTmp);
  if (lTmp >= 0) {
    if (dTmp < 0) {
      fprintf(stderr,"ERROR: Inclination must be > 0.\n");
      LineExit(lTmp);
    }
    body->dIncl = dTmp*RADDEG;
  } else {
    printf("ERROR: Semi-major Axis is not set.\n");
  }
}

void ReadLasc(char cFile[MAXLINES][LINE],BODY *body) {
  int lTmp=-1;
  double dTmp;

  AddOptionDouble(cFile,"dLasc",&dTmp,&lTmp);
  if (lTmp >= 0) {
    if (dTmp < 0) {
      fprintf(stderr,"ERROR: Longitude of Ascending Node must be > 0.\n");
      LineExit(lTmp);
    }
    body->dLasc = dTmp*RADDEG;
  } else {
    printf("ERROR: Semi-major Axis is not set.\n");
  }
}

void ReadArgPeri(char cFile[MAXLINES][LINE],BODY *body) {
  int lTmp=-1;
  double dTmp;

  AddOptionDouble(cFile,"dIncl",&dTmp,&lTmp);
  if (lTmp >= 0) {
    body->dArgPeri = dTmp*RADDEG;
  } else {
    printf("ERROR: Semi-major Axis is not set.\n");
  }
}

void ReadMeanAnom(char cFile[MAXLINES][LINE],BODY *body) {
  int lTmp=-1;
  double dTmp;

  AddOptionDouble(cFile,"dMeanAnom",&dTmp,&lTmp);
  if (lTmp >= 0) {
    body->dMeanAnom = dTmp*RADDEG;
  } else {
    printf("ERROR: Mean Anomaly is not set.\n");
  }
}

void ReadPeriDist(char cFile[MAXLINES][LINE],BODY *body) {
  int lTmp=-1;
  double dTmp;

  AddOptionDouble(cFile,"dPeriDist",&dTmp,&lTmp);
  if (lTmp >= 0) {
    body->dPeriDist = dTmp*AUM;
  } else {
    printf("ERROR: Perihelion Distance is not set.\n");
  }
}

void ReadTPeri(char cFile[MAXLINES][LINE],BODY *body) {
  int lTmp=-1;
  double dTmp;

  AddOptionDouble(cFile,"dTPeri",&dTmp,&lTmp);
  if (lTmp >= 0) {
    body->dTPeri = dTmp;
  } else {
    printf("ERROR: Time of Pericenter is not set.\n");
  }
}

void ReadEpoch(char cFile[MAXLINES][LINE],BODY *body) {
  int lTmp=-1;
  double dTmp;

  AddOptionDouble(cFile,"dEpoch",&dTmp,&lTmp);
  if (lTmp >= 0) {
    body->dEpoch = dTmp;
  } else {
    printf("ERROR: Epoch is not set.\n");
  }
}





void ReadInputCoord(char cFile[MAXLINES][LINE],OPTIONS *options) {
  int iLine;
  char cTmp[OPTLEN];

  AddOptionString(cFile,"sInputCoord",cTmp,&iLine);

  if (memcmp(sLower(cTmp),"k",1) == 0) {
    options->iInputCoord = KEPLERIAN;
    printf("Keplerian Input Elements\n");
  } else if (memcmp(sLower(cTmp),"co",2) == 0) {
    options->iInputCoord = COMETARY;
  } else if (memcmp(sLower(cTmp),"a",1) == 0) {
    options->iInputCoord = ASTEROIDAL;
  } else if (memcmp(sLower(cTmp),"ca",2) == 0) {
    options->iInputCoord = CARTESIAN;
  } else {
    fprintf(stderr,"ERROR: Unknown argument to sInputCoord: %s. Options are Keplerian, Cometary, or Asteroidal.\n",cTmp);
    LineExit(iLine);
  }
}

void ReadInputOrigin(char cFile[MAXLINES][LINE],OPTIONS *options) {
  int iLine;
  char cTmp[OPTLEN];

  AddOptionString(cFile,"sInputOrigin",cTmp,&iLine);

  if (memcmp(sLower(cTmp),"h",1) == 0) {
    options->iInputOrigin = HELIO;
  } else if (memcmp(sLower(cTmp),"b",1) == 0) {
    options->iInputCoord = BARY;
  } else {
    fprintf(stderr,"ERROR: Unknown argument to sInputOrigin: %s. Options are Heliocentric or Barycentric.\n",cTmp);
    LineExit(iLine);
  }
}

void ReadOutputCoord(char cFile[MAXLINES][LINE],OPTIONS *options) {
  int iLine;
  char cTmp[OPTLEN];

  AddOptionString(cFile,"sOutputCoord",cTmp,&iLine);

  if (memcmp(sLower(cTmp),"k",1) == 0) {
    options->iOutputCoord = KEPLERIAN;
    //printf("Keplerian Input Elements\n");
  } else if (memcmp(sLower(cTmp),"co",2) == 0) {
    options->iOutputCoord = COMETARY;
  } else if (memcmp(sLower(cTmp),"a",1) == 0) {
    options->iOutputCoord = ASTEROIDAL;
  } else if (memcmp(sLower(cTmp),"ca",2) == 0) {
    options->iOutputCoord = CARTESIAN;  
  } else {
    fprintf(stderr,"ERROR: Unknown argument to sOututCoord: %s. Options are Keplerian, Cometary, or Asteroidal.\n",cTmp);
    LineExit(iLine);
  }
}

void ReadOutputOrigin(char cFile[MAXLINES][LINE],OPTIONS *options) {
  int iLine;
  char cTmp[OPTLEN];

  AddOptionString(cFile,"sOutputOrigin",cTmp,&iLine);

  if (memcmp(sLower(cTmp),"h",1) == 0) {
    options->iOutputOrigin = HELIO;
    //printf("Keplerian Input Elements\n");
  } else if (memcmp(sLower(cTmp),"b",1) == 0) {
    options->iOutputCoord = BARY;
  } else {
    fprintf(stderr,"ERROR: Unknown argument to sOutputOrigin: %s. Options are Heliocentric or Barycentric.\n",cTmp);
    LineExit(iLine);
  }
}

void ReadMass(char cFile[MAXLINES][LINE],BODY *body,OPTIONS *options) {
  int lTmp=-1;
  double dTmp;

  AddOptionDouble(cFile,"dMass",&dTmp,&lTmp);
  if (lTmp >= 0) {
    if (dTmp < 0) {
      fprintf(stderr,"ERROR: Mass must be >= 0.\n");
      LineExit(lTmp);
    }
    body->dMass = dTmp*MEARTH;
  }
}

void ReadKeplerElems(char cFile[MAXLINES][LINE],BODY *body,OPTIONS *options) {

  ReadSemi(cFile,body);
  ReadEcc(cFile,body);
  ReadIncl(cFile,body);
  ReadLasc(cFile,body);
  ReadArgPeri(cFile,body);
  ReadMeanAnom(cFile,body);
}

void ReadCometElems(char cFile[MAXLINES][LINE],BODY *body,OPTIONS *options) {

  ReadSemi(cFile,body);
  ReadPeriDist(cFile,body);
  ReadIncl(cFile,body);
  ReadLasc(cFile,body);
  ReadArgPeri(cFile,body);
  ReadTPeri(cFile,body);
}

void ReadCartes(char cFile[MAXLINES][LINE],BODY *body,OPTIONS *options){
}

void ReadOptions(BODY *body,OPTIONS *options,char infile[]) {
  /* This function reads in the input data */

  char *cTmp,cFile[MAXLINES][LINE];
  FILE *fp;
  int iLine=0;
  size_t iLen;

  fp = fopen(infile,"r");

  /* Fill cFile with all the lines from the input file. */
  //getline(&cTmp,&iLen,fp);
  while (getline(&cTmp,&iLen,fp) != -1) {
    strcpy(cFile[iLine++],cTmp);
    //iLine++;
  }

  fclose(fp);

  ReadInputCoord(cFile,options);
  //exit(1);
  ReadInputOrigin(cFile,options);
  ReadOutputCoord(cFile,options);
  ReadOutputOrigin(cFile,options);

  ReadMass(cFile,body,options);
  if (options->iInputCoord == KEPLERIAN)
    ReadKeplerElems(cFile,body,options);
  else if (options->iInputCoord == COMETARY) 
    ReadCometElems(cFile,body,options);
  //ReadCartes(cFile,body,options);
}

/**************** Orbit Conversion Subroutines **********/

void MassFunctions(double *mass,double *ms,double *mu) {
  int j;
  double ksq;

  ksq=BIGG*mass[0];

  ms[0]=mass[0];
  for (j=1;j<=2;j++) {
      ms[j]=ms[j-1]+mass[j];
      mu[j]=(ksq*ms[j])/ms[j-1];
  }     
}

/* Astrocentric -> Barycentric Cartesian coordinates */
void hel_bar(double **hel,double **bar,double *m,double *ms,int P)
{
   int i,p;

   for(i= 0;i<3;i++)
      bar[0][i] = 0;
   for(p= 1;p<=P;p++)
      for(i=0;i<3;i++)
         bar[0][i] -= m[p]/ms[P]*hel[p][i];
   for(p= 1;p<=P;p++)
      for(i=0;i<3;i++)
         bar[p][i] = hel[p][i] + bar[0][i];
}

/* Barycentric -> Astrocentric Cartesian coordinates */
void bar_hel(double **bar,double **hel,int P)
{
   int i,p;

   for (p=1;p<=P;p++)
      for (i=0;i<3;i++) hel[p][i] = bar[p][i] - bar[0][i];
   for (i=0;i<3;i++) hel[0][i] = 0;

}

/* Calculate the total angular momentum vector */
double* invariable_plane(double **x,double **v,int P,double *m) {
  int p,i;
  double *dh,*h;
  double mg;

  dh=malloc(3*sizeof(double));
  h=malloc(3*sizeof(double));

  cross(h,x[0],v[0]);
  for (i=0;i<3;i++)
    h[i] *= m[0];   
  for(p=1;p<=P;p++) {
    cross(dh,x[p],v[p]);
    for(i=0;i<3;i++)
      h[i]+= m[p]*dh[i];
  }
  mg= pow(dot(h,h),0.5);
  for(i=0;i<3;i++) 
    h[i]/= mg;
  return h;
}

/* Rotate coordinates */
void rotate(double *z,double **x,int np)  {
    double phi, theta;
    int i,k;
    double **x1;
    
    x1=malloc((np+1)*sizeof(double*));
    for (i=0;i<=np;i++)
        x1[i]=malloc(3*sizeof(double));

    theta = atan2(z[1],z[0]);
    phi = atan2(sqrt(z[0]*z[0] + z[1]*z[1]),z[2]);

    /* Rotate about z-axis */
    for (i=0;i<=np;i++) {
      x1[i][0] = x[i][0]*cos(theta) + x[i][1]*sin(theta);
      x1[i][1] = -x[i][0]*sin(theta) + x[i][1]*cos(theta);
      x1[i][2] = x[i][2];
    }
    
    /* Rotate about new y-axis (z -> x, x -> y) */
    for (i=0;i<=np;i++) {
      x[i][0] = -x1[i][2]*sin(phi) + x1[i][0]*cos(phi);
      x[i][1] = x1[i][1];
      x[i][2] = x1[i][2]*cos(phi) + x1[i][0]*sin(phi);
    }

    free(x1);
}


double GetEccAnom(double m,double e) {
  double es,ec,w,wp,wpp,wppp,ecc,dx,lo,up,next;
  int iter;

  if(sin(m)>0)
    ecc= m+0.85*e;
  else ecc= m-0.85*e;
  lo = -2*PI;
  up = 2*PI;
  for(iter=1;iter<=32;iter++) {   
    es= e*sin(ecc);
    ec= e*cos(ecc);
    w= ecc-es-m;
    wp= 1-ec;
    wpp= es;
    wppp= ec;
    if(w>0)up= ecc;
    else lo= ecc;
    dx= -w/wp;
    dx= -w/(wp+dx*wpp/2);
    dx= -w/(wp+dx*wpp/2+dx*dx*wppp/6);
    next= ecc+dx;
    if(ecc==next)
        break;
    if((next>lo)&&(next<up))
      ecc= next;
    else 
      ecc= (lo+up)/2;
    if((ecc==lo)||(ecc==up))
        break;
    if(iter>30)
      printf("%4d %23.20f %e\n",iter,ecc,up-lo);
  }
  if(iter>32) {
    printf("Kepler soln failed\n");
    exit(1);
  }
    
  return ecc;
}

void cartes(BODY *body,double mu) {
  double a,e,m,cosi,sini,cos_lasc,sin_lasc,cos_aper,sin_aper;
  double ecc,sin_ecc,cos_ecc,l1,m1,n1,l2,m2,n2;
  double xi,eta,vel_scl;
 
  a= body->dSemi;
  e= body->dEcc;  
  m= body->dMeanAnom;
  cosi= cos(body->dIncl);
  sini= sin(body->dIncl);
  cos_lasc= cos(body->dLasc);
  sin_lasc= sin(body->dLasc);
  cos_aper= cos(body->dArgPeri);
  sin_aper= sin(body->dArgPeri);
 
  /*
   * Reduce mean anomoly to [0, 2*PI) 
   */
  m-= ((int)(m/(2*PI)))*2*PI;
  /*Solve Kepler's equation. */
  ecc = GetEccAnom(m,e);

  cos_ecc= cos(ecc);
  sin_ecc= sin(ecc);
    
  l1= cos_lasc*cos_aper-sin_lasc*sin_aper*cosi;
  m1= sin_lasc*cos_aper+cos_lasc*sin_aper*cosi;
  n1= sin_aper*sini;
  l2= -cos_lasc*sin_aper-sin_lasc*cos_aper*cosi;
  m2= -sin_lasc*sin_aper+cos_lasc*cos_aper*cosi;
  n2= cos_aper*sini;
    
  xi= a*(cos_ecc-e);
  eta= a*pow(1-e*e,0.5)*sin_ecc;
  body->x[0]= l1*xi+l2*eta;
  body->x[1]= m1*xi+m2*eta;
  body->x[2]= n1*xi+n2*eta;
  
  vel_scl= pow((mu*a)/dot(body->x,body->x),0.5);
  xi= -vel_scl*sin_ecc;
  eta= vel_scl*pow(1-e*e,0.5)*cos_ecc;
  body->v[0]= l1*xi+l2*eta;
  body->v[1]= m1*xi+m2*eta;
  body->v[2]= n1*xi+n2*eta;

}

/* Compute orbital elements from Cartesian coordinates */
void elems(BODY *body,double mu) {
   double hx,hy,hz,hsq,hxy,h,r,vsq,rdot;
   double sin_lasc,cos_lasc,sin_aperf,cos_aperf,sinf,cosf;
   double a,e,sin_ecc,cos_ecc,arg;
   int i;
/* 
 * Compute various intermediate quantities.
 */
   hx= (double)(body->x[1]*body->v[2] - body->x[2]*body->v[1]);
   hy= (double)(body->x[2]*body->v[0] - body->x[0]*body->v[2]);
   hz= (double)(body->x[0]*body->v[1] - body->x[1]*body->v[0]);
   hsq= hx*hx+hy*hy+hz*hz;
   h= sqrt(hsq);
   hxy= sqrt(hx*hx+hy*hy);
   r= sqrt((double)dot(body->x,body->x));
   vsq= (double)dot(body->v,body->v);
   rdot= (double)dot(body->x,body->v);
   a= 1/(2/r-vsq/mu);
   e= sqrt(1-hsq/(mu*a));

   sin_lasc= hx/hxy;
   cos_lasc= -hy/hxy;
   sin_aperf= body->x[2]*h/(hxy*r);
   cos_aperf= (body->x[0]*cos_lasc+body->x[1]*sin_lasc)/r;
   cosf= (hsq/(mu*r)-1)/e;
   if(fabs(cosf)>1) {
     if(cosf>0) 
       cosf= 1;
     else
       cosf= -1; 
   }
   sinf= sqrt(1-cosf*cosf);
   if(rdot<0)  
     sinf= -sinf;
   cos_ecc= (1-r/a)/e;
   if (fabs(cos_ecc)>1) {
      if(cos_ecc>0)
        cos_ecc= 1;
      else
        cos_ecc= -1;
   }
   sin_ecc = sqrt(1-cos_ecc*cos_ecc);
   if (rdot<0) 
     sin_ecc = -sin_ecc;

   body->dSemi=a;
   body->dEcc=e;
   body->dIncl=atan2(hxy/h,hz/h);
   body->dLasc= atan2(sin_lasc,cos_lasc);
   arg= atan2(sin_aperf,cos_aperf)-atan2(sinf,cosf);
   while(arg>M_PI) 
     arg-= 2*M_PI;
   while(arg<-M_PI) 
     arg+= 2*M_PI;
   body->dArgPeri = arg;
   body->dMeanAnom = atan2(sin_ecc,cos_ecc)-e*sin_ecc;
}

double dSemiToPeriod(double a,double ms) {
  return pow(4*PI*PI*a*a*a/(BIGG*ms),0.5);
}

void ConvertCometToKepler(BODY *body,double dMass) {
  double dPeriod;

  dPeriod = dSemiToPeriod(body->dSemi,dMass);

  body->dEcc = 1 - body->dPeriDist/body->dSemi;
  body->dMeanAnom = 2*PI/dPeriod*body->dEpoch - body->dTPeri;
}

void WriteOutput(BODY *body,OPTIONS *options) {
  if (options->iOutputCoord == CARTESIAN) {
    printf("Position = (%.3e,%.3e,%.3e)\n",body->x[0],body->x[1],body->x[2]);
    printf("Velocity = (%.3e,%.3e,%.3e)\n",body->v[0],body->v[1],body->v[2]);
    printf("Distance = %.3e\n",sqrt(dot(body->x,body->x)));
    printf("Velocity = %.3e\n",sqrt(dot(body->v,body->v)));
  }
}

int main(int argc, char *argv[]) {
  BODY body;
  OPTIONS options;
  double *mass,*ms,*mu,lng;
  double **x,**v,mstar;
  int numpl,numpts;
  int j,k;
  char name[16];
  FILE **fp,*ifp;

  body.x = malloc(3*sizeof(double));
  body.v = malloc(3*sizeof(double));

  if (argc != 2) {
      (void) fprintf(stderr,"Usage: %s file\n",argv[0]);
      exit(1);
  }

  ReadOptions(&body,&options,argv[1]);

  //printf("%d\n",options.iInputCoord);

  mu=malloc(2*sizeof(double));
  mass=malloc(2*sizeof(double));
  ms=malloc(2*sizeof(double));
  mu[0]=MSUN;
  mass[0]=MSUN;

  MassFunctions(mass,ms,mu);

  if (options.iInputCoord == COMETARY)
    ConvertCometToKepler(&body,(body.dMass+mass[0]));

  /* Convert to Cartesian astrocentric coordinates */
  cartes(&body,mu[1]);

  WriteOutput(&body,&options);
}
