/* 
   printtipsy
   Utility for examining contents of tipsy binary files
   (c) 2012-2021 James Wadsley

   Make:  gcc printtipsy.c -lm -o printtipsy
   Usage: printtipsy [[-av|-af] -a] [-p nprint (def 10)] [-i istart] [-e iend] [-n nth] [-gds] [-sphere x y z r] [-nc] [-prop] filename
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <rpc/xdr.h>
#include <rpc/rpc.h>
#include <assert.h>
#include <unistd.h>

#define MAXDIM 3
#define forever for(;;)

typedef float Real;

struct sphere { Real x,y,z, r2; };
struct sphere s;

struct gas_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real rho;
    Real temp;
    Real eps;
    Real metals ;
    Real phi ;
    } ;

struct gas_particle *gas_particles;

struct dark_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real eps;
    Real phi ;
    } ;

struct dark_particle *dark_particles;

struct star_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real metals ;
    Real tform ;
    Real eps;
    Real phi ;
    } ;

struct star_particle *star_particles;

struct dump {
    double time ;
    int nbodies ;
    int ndim ;
    int nsph ;
    int ndark ;
    int nstar ;
    } ;

#define Real float

int xdr_header(struct dump *header, XDR xdrlocal)
    {
    int pad;
  
    if(xdr_double(&xdrlocal, &header->time) != TRUE)
        return 0;
    if(xdr_int(&xdrlocal, &header->nbodies) != TRUE)
        return 0;
    if(xdr_int(&xdrlocal, &header->ndim) != TRUE)
        return 0;
    if(xdr_int(&xdrlocal, &header->nsph) != TRUE)
        return 0;
    if(xdr_int(&xdrlocal, &header->ndark) != TRUE)
        return 0;
    if(xdr_int(&xdrlocal, &header->nstar) != TRUE)
        return 0;
    if(xdr_int(&xdrlocal, &pad) != TRUE)
        return 0;
    return 1;
    }

double xcom[3],v2com[3],vcom[3],L[3],mass;
int nSph=-1,nDark=-1,nStar=-1,bReadAux=0;

int PrintTipsy(char *file,int bDark,int bGas,int bStar,int nth,int istart,int iend, int bCheck, int bProp, int bAux, int iAuxFormat, int bSphere)
    {
    int i,bStd=0;
    long offset;
    struct dump hread;
    struct gas_particle gp;
    struct dark_particle dp;
    struct star_particle sp;
    char stdtext[]="STD \0";
    XDR xdrread;
    FILE *fpread;
    struct stat fileinfo;
    int   fildes,status;
    long long size;
    int npartnat,npartstd;

    fildes = open(file, O_RDONLY);
    if (fildes < 0) {
        if (bAux) return 0;
        fprintf(stderr,"Can't open file: %s\n",file);
        exit(1);
        }
    status = fstat(fildes, &fileinfo);
    size = fileinfo.st_size;
    close(fildes);
	
    fpread = fopen( file, "r" );
    if (fpread == NULL) {
        if (bAux) return 0;
        fprintf(stderr,"Can't open file: %s\n",file);
        exit(1);
        }

    fread(&hread,28,1,fpread);
    if (hread.ndim !=3) {
        npartnat = *((int *) &hread);
        bStd = 1;
        rewind(fpread);
        xdrstdio_create(&xdrread, fpread, XDR_DECODE);
        xdr_header( &hread, xdrread );   
        if (bCheck && size != 32+sizeof(gp)*hread.nsph+sizeof(dp)*hread.ndark+sizeof(sp)*hread.nstar) {
            int DataSize=sizeof(float);
            union { float f; int i; } Data;
            // See if it isn't tipsy format but a tipsy array file
            rewind(fpread);
            xdr_int( &xdrread, &npartstd );
            bStd=-1;
            if (npartnat == (size-4)/4) {
                fprintf(stderr,"%d Particles in Native Binary array file %s\n",npartnat,file);
                bStd=0;
                }
            else if (npartstd == (size-4)/4) {
                fprintf(stderr,"%d Particles in STD Binary array file %s\n",npartstd,file);
                bStd=1;
                }
            if (bStd>=0) {
                int iendg = iend;
                bReadAux=1; // flag that its not a full tipsy file 
                if (istart) {
                    offset = DataSize*istart;
                    fseek(fpread, offset, SEEK_CUR );
                    }
                if (bGas && nSph > 0) {
                    if (iendg >= nSph) iendg = nSph;
                    if (iAuxFormat == 1)
                        printf("GAS  ");
                    else
                        printf("GAS        Float        Int \n");
                    }
                else 
                    if (iAuxFormat == 1)
                        printf("PARTICLES ");
                    else
                        printf("           Float        Int \n");
                for (i=istart;i<=iendg;++i) {
                    if (bStd) 
                        xdr_float(&xdrread,&(Data.f));
                    else 
                        fread(&Data,DataSize,1,fpread);
                    if ((i%nth)==0 && i>=istart && i<=iend) if (iAuxFormat == 1) printf("%g ",Data.f); else printf("%10d %g  %d \n",i,Data.f,Data.i);
                    }
                if (iAuxFormat == 1) printf("\n");
                
                if (bDark && nDark > 0) {
                    int iendd = (iend < nDark ? iend : nDark-1);
                    if (iAuxFormat == 1)
                        printf("DARK ");
                    else
                        printf("DARK       Float        Int \n");
                    rewind(fpread);
                    offset = DataSize*(nSph > 0 ? nSph : 0)+4;
                    fseek(fpread, offset, SEEK_CUR );
                    for (i=0;i<=iendd;++i) {
                        if (bStd) 
                            xdr_float(&xdrread,&(Data.f));
                        else 
                            fread(&Data,DataSize,1,fpread);
                        if ((i%nth)==0) if (iAuxFormat == 1) printf("%g ",Data.f); else printf("%10d %g  %d \n",i,Data.f,Data.i);
                        }
                    if (iAuxFormat == 1) printf("\n");
                    }
                if (bStar && nStar > 0) {
                    int iends = (iend < nStar ? iend : nStar-1);
                    if (iAuxFormat == 1)
                        printf("STAR ");
                    else
                        printf("STAR       Float        Int \n");
                    rewind(fpread);
                    offset = DataSize*((nSph > 0 ? nSph : 0)+(nDark > 0 ? nDark : 0))+4;
                    fseek(fpread, offset, SEEK_CUR );
                    for (i=0;i<=iends;++i) {
                        if (bStd) 
                            xdr_float(&xdrread,&(Data.f));
                        else 
                            fread(&Data,DataSize,1,fpread);
                        if ((i%nth)==0) if (iAuxFormat == 1) printf("%g ",Data.f); else printf("%10d %g  %d \n",i,Data.f,Data.i);
                        }
                    if (iAuxFormat == 1) printf("\n");
                    }
                return (bStd ? npartstd : npartnat);
                }
            // No idea what it is
            fprintf(stderr,"Std file size nuts: %lld != %lld\n",(long long int)size,(long long int) (32+sizeof(gp)*hread.nsph+sizeof(dp)*hread.ndark+sizeof(sp)*hread.nstar));
            exit(1);
            }
        }
    else if (bCheck && size == 32+sizeof(gp)*hread.nsph+sizeof(dp)*hread.ndark+sizeof(sp)*hread.nstar) {
        /* Read 4 byte pad */
        assert(sizeof(int)==4);
        fread(&i,4,1,fpread);
        }
    else if (bCheck && size != 28+sizeof(gp)*hread.nsph+sizeof(dp)*hread.ndark+sizeof(sp)*hread.nstar) {
        fprintf(stderr,"File size nuts: %lld != %lld or %lld\n",(long long int) size,(long long int) (28+sizeof(gp)*hread.nsph+sizeof(dp)*hread.ndark+sizeof(sp)*hread.nstar),(long long int) (32+sizeof(gp)*hread.nsph+sizeof(dp)*hread.ndark+sizeof(sp)*hread.nstar));
        exit(1);
        }
           
    fprintf(stderr,"  Time %f Particles in %sTIPSY file: %d gas, %d dark, %d stars\n",hread.time,stdtext+(bStd ? 0 : 4),hread.nsph,hread.ndark,hread.nstar);
//	fprintf(stderr,"%ld %ld %ld %ld\n",size,sizeof(gp)*hread.nsph,sizeof(dp)*hread.ndark,sizeof(sp)*hread.nstar);

    if (bProp) {
        xcom[0]=0;	    xcom[1]=0;	    xcom[2]=0;
        vcom[0]=0;	    vcom[1]=0;	    vcom[2]=0;
        v2com[0]=0;	    v2com[1]=0;	    v2com[2]=0;
        mass = 0;
        L[0]=0;	    L[1]=0;	    L[2]=0;
        }
	    

    /*
    ** Read Stuff!
    */
    nSph = hread.nsph;
    nDark = hread.ndark;
    nStar = hread.nstar;
    if (bGas && hread.nsph) {
        if (istart) {
            if (istart > hread.nsph) offset = sizeof(gp)*hread.nsph;
            else offset = sizeof(gp)*istart;
            fseek(fpread, offset, SEEK_CUR );
            }
        printf("       GAS     mass         x          y          z             vx         vy         vz         dens         temp      eps          Z        phi  \n");
        for (i=istart;i<hread.nsph;++i) {
            if (i > iend && !bProp) {
                offset = sizeof(gp)*(hread.nsph-i);
                fseek(fpread, offset, SEEK_CUR );
                break;
                }
            if (bStd) 
                xdr_vector(&xdrread,(char *) &gp, (sizeof(gp)/sizeof(Real)),
                    sizeof(Real),(xdrproc_t) xdr_float);
            else 
                fread(&gp,sizeof(struct gas_particle),1,fpread);

            if (bSphere) {
                Real dx, dy, dz;	
                dx = gp.pos[0]-s.x; dy = gp.pos[1]-s.y; dz=gp.pos[2]-s.z; 
                if (dx*dx+dy*dy+dz*dz > s.r2) goto skipgp;
                }
            if (bProp) {
                xcom[0]+=gp.mass*gp.pos[0];
                xcom[1]+=gp.mass*gp.pos[1];
                xcom[2]+=gp.mass*gp.pos[2];
                vcom[0]+=gp.mass*gp.vel[0];
                vcom[1]+=gp.mass*gp.vel[1];
                vcom[2]+=gp.mass*gp.vel[2];
                v2com[0]+=gp.mass*gp.vel[0]*gp.vel[0];
                v2com[1]+=gp.mass*gp.vel[1]*gp.vel[1];
                v2com[2]+=gp.mass*gp.vel[2]*gp.vel[2];
                mass += gp.mass;
                L[0]+=gp.mass*(gp.pos[1]*gp.vel[2]-gp.pos[2]*gp.vel[1]);
                L[1]+=gp.mass*(gp.pos[2]*gp.vel[0]-gp.pos[1]*gp.vel[2]);
                L[2]+=gp.mass*(gp.pos[0]*gp.vel[1]-gp.pos[1]*gp.vel[0]);
                }
            if ((i%nth)==0 && i>=istart && i<=iend) {
                printf("%10d %10g  %10g %10g %10g  %10g %10g %10g  %10g %10g %10g %10g %10g\n",i,gp.mass,gp.pos[0],gp.pos[1],gp.pos[2],gp.vel[0],gp.vel[1],gp.vel[2],gp.rho,gp.temp,gp.eps,gp.metals,gp.phi);
                if (nth > 1) {
                    int iadd;
                    if (i+nth < hread.nsph) iadd=nth-1;
                    else iadd = hread.nsph-i-1;
                    offset = sizeof(gp)*iadd;
                    i+=iadd;
                    fseek(fpread, offset, SEEK_CUR );
                    }
                }
        skipgp:
            1;
            }
        }
    else {
        offset = sizeof(gp)*hread.nsph;
        fseek(fpread, offset, SEEK_CUR );
        }
	    
    if (bDark && hread.ndark) {
        if (istart) {
            if (istart > hread.ndark) offset = sizeof(dp)*hread.ndark;
            else offset = sizeof(dp)*istart;
            fseek(fpread, offset, SEEK_CUR );
            }
        printf("      DARK     mass         x          y          z             vx         vy         vz         eps          phi  \n");

        for (i=istart;i<hread.ndark;++i) {
            if (i > iend && !bProp) {
                offset = sizeof(dp)*(hread.ndark-i);
                fseek(fpread, offset, SEEK_CUR );
                break;
                }
            if (bStd) 
                xdr_vector(&xdrread,(char *) &dp, (sizeof(dp)/sizeof(Real)),
                    sizeof(Real),(xdrproc_t) xdr_float);
            else 
                fread(&dp,sizeof(struct dark_particle),1,fpread);

            if (bSphere) {
                Real dx, dy, dz;	
                dx = dp.pos[0]-s.x; dy = dp.pos[1]-s.y; dz=dp.pos[2]-s.z; 
                if (dx*dx+dy*dy+dz*dz > s.r2) goto skipdp;
                }

            if (bProp) {
                xcom[0]+=dp.mass*dp.pos[0];
                xcom[1]+=dp.mass*dp.pos[1];
                xcom[2]+=dp.mass*dp.pos[2];
                vcom[0]+=dp.mass*dp.vel[0];
                vcom[1]+=dp.mass*dp.vel[1];
                vcom[2]+=dp.mass*dp.vel[2];
                v2com[0]+=dp.mass*dp.vel[0]*dp.vel[0];
                v2com[1]+=dp.mass*dp.vel[1]*dp.vel[1];
                v2com[2]+=dp.mass*dp.vel[2]*dp.vel[2];
                mass += dp.mass;
                L[0]+=dp.mass*(dp.pos[1]*dp.vel[2]-dp.pos[2]*dp.vel[1]);
                L[1]+=dp.mass*(dp.pos[2]*dp.vel[0]-dp.pos[1]*dp.vel[2]);
                L[2]+=dp.mass*(dp.pos[0]*dp.vel[1]-dp.pos[1]*dp.vel[0]);
                }
            if ((i%nth)==0 && i>=istart && i<=iend) {
                printf("%10d %10g  %10g %10g %10g  %10g %10g %10g  %10g %10g\n",i,dp.mass,dp.pos[0],dp.pos[1],dp.pos[2],dp.vel[0],dp.vel[1],dp.vel[2],dp.eps,dp.phi);
                if (nth > 1) {
                    int iadd;
                    if (i+nth < hread.ndark) iadd=nth-1;
                    else iadd = hread.ndark-i-1;
                    offset = sizeof(dp)*iadd;
                    i+=iadd;
                    fseek(fpread, offset, SEEK_CUR );
                    }
                }
        skipdp:
            1;
            }
        }
    else {
        offset = sizeof(dp)*hread.ndark;
        fseek(fpread, offset, SEEK_CUR );
        }

    if (bStar && hread.nstar) {
        if (istart) {
            if (istart > hread.nstar) offset = sizeof(sp)*hread.nstar;
            else offset = sizeof(sp)*istart;
            fseek(fpread, offset, SEEK_CUR );
            }
        printf("      STAR     mass         x          y          z             vx         vy         vz           Z          tform     eps         phi  \n");
        for (i=istart;i<hread.nstar;++i) {
            if (i > iend && !bProp) {
                offset = sizeof(sp)*(hread.nstar-i);
                fseek(fpread, offset, SEEK_CUR );
                break;
                }
            if (bStd) 
                xdr_vector(&xdrread,(char *) &sp, (sizeof(sp)/sizeof(Real)),
                    sizeof(Real),(xdrproc_t) xdr_float);
            else 
                fread(&sp,sizeof(struct star_particle),1,fpread);
            
            if (bSphere) {
                Real dx, dy, dz;	
                dx = dp.pos[0]-s.x; dy = dp.pos[1]-s.y; dz=dp.pos[2]-s.z; 
                if (dx*dx+dy*dy+dz*dz > s.r2) goto skipsp;
                }

            if (bProp) {
                xcom[0]+=sp.mass*sp.pos[0];
                xcom[1]+=sp.mass*sp.pos[1];
                xcom[2]+=sp.mass*sp.pos[2];
                vcom[0]+=sp.mass*sp.vel[0];
                vcom[1]+=sp.mass*sp.vel[1];
                vcom[2]+=sp.mass*sp.vel[2];
                v2com[0]+=sp.mass*sp.vel[0]*sp.vel[0];
                v2com[1]+=sp.mass*sp.vel[1]*sp.vel[1];
                v2com[2]+=sp.mass*sp.vel[2]*sp.vel[2];
                mass += sp.mass;
                L[0]+=sp.mass*(sp.pos[1]*sp.vel[2]-sp.pos[2]*sp.vel[1]);
                L[1]+=sp.mass*(sp.pos[2]*sp.vel[0]-sp.pos[1]*sp.vel[2]);
                L[2]+=sp.mass*(sp.pos[0]*sp.vel[1]-sp.pos[1]*sp.vel[0]);
                }
            if ((i%nth)==0 && i>=istart && i<=iend) {
                printf("%10d %10g  %10g %10g %10g  %10g %10g %10g  %10g %10g %10g %10g\n",i,sp.mass,sp.pos[0],sp.pos[1],sp.pos[2],sp.vel[0],sp.vel[1],sp.vel[2],sp.metals,sp.tform,sp.eps,sp.phi);
                if (nth > 1) {
                    int iadd;
                    if (i+nth < hread.nstar) iadd=nth-1;
                    else iadd = hread.nstar-i-1;
                    offset = sizeof(sp)*iadd;
                    i+=iadd;
                    fseek(fpread, offset, SEEK_CUR );
                    }
                }
        skipsp:
            1;
            }
        }
    else {
/* Not necessary */
/*  offset = sizeof(sp)*hread.nstar;
    fseek(fpread, offset, SEEK_CUR );*/
        }

    fclose(fpread);
	
    return (hread.nbodies);
    }

void usage() {
    fprintf(stderr,"Usage: printtipsy [-af] [-av] [-a CSV] [-p nprint (def 10)] [-i istart] [-e iend] [-n nth] [-gds] [-sphere x y z r] [-nc] [-prop] filename\n");
    fprintf(stderr,"       -a Try to open files with CSV list of added extensions, e.g. -a HI,HeI,HeII");
    exit(1);
    }

int main(int argc, char **argv) {
    int i,bGas=0,bDark=0,bStar=0,bCheck=1,bProp=0,bSphere=0;
    int iAuxFormat=0;
    char *auxExt = NULL;
    int nth=1,istart=0,iend=0,nprint=0;
    char *pch, *file=NULL;
    
    i=1;
    while(i < argc) {
        pch = argv[i];
        if (!strcmp(argv[i],"-n")) {
            ++i;
            if (i >= argc) usage();
            nth = atoi(argv[i]);
            ++i;
            }
        else if (!strcmp(argv[i],"-sphere")) {
            bSphere =  1;
            ++i;
            if (i >= argc) usage();
            s.x = atof(argv[i]);
            ++i;
            if (i >= argc) usage();
            s.y = atof(argv[i]);
            ++i;
            if (i >= argc) usage();
            s.z = atof(argv[i]);
            ++i;
            if (i >= argc) usage();
            s.r2 = atof(argv[i]);
            s.r2 = s.r2*s.r2;
            ++i;
            }
        else if (!strcmp(argv[i],"-af")) {
            iAuxFormat |= 1;
            ++i;
            }
        else if (!strcmp(argv[i],"-av")) {
            iAuxFormat |= 2;
            ++i;
            }
        else if (!strcmp(argv[i],"-a")) {
            ++i;
            if (i >= argc) usage();
            auxExt = argv[i];
            ++i;
            }
        else if (!strcmp(argv[i],"-prop")) {
            bProp=1;
            ++i;
            }
        else if (!strcmp(argv[i],"-nc")) {
            bCheck=0;
            ++i;
            }
        else if (!strcmp(argv[i],"-i")) {
            ++i;
            if (i >= argc) usage();
            istart = atoi(argv[i]);
            ++i;
            }
        else if (!strcmp(argv[i],"-e")) {
            ++i;
            if (i >= argc) usage();
            iend = atoi(argv[i]);
            ++i;
            }
        else if (!strcmp(argv[i],"-p")) {
            ++i;
            if (i >= argc) usage();
            nprint = atoi(argv[i]);
            ++i;
            }
        else if (*pch == '-') {
            pch++;
            while (*pch) {
                if (*pch == 'g') bGas = 1;
                else if (*pch == 'd') bDark = 1;
                else if (*pch == 's') bStar = 1;
                pch++;
                }
            i++;
            }
        else if (file == NULL) {
            file = argv[i];
            i++;
            }
        else {
            usage();
            }
        }
    
    if (file == NULL) usage();
    if (!bGas && !bStar && !bDark) {
        bGas = 1; bStar = 1; bDark = 1;
        }
    if (auxExt != NULL && iAuxFormat == 0) iAuxFormat = 1;
    if (!nprint) nprint=10;
    if (!iend) iend = istart+(nprint-1)*nth;

    PrintTipsy( file, bDark,bGas,bStar, nth, istart, iend, bCheck, bProp, 0, iAuxFormat, bSphere );

    if (bProp) {
        fprintf(stderr,"Total Mass %g\nCOM position %g %g %g\nCOM velocity %g %g %g\nAngular Momentum  %g %g %g\n",mass,xcom[0]/mass,xcom[1]/mass,xcom[2]/mass,vcom[0]/mass,vcom[1]/mass,vcom[2]/mass,L[0]/mass,L[1]/mass,L[2]/mass);
        fprintf(stderr,"velocity dispersion %g %g %g  %g\n",sqrt(v2com[0]/mass-vcom[0]*vcom[0]/mass/mass),sqrt(v2com[1]/mass-vcom[1]*vcom[1]/mass/mass),sqrt(v2com[2]/mass-vcom[2]*vcom[2]/mass/mass),sqrt(v2com[0]/mass-vcom[0]*vcom[0]/mass/mass+v2com[1]/mass-vcom[1]*vcom[1]/mass/mass+v2com[2]/mass-vcom[2]*vcom[2]/mass/mass));
        }

    if (auxExt != NULL) {
        char allExt[] = "ss,tipsy,iord,den,denRFC,col,pot,amag,rung,mass,dt,SPHdt,soft,denu,pres,temperature,GasDensity,u,MassHot,uHot,Tinc,uDot,HI,HeI,HeII,H2,correL,lw,radTens,radKappa,radKappa0,radKappa1,radKappa2,radKappa3,radKappa4,radKappa5,radKappa6,radKappa7,radKappa8,radKappa9,radFlux,radFlux0,radFlux1,radFlux2,radFlux3,radFlux4,radFlux5,radFlux6,radFlux7,radFlux8,radFlux9,radTau,radTau0,radTau1,radTau2,radTau3,radTau4,radTau5,radTau6,radTau7,radTau8,radTau9,radIType,eDot,eCool,eHeat,BSw,alpha,divv,dvxdx,dvydx,dvzdx,dvxdy,dvydy,dvzdy,dvxdz,dvydz,dvzdz,divv_dens,divvdot,dvds,vsigmax,rcd,snorm,sfull,dvdsonsfull,alphaloc,alphanoise,area,divvt,divvcorr,ia,c,mumax,dch,dcx,uDotHydro,PdVRFC,uDotPdV,uDotAV,uDotDiff,uHotDot,uHotDotDiff,uHotDotPdV,uHotDotConv,uHotDotFB,Metals,Metalsdot,ST,divrhov,sigma2,igasorder,timeform,massform,coolontime,OxMassFrac,FeMassFrac,OxMassFracdot,FeMassFracdot,uDotFB,tcooloff,tcoolyr,tdynyr,rsnddyn,ljeans,ijeans,tCoolAgain,mStar,SPHH,smoothlength,pos,vel,acc,accg,accRFC,curl,norm,vpred,gradrho,accp,angmom";
        char auxFile[256];
        char *start,*next;

        if (auxExt[0]=='\0' || !strcmp(auxExt,"all") || !strcmp(auxExt,".")) {
            auxExt = allExt;
            }

        start = auxExt;
        for (;;) {
            if (start[0]=='\0') break;
            for (next=start;next[0]!='\0';next++)
                if (next[0]==',') { next[0]='\0'; next++; break;}
            // printf("Trying %s.%s\n",file, start);
            snprintf(auxFile, 256, "%s.%s", file, start);
            PrintTipsy( auxFile, bDark,bGas,bStar, nth, istart, iend, bCheck, bProp, 1, iAuxFormat, bSphere );
            start = next;
            }
        }
    
    return 0;
    }

