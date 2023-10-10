/* 
   tipsysussheader
   Utility for overwriting tipsy headers (native 32 byte headers only)
   Supports 40 bit extended headers (pkdgrav3, up to 1 T particles)
   (c) 2023 James Wadsley

   Usage: tipsysussheader npsh ndark nstar [file] [time]
   will use time from file if exists if not specified
   will overwrite file header if file exists or creates just header in file otherwise

   gcc -I/usr/include/tirpc tipsysussheader.c -ltirpc -lm -o tipsysussheader
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdint.h>

/* pkdgrav3 ( pad bytes required as well ) */
struct dump {
    double time ;
    unsigned int nbodies ;
    unsigned int ndim ;
    unsigned int nsph ;
    unsigned int ndark ;
    unsigned int nstar ;
    unsigned int npad ;
} ;

int main(int argc, char **argv) {
    FILE *fpread, *fpwrite;
    uint64_t pkd3N, pkd3Dark, pkd3Sph, pkd3Star;
    double time=0.;
    char *file="(none)";
    int bNewFile=0;
    struct dump hwrite;

    if (argc<4) {
        printf("Usage: tipsysussheader nsph ndark nstar [file] [time]\n");
        exit(1);
        }
    
    pkd3Sph=atol(argv[1]);
    pkd3Dark=atol(argv[2]);
    pkd3Star=atol(argv[3]);
    pkd3N=pkd3Sph+pkd3Dark+pkd3Star;

    hwrite.nbodies=pkd3N & 0xffffffff;
    hwrite.ndim   =3;
    hwrite.nsph   =pkd3Sph & 0xffffffff;
    hwrite.ndark  =pkd3Dark & 0xffffffff;
    hwrite.nstar  =pkd3Star & 0xffffffff;
    hwrite.npad   =(pkd3N >> 32) + ((pkd3Sph >> 24)& 0x0000ff00) +
        ((pkd3Dark >> 16) & 0x00ff0000) + ((pkd3Star >> 8) & 0xff000000);

    if (argc>=5) {
    
        int   fildes,status;
        file=argv[4];
        
        fildes = open(file, O_RDONLY);
        if (fildes < 0) {
            bNewFile=1;
            fprintf(stderr,"Can't open file: %s\n",file); /* assume does not exist */
            time=0.;
            }
        else {
            int   fildes,status;
            struct stat fileinfo;
            long long size;
            struct dump hread;
            status = fstat(fildes, &fileinfo);
            size = fileinfo.st_size;
            close(fildes);
            fpread = fopen( file, "r" );
            if (fpread == NULL) {
                fprintf(stderr,"Can't open file: %s\n",file);
                exit(1);
                }
            fread(&hread,28,1,fpread);
            fclose(fpread);
            
            time = hread.time;
            fprintf(stderr,"READ: %s  time=%lf n=%u ndim=%u sph=%u dark=%u star=%u pad=%u\n",file,hread.time,hread.nbodies,hread.ndim,hread.nsph,hread.ndark,hread.nstar,hread.npad);
            }
        
        if (argc==6) {
            time=atof(argv[5]);
            fprintf(stderr,"command line  time=%lf\n",time);
            }
        }
    hwrite.time = time;

    /* all header info collected */
        
    fprintf(stderr,"WRITE: %s  time=%lf n=%u ndim=%u sph=%u dark=%u star=%u pad=%u\n",file,hwrite.time,hwrite.nbodies,hwrite.ndim,hwrite.nsph,hwrite.ndark,hwrite.nstar,hwrite.npad);

    if (argc>=5) {
        /* Write to file */
        if (bNewFile) {
            fpwrite = fopen( file, "w" );
            }
        else {
            fpwrite = fopen( file, "r+" );
            /* Overwrite header only -- filepointer = 0 as needed */
            }
        if (fpwrite==NULL) {
            fprintf(stderr,"Error could not open %s for write\n",file);
            exit(2);
            }

        fwrite(&hwrite,sizeof(struct dump),1,fpwrite);
        fclose(fpwrite);
        }
    }
