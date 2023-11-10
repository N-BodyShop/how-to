/* 
   printtipsy
   Utility for examining contents of tipsy binary files
   (c) 2012-2023 James Wadsley

   Usage: arrayascii2bin filename

   Make:  gcc arrayascii2bin.c -lm -o arrayascii2bin
   Note: May need tirpc to get xdr functions:
   Make2: gcc -I/usr/include/tirpc printtipsy.c -ltirpc -lm -o printtipsy
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

typedef float Real;


#define Real float

int nSph=-1,nDark=-1,nStar=-1,bReadAux=0;

//  ReadArray( file, filew, bDark,bGas,bStar, bReplace, iWriteFormat, bInt, bCreate, iReplace, fReplace);
int ReadArray(char *file, char *filew, int bDark,int bGas,int bStar,int bReplace,int iWriteFormat,int bInt, int bCreate,int iReplace, float fReplace)
    {
    int i,iReadFormat=0;
    long offset;
    XDR xdrread,xdrwrite;
    FILE *fpread,*fpwrite;
    struct stat fileinfo;
    int   fildes,status;
    long long size;
    int npart;

    fildes = open(file, O_RDONLY);
    if (fildes < 0) {
        fprintf(stderr,"Can't open file: %s\n",file);
        exit(1);
        }
    status = fstat(fildes, &fileinfo);
    size = fileinfo.st_size;
    close(fildes);
	
    fpread = fopen( file, "r" );
    if (fpread == NULL) {
        fprintf(stderr,"Can't open file: %s\n",file);
        exit(1);
        }

    fread(&npart,4,1,fpread);
    if (npart == (size-4)/4) {
        fprintf(stderr,"%d Particles in Native Binary array file %s\n",npart,file);
        iReadFormat=0;
        }
    else {
        rewind(fpread);
        xdrstdio_create(&xdrread, fpread, XDR_DECODE);
        xdr_int( &xdrread, &npart );
        if (npart == (size-4)/4) {
            fprintf(stderr,"%d Particles in STD Binary array file %s\n",npart,file);
            iReadFormat=1;
            }
        else {
            const int nChar=20;
            if (size<nChar) {
                fprintf(stderr,"File is very short -- truncated? Aborting.\n");
                exit(2);
            }
            rewind(fpread);
            char readtest[nChar+1];
            fscanf(fpread, "%20s", readtest);
            int bNumber=0,nOK=0;
            // Check for whitespace or numbers in 1st chars
            for (i=0;i<nChar;i++) {
                switch(readtest[i]) {
                case 48: // 0
                case 49:
                case 50:
                case 51:
                case 52:
                case 53:
                case 54:
                case 55:
                case 56:
                case 57: // 9
                    bNumber=1;
                case 9: // tab
                case 10: // linefeed
                case 13: // return
                case 32: // space
                case 43: // +
                case 44: // ,
                case 45: // -
                case 46: // .
                case 101: // e
                case 69: // E
                    nOK++;
                    }
                }
            if (nOK==nChar && bNumber) { // ascii text numbers
                rewind(fpread);
                fscanf(fpread, "%d", &npart);
                fprintf(stderr,"%d Particles in ascii text array file %s\n",npart,file);
                iReadFormat=-1;
                }
            else {
                fprintf(stderr,"Unknown format. Aborting\n");
                exit(1);
                }
            }
        }

    fpwrite = fopen( filew, "w" );
    if (fpwrite == NULL) {
        fprintf(stderr,"Can't open file for writing: %s\n",filew);
        exit(1);
        }

    // Write new headers
    switch(iWriteFormat) {
    case 0:
        fwrite(&npart,4,1,fpwrite);
        break;
    case 1:
        xdrstdio_create(&xdrwrite, fpwrite, XDR_ENCODE);
        xdr_int( &xdrwrite, &npart );
        break;
    case -1:
        fprintf(fpwrite,"%d\n",npart);
        break;
        }

    int DataSize=sizeof(float);
    int iRead=0,nChunk0=1024,nChunk;
    nChunk=nChunk0;
    union { float f; int i; } Data[nChunk0];

    if (bReplace) {
        iWriteFormat=-2;
        if (bInt) {
            for (i=0;i<nChunk;i++) {
                Data[i].i = iReplace;
                }
            }
        else { 
            for (i=0;i<nChunk;i++) {
                Data[i].f = fReplace;
                }
            }
        }
           
        
    for (;;) {
        if (iRead+nChunk0 > npart) {
            nChunk=npart-iRead;
            }
        switch(iReadFormat) {
        case 0:
            fread(&Data[0],DataSize*nChunk,1,fpread);
            break;
        case 1:
            if (bInt) {
                for (i=0;i<nChunk;i++) {
                    xdr_int(&xdrread,&(Data[i].i));
                    }
                }
            else {
                for (i=0;i<nChunk;i++) {
                    xdr_float(&xdrread,&(Data[i].f));
                    }
                }
            break;
        case -1:
            if (bInt) {
                for (i=0;i<nChunk;i++) {
                    fscanf(fpread,"%d",&(Data[i].i));
                    }
                }
            else {
                for (i=0;i<nChunk;i++) {
                    fscanf(fpread,"%f",&(Data[i].f));
                    }
                }
        case -2:
            break;
        default:
            exit(1);
            }
        switch(iWriteFormat) {
        case 0:
            fwrite(&Data[0],DataSize*nChunk,1,fpwrite);
            break;
        case 1:
            if (bInt) {
                for (i=0;i<nChunk;i++) {
                    xdr_int(&xdrwrite,&(Data[i].i));
                    }
                }
            else {
                for (i=0;i<nChunk;i++) {
                    xdr_float(&xdrwrite,&(Data[i].f));
                    }
                }
            break;
        case -1:
            if (bInt) {
                for (i=0;i<nChunk;i++) {
                    fprintf(fpwrite,"%d",(Data[i].i));
                    }
                }
            else {
                for (i=0;i<nChunk;i++) {
                    fprintf(fpwrite,"%f",(Data[i].f));
                    }
                }
            break;
        default:
            exit(1);
            }

        iRead+=nChunk;
        if (iRead >= npart) break;
        }

    fclose(fpread);
	
    return npart;
    }

void usage() {
    fprintf(stderr,"Usage: arrayascii2binary [-gds]\n");
    exit(1);
    }

int main(int argc, char **argv) {
    int i,bGas=0,bDark=0,bStar=0;
    char defaultw[]="arrayout", *filew=defaultw;
    char *pch, *file=NULL;
    float fReplace;
    int iReplace;
    int bReplace=0,bCreate=0,iWriteFormat=0,bInt=0;
    
    i=1;
    while(i < argc) {
        pch = argv[i];
        if (!strcmp(argv[i],"-r")) {
            ++i;
            bReplace=1;
            if (i >= argc) usage();
            if (bInt)
                iReplace = atoi(argv[i]);
            else 
                fReplace = atof(argv[i]);
            ++i;
            }
        if (!strcmp(argv[i],"-w")) {
            ++i;
            if (i >= argc) usage();
            filew = argv[i];
            ++i;
            }
        else if (!strcmp(argv[i],"-std")) {
            iWriteFormat=1; // xdr
            ++i;
            }
        else if (!strcmp(argv[i],"-a")) {
            iWriteFormat=-1; // ascii
            ++i;
            }
        else if (!strcmp(argv[i],"-i")) {
            bInt=1;
            ++i;
            if (bReplace) {
                fprintf(stderr,"set int mode before setting replace value\n");
                exit(1);
                }
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

    ReadArray( file, filew, bDark,bGas,bStar, bReplace, iWriteFormat, bInt, bCreate, iReplace, fReplace);

    return 0;
    }

