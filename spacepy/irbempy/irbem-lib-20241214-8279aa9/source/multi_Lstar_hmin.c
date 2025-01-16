/*
An MPI program to genereate multiple Lstar/hmin values
From time/loc parameters from a binary file

Example compile command:
mpicc -s -O3 -L/home/tpo26062/lib -o multi_Lstar_hmin.mpix multi_Lstar_hmin.c -lonera_desp_lib_linux_x86_64
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <mpi.h>

#define MY_MPI_OP_EXIT (0)
#define MY_MPI_OP_GO (1)
#define MY_MPI_TAG_OP       (0)
#define MY_MPI_TAG_COMMON   (1)
#define MY_MPI_TAG_VARINPUT (2)
#define MY_MPI_TAG_OUTPUT   (3)

/* from onera_desp_lib.h */
void drift_bounce_orbit2_1_(int *kext, int *options,
			    int *sysaxes, int *iyear,
			    int *idoy, double * UT,
			    double *x1, double *x2, double *x3, double *alpha,
			    double *maginput, double *R0,
			    double *Lm,
			    double *Lstar, 
			    double *Blocal,
			    double *Bmin, 
			    double *Bmir, 
			    double *J,
			    double *posit,
			    int *ind, double *hmin, double *hmin_lon);

void mpi_cleanup() {
  int mpi_slave, mpi_size=1, mpi_op_code;

  MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);	/* get number of processes */  
  printf("%s: cleaning up %d nodes\n",__func__,mpi_size);

  if (mpi_size>1) {
    /* send op code exit to slaves, i.e., slaves die */
    mpi_op_code = MY_MPI_OP_EXIT;
    for (mpi_slave=1; mpi_slave < mpi_size; mpi_slave++) {
      MPI_Send(&mpi_op_code,1,MPI_INT,mpi_slave,MY_MPI_TAG_OP,MPI_COMM_WORLD);
    }
  }
  
  printf("%s: finalizing\n",__func__);
  MPI_Finalize();  
}

void cleanup(int exit_code){
  /* cleans up, including sending end op code 0 to mpi children */
  /* if exit_code !=0, exits */
  mpi_cleanup(); /* kill slaves as needed */
  if(exit_code) {
    exit(exit_code);
  }
}

int mpi_init(int *argc,char **argv[]) {
  int mpi_rank=0;
  MPI_Init (argc, argv);	/* starts MPI */  
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);	/* get current process id */  
  printf("%s: MPI initialized, rank=%d\n",__func__,mpi_rank);
  return(mpi_rank);
}

void print_usage() {
  printf("multi_Lstar_hmin.mpix infile outfile\n");
}


#define FILENAME_LEN (2048)
#define HEADER_LEN (100)
#define HEADER_MAGIC ("multi_Lstar_hmin")
int main(int argc, char *argv[]) {

  /* internal variables */
  char infilename[FILENAME_LEN],outfilename[FILENAME_LEN];
  char header[HEADER_LEN];
  int maginput_varies=0,maginput_size=0;
  long int first_step,istep;
  FILE *infilep,*outfilep;
  int32_t ntimes,nchunk;
  int mpi_master, mpi_slave, mpi_size=1, mpi_op_code,mpi_rank=0;
  MPI_Status status;

  /* time variable arguments to drift_bounce_orbit_2_1 */
  int32_t *iyear,*idoy; /* inputs */
  double *UT,*x1,*x2,*x3,*alpha, *maginput; /* inputs */
  double *Lm,*Lstar,*Bmin,*Bmir,*J,*hmin,*hmin_lon; /* outputs */
  /* note: returned values Blocal, posit, and ind are ignored */

  /* time-fixed arguments to drift_bounce_orbit2_1 */
  int32_t kext, options[5], sysaxes;
  double R0,Blocal[1000*25],posit[3*1000*25];
  int32_t ind[25];


  mpi_master = (mpi_init(&argc,&argv) == 0); /* initialize mpi */

  if (mpi_master) {
    if (sizeof(int32_t) != sizeof(int)) {
      perror("main: int is not 32-bit");
      cleanup(-1);
    }
    /* read arguments */
    printf("%s: reading %d args\n",__func__,argc); /* debug */
    if (argc != 3) {
      fprintf(stderr,"%s: wrong syntax\n",__func__);
      print_usage();
      cleanup(-1);
    }

    if ((strlen(argv[1])>=FILENAME_LEN) || (strlen(argv[2])>=FILENAME_LEN)) {
      fprintf(stderr,"%s: maximum file name length is %i\n",__func__,FILENAME_LEN-1);
      cleanup(-1);
    }

    strcpy(infilename,argv[1]);
    strcpy(outfilename,argv[2]);
    printf("Master: infilename = <%s>, outfilename = <%s>\n",infilename,outfilename);
    
    MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);	/* get number of processes */  
    printf("Master: %i processes\n",mpi_size);

    if (!(infilep = fopen(infilename,"rb"))) {
      perror("main:Unable to create file");
      fprintf(stderr,"%s: file name: %s\n",__func__,infilename);
      cleanup(-1);
    }

    /* check header */
    fread(header, 1, HEADER_LEN, infilep); /* header */
    if (strncmp(header,HEADER_MAGIC,strlen(HEADER_MAGIC))!=0) {
      fclose(infilep);
      perror("main:Header does not start with magic string");
      fprintf(stderr,"%s: magic=<%s> file name: %s\n",__func__,HEADER_MAGIC,infilename);
      cleanup(-1);
    }
    maginput_varies = header[strlen(HEADER_MAGIC)] & 1;
    printf("Master: maginput_varies=%i\n",maginput_varies);
    
    /* load data, allocate large arrays */
    
    fread(&ntimes, sizeof(ntimes), 1, infilep); /* number of samples */
    printf("Master: ntimes=%i\n",ntimes);

    /* send GO code to slaves - MPI */
    mpi_op_code = MY_MPI_OP_GO;
    nchunk = ntimes / mpi_size; /* number of points per slave */
    for (mpi_slave = 1; mpi_slave < mpi_size; mpi_slave++) {
      printf("Master: Sending GO code to slave %i\n",mpi_slave);
      MPI_Send(&mpi_op_code,1,MPI_INT,mpi_slave,MY_MPI_TAG_OP,MPI_COMM_WORLD);
      printf("Master: Sending maginput_varies=%i to slave %i\n",maginput_varies,mpi_slave);
      MPI_Send(&maginput_varies,1,MPI_INT,mpi_slave,MY_MPI_TAG_COMMON,MPI_COMM_WORLD);
      printf("Master: Sending nchunk=%i to slave %i\n",nchunk,mpi_slave);
      MPI_Send(&nchunk,1,MPI_INT,mpi_slave,MY_MPI_TAG_COMMON,MPI_COMM_WORLD);
    }


  } else {
    MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);	/* get slave process id */  

    /* wait for GO op code or abort */
    printf("Slave %i: waiting for op code\n",mpi_rank);
    MPI_Recv(&mpi_op_code,1,MPI_INT,0,MY_MPI_TAG_OP,MPI_COMM_WORLD,&status);
    printf("Slave %i: received op code %i\n",mpi_rank,mpi_op_code);
    if (mpi_op_code != MY_MPI_OP_GO) { /* early abort */
      MPI_Finalize();
      exit(-1);
    }

    /* get maginput_varies via MPI */
    MPI_Recv(&maginput_varies,1,MPI_INT,0,MY_MPI_TAG_COMMON,MPI_COMM_WORLD,&status);
    printf("Slave %i: received maginput_varies %i\n",mpi_rank,maginput_varies);
    /* get ntimes via MPI */
    MPI_Recv(&ntimes,1,MPI_INT,0,MY_MPI_TAG_COMMON,MPI_COMM_WORLD,&status);
    printf("Slave %i: received ntimes %i\n",mpi_rank,ntimes);

  }
  
  maginput_size = (maginput_varies?ntimes:1)*25;
    
  /* allocate time series of parameters */
  iyear = malloc(ntimes*sizeof(*iyear));
  idoy = malloc(ntimes*sizeof(*idoy));
  UT = malloc(ntimes*sizeof(*UT));
  x1 = malloc(ntimes*sizeof(*x1));
  x2 = malloc(ntimes*sizeof(*x2));
  x3 = malloc(ntimes*sizeof(*x3));
  alpha = malloc(ntimes*sizeof(*alpha));
  maginput = malloc(maginput_size*sizeof(*maginput));

  /* allocate output variables */
  Lm = malloc(ntimes*sizeof(*Lm));
  Lstar = malloc(ntimes*sizeof(*Lstar));
  Bmin = malloc(ntimes*sizeof(*Bmin));
  Bmir = malloc(ntimes*sizeof(*Bmir));
  J = malloc(ntimes*sizeof(*J));
  hmin = malloc(ntimes*sizeof(*hmin));
  hmin_lon = malloc(ntimes*sizeof(*hmin_lon));

  
  if (mpi_master) {
    /* read common parameters: kext, options, sysaxes, R0 */
    fread(&kext, sizeof(kext), 1, infilep);
    printf("Master: kext=%i\n",kext);
    fread(options, sizeof(*options), 5, infilep);
    fread(&sysaxes, sizeof(sysaxes), 1, infilep);
    printf("Master: sysaxes=%i\n",sysaxes);
    fread(&R0, sizeof(R0), 1, infilep);
    printf("Master: R0=%g\n",R0);
    
    /* read time series params */
    printf("Master: Reading VARINPUT data\n");
    fread(maginput, sizeof(*maginput), maginput_size, infilep);
    fread(iyear, sizeof(*iyear), ntimes, infilep);
    fread(idoy, sizeof(*idoy), ntimes, infilep);
    fread(UT, sizeof(*UT), ntimes, infilep);
    fread(x1, sizeof(*x1), ntimes, infilep);
    fread(x2, sizeof(*x2), ntimes, infilep);
    fread(x3, sizeof(*x3), ntimes, infilep);
    fread(alpha, sizeof(*alpha), ntimes, infilep);
    fclose(infilep);
    printf("Master: Done reading from %s\n",infilename);
    
    /* prep out file */
    
    if (!(outfilep = fopen(outfilename,"wb"))) {
      perror("main:Unable to create file");
      fprintf(stderr,"%s: file name: %s\n",__func__,outfilename);
      cleanup(-1);
    }

    printf("Master: Dispatching to %i slaves\n",mpi_size-1);
    /* send GO op code, break up job & send to slaves */
    mpi_op_code = MY_MPI_OP_GO;
    first_step = 0;
    for (mpi_slave = 1; mpi_slave < mpi_size; mpi_slave++) {
      /* send GO code to slaves - MPI */
      printf("Master: Sending GO code to slave %i\n",mpi_slave);
      MPI_Send(&mpi_op_code,1,MPI_INT,mpi_slave,MY_MPI_TAG_OP,MPI_COMM_WORLD);

      /* send COMMON data */
      printf("Master: Sending COMMON data to slave %i\n",mpi_slave);
      MPI_Send(&kext,1,MPI_INT,mpi_slave,MY_MPI_TAG_COMMON,MPI_COMM_WORLD);
      MPI_Send(options,5,MPI_INT,mpi_slave,MY_MPI_TAG_COMMON,MPI_COMM_WORLD);
      MPI_Send(&sysaxes,1,MPI_INT,mpi_slave,MY_MPI_TAG_COMMON,MPI_COMM_WORLD);
      MPI_Send(&R0,1,MPI_DOUBLE,mpi_slave,MY_MPI_TAG_COMMON,MPI_COMM_WORLD);

      /* send VARINPUT data */
      printf("Master: Sending VARINPUT data to slave %i, first_step=%i, nchunk=%i\n",mpi_slave,first_step,nchunk);
      if (maginput_varies) {
	MPI_Send(maginput+first_step*25,nchunk*25,MPI_DOUBLE,mpi_slave,MY_MPI_TAG_VARINPUT,MPI_COMM_WORLD);
      } else {
	MPI_Send(maginput,25,MPI_DOUBLE,mpi_slave,MY_MPI_TAG_VARINPUT,MPI_COMM_WORLD);
      }
      MPI_Send(iyear+first_step,nchunk,MPI_INT,mpi_slave,MY_MPI_TAG_VARINPUT,MPI_COMM_WORLD);
      MPI_Send(idoy+first_step,nchunk,MPI_INT,mpi_slave,MY_MPI_TAG_VARINPUT,MPI_COMM_WORLD);
      MPI_Send(UT+first_step,nchunk,MPI_DOUBLE,mpi_slave,MY_MPI_TAG_VARINPUT,MPI_COMM_WORLD);
      MPI_Send(x1+first_step,nchunk,MPI_DOUBLE,mpi_slave,MY_MPI_TAG_VARINPUT,MPI_COMM_WORLD);
      MPI_Send(x2+first_step,nchunk,MPI_DOUBLE,mpi_slave,MY_MPI_TAG_VARINPUT,MPI_COMM_WORLD);
      MPI_Send(x3+first_step,nchunk,MPI_DOUBLE,mpi_slave,MY_MPI_TAG_VARINPUT,MPI_COMM_WORLD);
      MPI_Send(alpha+first_step,nchunk,MPI_DOUBLE,mpi_slave,MY_MPI_TAG_VARINPUT,MPI_COMM_WORLD);

      first_step += nchunk;
      printf("Master: Done sending VARINPUT data to slave %i, new first_step=%i\n",mpi_slave,first_step);
    }
    printf("Master: Keeping %i for myself (vs %i for slaves)\n",ntimes-first_step,nchunk);

  } else {
    /* wait for GO op code or abort */
    printf("Slave %i: waiting for op code\n",mpi_rank);
    MPI_Recv(&mpi_op_code,1,MPI_INT,0,MY_MPI_TAG_OP,MPI_COMM_WORLD,&status);
    printf("Slave %i: received op code %i\n",mpi_rank,mpi_op_code);
    if (mpi_op_code != MY_MPI_OP_GO) { /* early abort */
      MPI_Finalize();
      exit(-1);
    }

    /* receive common params */
    MPI_Recv(&kext,1,MPI_INT,0,MY_MPI_TAG_COMMON,MPI_COMM_WORLD,&status);
    printf("Slave %i: received kext=%i\n",mpi_rank,kext);
    MPI_Recv(options,5,MPI_INT,0,MY_MPI_TAG_COMMON,MPI_COMM_WORLD,&status);
    MPI_Recv(&sysaxes,1,MPI_INT,0,MY_MPI_TAG_COMMON,MPI_COMM_WORLD,&status);
    printf("Slave %i: received sysaxes=%i\n",mpi_rank,sysaxes);
    MPI_Recv(&R0,1,MPI_DOUBLE,0,MY_MPI_TAG_COMMON,MPI_COMM_WORLD,&status);
    printf("Slave %i: received R0=%g\n",mpi_rank,R0);

    /* receive time series params */
    MPI_Recv(maginput,maginput_size,MPI_DOUBLE,0,MY_MPI_TAG_VARINPUT,MPI_COMM_WORLD,&status);
    MPI_Recv(iyear,ntimes,MPI_INT,0,MY_MPI_TAG_VARINPUT,MPI_COMM_WORLD,&status);
    MPI_Recv(idoy,ntimes,MPI_INT,0,MY_MPI_TAG_VARINPUT,MPI_COMM_WORLD,&status);
    MPI_Recv(UT,ntimes,MPI_DOUBLE,0,MY_MPI_TAG_VARINPUT,MPI_COMM_WORLD,&status);
    MPI_Recv(x1,ntimes,MPI_DOUBLE,0,MY_MPI_TAG_VARINPUT,MPI_COMM_WORLD,&status);
    MPI_Recv(x2,ntimes,MPI_DOUBLE,0,MY_MPI_TAG_VARINPUT,MPI_COMM_WORLD,&status);
    MPI_Recv(x3,ntimes,MPI_DOUBLE,0,MY_MPI_TAG_VARINPUT,MPI_COMM_WORLD,&status);
    MPI_Recv(alpha,ntimes,MPI_DOUBLE,0,MY_MPI_TAG_VARINPUT,MPI_COMM_WORLD,&status);

    printf("Slave %i: done receiving run parameters\n",mpi_rank);

    first_step = 0; /* slaves always to 0 to ntimes-1 */
  }
    
  /* run loop */

  printf("Master/Slave %i: running %i cases, %i to %i\n",mpi_rank,ntimes-first_step,first_step,ntimes-1);

  for (istep = first_step; istep < ntimes; istep++) {
    printf("Master/Slave %i: Running case %i / %i\n",mpi_rank,istep,ntimes);
    drift_bounce_orbit2_1_(&kext,options,&sysaxes, iyear+istep,idoy+istep,UT+istep,x1+istep,x2+istep,x3+istep,alpha+istep,
			   maginput+istep*25*maginput_varies, 
			   &R0,Lm+istep,Lstar+istep,Blocal,Bmin+istep,Bmir+istep,J+istep,posit,ind,hmin+istep,hmin_lon+istep);
  }

  printf("Master/Slave %i: Done running %i cases\n",mpi_rank,ntimes-first_step);

  /* collect */
  if (mpi_master) {
    /* receive results from slaves via MPI */

    first_step = 0;
    for (mpi_slave = 1; mpi_slave < mpi_size; mpi_slave++) {
      printf("Master: receiving %i results from %i, first_step = %i\n",nchunk,mpi_slave,first_step);
      MPI_Recv(Lm+first_step,nchunk,MPI_DOUBLE,mpi_slave,MY_MPI_TAG_OUTPUT,MPI_COMM_WORLD,&status);
      MPI_Recv(Lstar+first_step,nchunk,MPI_DOUBLE,mpi_slave,MY_MPI_TAG_OUTPUT,MPI_COMM_WORLD,&status);
      MPI_Recv(Bmin+first_step,nchunk,MPI_DOUBLE,mpi_slave,MY_MPI_TAG_OUTPUT,MPI_COMM_WORLD,&status);
      MPI_Recv(Bmir+first_step,nchunk,MPI_DOUBLE,mpi_slave,MY_MPI_TAG_OUTPUT,MPI_COMM_WORLD,&status);
      MPI_Recv(J+first_step,nchunk,MPI_DOUBLE,mpi_slave,MY_MPI_TAG_OUTPUT,MPI_COMM_WORLD,&status);
      MPI_Recv(hmin+first_step,nchunk,MPI_DOUBLE,mpi_slave,MY_MPI_TAG_OUTPUT,MPI_COMM_WORLD,&status);
      MPI_Recv(hmin_lon+first_step,nchunk,MPI_DOUBLE,mpi_slave,MY_MPI_TAG_OUTPUT,MPI_COMM_WORLD,&status);
      first_step += nchunk;
      printf("Master: received %i results from %i, new first_step=%i\n",nchunk,mpi_slave,first_step);
    }

    /* write */
    printf("Master: writing data to%s\n",outfilename);
    fwrite(&ntimes,sizeof(ntimes),1,outfilep);
    fwrite(Lm,sizeof(*Lm),ntimes,outfilep);
    fwrite(Lstar,sizeof(*Lm),ntimes,outfilep);
    fwrite(Bmin,sizeof(*Lm),ntimes,outfilep);
    fwrite(Bmir,sizeof(*Lm),ntimes,outfilep);
    fwrite(J,sizeof(*Lm),ntimes,outfilep);
    fwrite(hmin,sizeof(*Lm),ntimes,outfilep);
    fwrite(hmin_lon,sizeof(*Lm),ntimes,outfilep);

    fclose(outfilep);
    printf("Master: Done\n");
  } else {
    /* send results via MPI, OUTPUT channel */
    printf("Slave %i: sending results\n",mpi_rank);
    MPI_Send(Lm,ntimes,MPI_DOUBLE,0,MY_MPI_TAG_OUTPUT,MPI_COMM_WORLD);
    MPI_Send(Lstar,ntimes,MPI_DOUBLE,0,MY_MPI_TAG_OUTPUT,MPI_COMM_WORLD);
    MPI_Send(Bmin,ntimes,MPI_DOUBLE,0,MY_MPI_TAG_OUTPUT,MPI_COMM_WORLD);
    MPI_Send(Bmir,ntimes,MPI_DOUBLE,0,MY_MPI_TAG_OUTPUT,MPI_COMM_WORLD);
    MPI_Send(J,ntimes,MPI_DOUBLE,0,MY_MPI_TAG_OUTPUT,MPI_COMM_WORLD);
    MPI_Send(hmin,ntimes,MPI_DOUBLE,0,MY_MPI_TAG_OUTPUT,MPI_COMM_WORLD);
    MPI_Send(hmin_lon,ntimes,MPI_DOUBLE,0,MY_MPI_TAG_OUTPUT,MPI_COMM_WORLD);
    printf("Slave %i: Done\n",mpi_rank);
  }


  /* free */
  free(iyear);
  free(idoy);
  free(UT);
  free(x1);
  free(x2);
  free(x3);
  free(alpha);
  free(maginput);
  free(Lm);
  free(Lstar);
  free(Bmin);
  free(Bmir);
  free(J);
  free(hmin);
  free(hmin_lon);

  /* finalize */
  MPI_Finalize();
  return(0);
}
