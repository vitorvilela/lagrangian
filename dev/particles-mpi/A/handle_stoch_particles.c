#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "uthash.h"
#include "lag_map.h"
#include "handle_stoch_particles.h"
#include "mpi.h"




notionalParticle_t *lagrangian = NULL;

idStack_t *idstack;


idStack_t* createStack(unsigned int MAX_NP_PROC) {
  
  idStack_t *S;
  S = malloc(MAX_NP_PROC*sizeof(unsigned int));
  S->size = 0;
  S->capacity = MAX_NP_PROC;
  
  return S;
  
}


void pop(idStack_t *S) {
  
  if(S->size == 0)
    printf("Stack is Empty");

  else {    
    S->size--;
    return S->id[S->size+1];    
  }
  
}

unsigned int top(idStack_t *S) {
  
  if(S->size == 0)
    printf("Stack is Empty");

  else {    
    return S->id[S->size-1];    
  }
  
}


void push(idStack_t *S, unsigned int id) {
  
  if(S->size == S->capacity) 
    printf("Stack is Full");
  
  else 
    S->id[S->size++] = id;
   
}


void U_INITIALIZE_PARTICLES(const int proc_id, const double a1, const double a2, const double a3, const double b1, const double b2, const double b3, const double A[]) {
  

  int MAX_NP_PROC = 10;
  
  long int idum3 = -3213853;
  long int *ptr_idum3 = &idum3;

  int p, sid;
    
  //Initialize idstack **********************
     
  idstack = createStack(MAX_NP_PROC);
  for(sid=MAX_NP_PROC; sid>=1; sid--) {
    push(idstack, sid);
    printf("\n\npush(%d) - proc(%d)", sid, proc_id);
  }  
  //*****************************************


  
  volume_t *ptr_volume;
  
  id_np_t *ptr_id_p;
  
  notionalParticle_t *ptr_particle;
      
  const double lx = b1-a1;
  const double ly = b2-a2;
  const double lz = b3-a3;
  
  const double vol_proc = lx*ly*lz;
    
  const double vol_total = (A[1]-A[0])*(A[3]-A[2])*(A[5]-A[4]);

  const int np_proc = (int)(np_total*(vol_proc/vol_total));
  
  printf("\n\nvol_proc(%d) = %f \n\n", proc_id, vol_proc); 

  printf("\n\nnp_proc(%d) = %d de np_total (%d)\n\n", proc_id, np_proc, np_total); 
  
  for(p=1; p<=np_proc; p++) {
       
    printf("\n\np(%d) - proc(%d)", p, proc_id);
    
    ptr_particle = malloc(sizeof(notionalParticle_t));        
    ptr_particle->id = p;
    HASH_ADD(hh, lagrangian, id, sizeof(int), ptr_particle);
    
    //remove from idstack **************
    sid = top(idstack);
    pop(idstack);
    printf("\n\npop(%d) - proc(%d)", sid, proc_id);
    //**********************************
    
    ptr_particle->proc = proc_id;
    ptr_particle->new_proc = proc_id;
            
    ptr_particle->posixo = a1 + ran2(ptr_idum3)*lx;
    ptr_particle->posiyo = a2 + ran2(ptr_idum3)*ly;
    ptr_particle->posizo = a3 + ran2(ptr_idum3)*lz;
 
    ptr_particle->posix = ptr_particle->posixo;
    ptr_particle->posiy = ptr_particle->posiyo;
    ptr_particle->posiz = ptr_particle->posizo;
           
    ptr_volume = SEARCH_PARTICLE(A, ptr_particle->posixo, ptr_particle->posiyo, ptr_particle->posizo, npl);
        
    ptr_particle->index.igx = ptr_volume->index.igx;
    ptr_particle->index.igy = ptr_volume->index.igy;
    ptr_particle->index.igz = ptr_volume->index.igz;
       
    ptr_particle->new_index.igx = ptr_volume->index.igx;
    ptr_particle->new_index.igy = ptr_volume->index.igy;
    ptr_particle->new_index.igz = ptr_volume->index.igz;
       
    ptr_particle->treact = 0.0;
    ptr_particle->dwall = 0.0;
    ptr_particle->omegap = 0.0;
    ptr_particle->rhop = 0.0;
    ptr_particle->weight = 0.0; 
     
    ptr_id_p = malloc(sizeof(id_np_t)); 
    ptr_id_p->id = p;      
    HASH_ADD_INT(ptr_volume->id_np, id, ptr_id_p); 
     
  }
  
}


void PRINT_PARTICLES_VTK(const int *ct, const int *proc_id) {
    
  char string[40];
  
  int part_counter = 0;
  
  notionalParticle_t *ptr_particle, *tmp_particle;
  
  part_counter = HASH_COUNT(lagrangian);
  
  sprintf(string, "../output/particles_proc%i_ct%i.vtk", *proc_id, *ct);
    
  FILE *pvtk;
  pvtk = fopen(string, "w");
  fprintf(pvtk, "# vtk DataFile Version 2.0\nPoints\nASCII\nDATASET POLYDATA\n");
  fprintf(pvtk, "Points %i float\n", part_counter);
       
  HASH_ITER(hh, lagrangian, ptr_particle, tmp_particle) {
    fprintf(pvtk, "%1.2e %1.2e %1.2e \n", ptr_particle->posixo, ptr_particle->posiyo, ptr_particle->posizo);
  }
  
  fprintf(pvtk, "POINT_DATA %d\n", part_counter);
  fprintf(pvtk, "SCALARS zeta float\n");
  fprintf(pvtk, "LOOKUP_TABLE default\n");
    
  HASH_ITER(hh, lagrangian, ptr_particle, tmp_particle) {
    fprintf(pvtk, "%f \n", ptr_particle->zetap[0]);
  }
  
  fclose(pvtk);
     
}


void INITIALIZE_PARTICLES_ZETA(const int k) {
    
  int proc_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  
  notionalParticle_t *ptr_particle, *tmp_particle;
  
  int c;
  
  if(k == 0) {
    
    HASH_ITER(hh, lagrangian, ptr_particle, tmp_particle) { 
      for(c=0; c<n_Yi; c++) {
	ptr_particle->zetap[c] = 0.0;
	printf("\nP %d - Zeta %f \n", ptr_particle->id, ptr_particle->zetap[c]);
      }
    }
  }
  
  else if(k == 1) {
    
//     const double zeta1_mix = 0.0;
//     const double zeta2_mix = 1.0;
//     const double targ_mix = 0.5;
//     const double mixlgth = 0.1;
    
    HASH_ITER(hh, lagrangian, ptr_particle, tmp_particle) { 
      for(c=0; c<n_Yi; c++) {
	ptr_particle->zetap[c] = 0.5*(zeta1_mix+zeta2_mix) + 0.5*(zeta1_mix-zeta2_mix) * tanh((2.0*(ptr_particle->posiyo))/mixlgth - 2.0*targ_mix/mixlgth);
	printf("\nP %d (%i) - Zeta %f \n", ptr_particle->id, proc_id, ptr_particle->zetap[c]);
      }
    }
  }
  
}


void TRANSPORT_STOCH_PARTICLES(const double A[], const double a1, const double a2, const double a3, const double b1, const double b2, const double b3, const int* ct_seed, const int proc_id) {
     
  char in[1];
  
  printf("\nTRANSPORT_STOCH_PARTICLES()\n");
  
  printf("\nct = %i\n", *ct_seed);
 
  int i, j;
  
  int transf_counter = 0;
  
  int id;
  
  int idum3 = -3213853;
   
  int *ptr_idum3;
    
  if(*ct_seed == 1)
    ptr_idum3 = &idum3;
  else
    *ptr_idum3 = ct_seed;
    
//   printf("\n*ptr_idum3  = %i\n", *ptr_idum3);

  MPI_Status status;
  MPI_Request R;
  
  MPI_Datatype MPI_Particle_type;
  MPI_Datatype type[ITENS] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
  MPI_Aint disp[ITENS], address[ITENS];
  
  int blocklen[ITENS] = {1, 1, 1, 1, n_Yi, 1};
  

    
  volume_t *ptr_volume, *ptr_old_vol;
   
  notionalParticle_t *iter_particle, *tmp_particle, *ptr_particle;
  
  bufferParticle_t *ptr_buffer;
  
  id_np_t *add_id, *rm_id;
    
  
  HASH_ITER(hh, lagrangian, iter_particle, tmp_particle) {
        
//     printf("\niter_particle->id = %i  -   id %i\n", iter_particle->id, *proc_id);
    
    // Transport - Convection + Weiner process
    
    double gasx = gasdev(ptr_idum3);
    double gasy = gasdev(ptr_idum3);
    double gasz = gasdev(ptr_idum3);
    double factor = sqrt((dt)*(2.0*Ni/Sc));
    
//     printf("\n\tsqrt(4) = %e\n", sqrt(4));
//     printf("\n\tfactor = %e\n", factor);
//     printf("\n\tgasx = %e\n", gasx);
//     printf("\n\tgasy = %e\n", gasy);
//     printf("\n\tgasz = %e\n", gasz);
//     
//     printf("\n\tU*(*dt) = %e\n", U*(*dt));
//     printf("\n\tV*(*dt) = %e\n", V*(*dt));
//     printf("\n\tW*(*dt) = %e\n", W*(*dt));
    
    iter_particle->posix = iter_particle->posixo + U*(dt) + factor*gasx;
    iter_particle->posiy = iter_particle->posiyo + V*(dt) + factor*gasy;      
    iter_particle->posiz = iter_particle->posizo + W*(dt) + factor*gasz;
	    
    

    
    // Wall Boundary Conditions
	
    if(iter_particle->posix < A[0])
      iter_particle->posix = A[0] + (A[0] - iter_particle->posix);
   
    if(iter_particle->posix > A[1]) 
      iter_particle->posix = A[1] - (iter_particle->posix - A[1]);
    
    if(iter_particle->posiy < A[2])
      iter_particle->posiy = A[2] + (A[2] - iter_particle->posiy);
	
    if(iter_particle->posiy > A[3])
      iter_particle->posiy = A[3] - (iter_particle->posiy - A[3]);
      
    if(iter_particle->posiz < A[4])
      iter_particle->posiz = A[4] + (A[4] - iter_particle->posiz);
	
    if(iter_particle->posiz > A[5])
      iter_particle->posiz = A[5] - (iter_particle->posiz - A[5]);     

    //Counter for Transport through process in X direction
    //Just two process, so the particle knows its new process  
          
    if(proc_id == 0) {
                
      if(iter_particle->posix > b1) {
	 printf("\nb1 %f\n", b1);
	printf("\nif(iter_particle->posix > *b1)   %f\n" , iter_particle->posix);
	transf_counter++;
	iter_particle->new_proc = 1;
	printf("\niter_particle->id = %i\n", iter_particle->id);
	printf("\ntransf_counter = %i\n", transf_counter);
      }
      
    }
         
//     printf("\n\titer_particle->posixo = %f\n", iter_particle->posixo);
//     printf("\n\titer_particle->posix = %f\n", iter_particle->posix);
//     
//     printf("\n\titer_particle->posiyo = %f\n", iter_particle->posiyo);
//     printf("\n\titer_particle->posiy = %f\n", iter_particle->posiy);
//     
//     printf("\n\titer_particle->posizo = %f\n", iter_particle->posizo);
//     printf("\n\titer_particle->posiz = %f\n", iter_particle->posiz);
    
    
    // Change particle's index if it remains in the same process, otherwise
    // it will be done after process transport
	
    if(iter_particle->proc == iter_particle->new_proc) {
      
      ptr_volume = SEARCH_PARTICLE(A, iter_particle->posix, iter_particle->posiy, iter_particle->posiz, npl);
	
      iter_particle->new_index.igx = ptr_volume->index.igx;
      iter_particle->new_index.igy = ptr_volume->index.igy;
      iter_particle->new_index.igz = ptr_volume->index.igz;
    
      printf("\n\t\titer_particle->index.igx = %i  -   id %i , proc_id %i\n", iter_particle->index.igx, iter_particle->id, proc_id);
      printf("\n\t\titer_particle->new_index.igx = %i  -   id %i , proc_id %i\n", iter_particle->new_index.igx, iter_particle->id, proc_id);
//     
      // Mapping 
    
      if( (iter_particle->new_index.igx != iter_particle->index.igx) || (iter_particle->new_index.igy != iter_particle->index.igy) || (iter_particle->new_index.igz != iter_particle->index.igz) ) {
            
	id = iter_particle->id;
	
// 	printf("\nIndex changed!, id %i , proc_id %i\n", id, *proc_id);
      
	add_id = malloc(sizeof(id_np_t));
// 	printf("\nadd_id = malloc(sizeof(id_np_t))  -  proc_id %i \n", *proc_id);
	add_id->id = id;
// 	printf("\nadd_id->id = id  -  proc_id %i \n", *proc_id);
	HASH_ADD_INT(ptr_volume->id_np, id, add_id);
// 	printf("\nHASH_ADD_INT(ptr_volume->id_np, id, add_id)  -  proc_id %i \n", *proc_id);
            
	ptr_old_vol = SEARCH_PARTICLE(A, iter_particle->posixo, iter_particle->posiyo, iter_particle->posizo, npl);
// 	printf("\nptr_old_vol = SEARCH_PARTICLE(A, iter_particle->posixo, iter_particle->posiyo, iter_particle->posizo, *npl)  -  proc_id %i \n", *proc_id);
	HASH_FIND_INT(ptr_old_vol->id_np, &id, rm_id);
// 	printf("\nHASH_FIND_INT(ptr_old_vol->id_np, id, rm_id) -  proc_id %i\n", *proc_id);
	HASH_DEL(ptr_old_vol->id_np, rm_id);
// 	printf("\nHASH_DEL(ptr_old_vol->id_np, rm_id)   -  proc_id %i\n", *proc_id);
	free(rm_id); //NECESSÁRIO?
// 	printf("\nfree(rm_id)   -  proc_id %i\n", *proc_id);
	rm_id = NULL; 
// 	printf("\nrm_id = NULL   -  proc_id %i\n", *proc_id);
     
      }
      
      
      iter_particle->posixo = iter_particle->posix;
      iter_particle->posiyo = iter_particle->posiy;
      iter_particle->posizo = iter_particle->posiz;
    
      iter_particle->index.igx = iter_particle->new_index.igx;
      iter_particle->index.igy = iter_particle->new_index.igy;
      iter_particle->index.igz = iter_particle->new_index.igz;
    
    }
  }
  
//   printf("\n\n\nOK HERE - proc_id %i\n\n\n", *proc_id);
  
       
  if(proc_id == 0) {
    
     printf("\ntransf_counter : %d, proc_id %d\n", transf_counter, proc_id);

     MPI_Send(&transf_counter, 1, MPI_INT, 1, FLAG1, MPI_COMM_WORLD);
    
     printf("\nMPI_Send(&transf_counter, 1, MPI_INT, 1, FLAG1, MPI_COMM_WORLD) -  proc_id %d\n", proc_id);
      
    if(transf_counter != 0) {
                      
      ptr_buffer = malloc(transf_counter*sizeof(bufferParticle_t));
     
      MPI_Address(&ptr_buffer[0].x, &address[0]);
      MPI_Address(&ptr_buffer[0].y, &address[1]);
      MPI_Address(&ptr_buffer[0].z, &address[2]);
      MPI_Address(&ptr_buffer[0].rhop, &address[3]);
      MPI_Address(&ptr_buffer[0].zetap, &address[4]);
      MPI_Address(&ptr_buffer[0].id, &address[5]);
    
      disp[0] = address[0] - address[0];
      disp[1] = address[1] - address[0];
      disp[2] = address[2] - address[0];
      disp[3] = address[3] - address[0];
      disp[4] = address[4] - address[0];
      disp[5] = address[5] - address[0];
          
      MPI_Type_create_struct(ITENS, blocklen, disp, type, &MPI_Particle_type);
      MPI_Type_commit(&MPI_Particle_type);  
          
      int tp = 0;
              
      HASH_ITER(hh, lagrangian, iter_particle, tmp_particle) {
  
	if(iter_particle->proc != iter_particle->new_proc) {
      
	  printf("\ntp : %d\n", tp);
	  ptr_buffer[tp++].x = iter_particle->posix;
	  ptr_buffer[tp++].y = iter_particle->posiy;
	  ptr_buffer[tp++].z = iter_particle->posiz;
	  ptr_buffer[tp++].rhop = iter_particle->rhop;
	  ptr_buffer[tp++].id = iter_particle->id;  
	  for(i=0; i<n_Yi; i++)
	    ptr_buffer[tp++].zetap[i] = iter_particle->zetap[i];
	    
	  id = iter_particle->id;  
	  
	  push(idstack, id);
	
	  //verificar existencia e disponibilidade de ptr_old_vol e rm_id 
	
	  ptr_old_vol = SEARCH_PARTICLE(A, iter_particle->posixo, iter_particle->posiyo, iter_particle->posizo, npl);
	  printf("\nptr_old_vol = SEARCH_PARTICLE(A, iter_particle->posixo, iter_particle->posiyo, iter_particle->posizo, *npl)\n");
	  HASH_FIND_INT(ptr_old_vol->id_np, &id, rm_id);
	  printf("\nHASH_FIND_INT(ptr_old_vol->id_np, id, rm_id)\n");
	  HASH_DEL(ptr_old_vol->id_np, rm_id);
	  printf("\nHASH_DEL(ptr_old_vol->id_np, rm_id)\n");
	  free(rm_id);
	  printf("\nfree(rm_id)\n");
	  rm_id = NULL; 
	  printf("\nrm_id = NULL\n");
	  
	  
	  //DELETAR TBM PARTICULA DE LAGRANGIAN
       
	  ptr_particle = iter_particle;
	  printf("\nptr_particle = iter_particle\n");
	  HASH_DEL(lagrangian, ptr_particle);
	  printf("\nHASH_DEL(lagrangian, ptr_particle)\n");
	  free(ptr_particle);
	  printf("\nfree(ptr_particle)\n");
	  ptr_particle = NULL;
	  printf("\nptr_particle = NULL\n");
	  
	}
      
      
	//Send buffer from proc 0 to proc 1
      
	MPI_Isend(&ptr_buffer, transf_counter, MPI_Particle_type, 1, FLAG2, MPI_COMM_WORLD, &R);
	MPI_Wait(&R, &status); 
            
	free(ptr_buffer);
	ptr_buffer = NULL;
  
      }
    }
            
  }
    
//   printf("\n\n\nSTOP HERE - proc_id %i\n\n\n", *proc_id);
 
    
  if(proc_id == 1) {
      
      MPI_Recv(&transf_counter, 1, MPI_INT, 0, FLAG1, MPI_COMM_WORLD, &status);
      printf("\nMPI_Recv(&transf_counter, 1, MPI_INT, 0, FLAG1, MPI_COMM_WORLD, &status)  -  proc_id %d\n", proc_id);

      printf("\ntransf_counter : %d, proc_id %d\n", transf_counter, proc_id);
           

      if(transf_counter != 0) {

	ptr_buffer = malloc(transf_counter*sizeof(bufferParticle_t));
	
	//Verificar se é possível e melhor criar tipo apenas uma vez para todos os processos em main()
	
	MPI_Address(&ptr_buffer[0].x, &address[0]);
	MPI_Address(&ptr_buffer[0].y, &address[1]);
	MPI_Address(&ptr_buffer[0].z, &address[2]);
	MPI_Address(&ptr_buffer[0].rhop, &address[3]);
	MPI_Address(&ptr_buffer[0].zetap, &address[4]);
	MPI_Address(&ptr_buffer[0].id, &address[5]);
    
	disp[0] = address[0] - address[0];
	disp[1] = address[1] - address[0];
	disp[2] = address[2] - address[0];
	disp[3] = address[3] - address[0];
	disp[4] = address[4] - address[0];
	disp[5] = address[5] - address[0];
	
	MPI_Type_create_struct(ITENS, blocklen, disp, type, &MPI_Particle_type);
	MPI_Type_commit(&MPI_Particle_type); 
	
	MPI_Recv(ptr_buffer, transf_counter, MPI_Particle_type, 0, FLAG2, MPI_COMM_WORLD, &status);
            
	for(i=0; i<transf_counter; i++) {
	
	  printf("\ni = \n", i);
	
	  ptr_volume = SEARCH_PARTICLE(A, ptr_buffer[i].x, ptr_buffer[i].y, ptr_buffer[i].z, npl);
	  add_id = malloc(sizeof(id_np_t));
	  printf("\nadd_id = malloc(sizeof(id_np_t))\n");
	
	  id = top(idstack);
	  add_id->id = id;
	  pop(idstack);
	  printf("\n\ntake id (%d) from stack - proc(%d)", id, proc_id);
	
	  HASH_ADD_INT(ptr_volume->id_np, id, add_id);
	  printf("\nHASH_ADD_INT(ptr_volume->id_np, id, add_id)\n");
	
	  ptr_particle = malloc(sizeof(notionalParticle_t));        
	  ptr_particle->id = id;
	  HASH_ADD(hh, lagrangian, id, sizeof(int), ptr_particle);
	
	  ptr_particle->proc = proc_id;
	  ptr_particle->new_proc = proc_id;
	
	  ptr_particle->posixo = ptr_buffer[i].x;
	  ptr_particle->posiyo = ptr_buffer[i].y;
	  ptr_particle->posizo = ptr_buffer[i].z;
	
	  ptr_particle->posix = ptr_particle->posixo;
	  ptr_particle->posiy = ptr_particle->posiyo;
	  ptr_particle->posiz = ptr_particle->posizo;
	
	  ptr_particle->index.igx = ptr_volume->index.igx;
	  ptr_particle->index.igy = ptr_volume->index.igy;
	  ptr_particle->index.igz = ptr_volume->index.igz;
       
	  ptr_particle->new_index.igx = ptr_volume->index.igx;
	  ptr_particle->new_index.igy = ptr_volume->index.igy;
	  ptr_particle->new_index.igz = ptr_volume->index.igz;
	
	  ptr_particle->treact = 0.0;
	  ptr_particle->dwall = 0.0;
	  ptr_particle->omegap = 0.0;
	  ptr_particle->weight = 0.0;
	
	  ptr_particle->rhop = ptr_buffer[i].rhop;
	  
	  for(j=0; j<n_Yi; j++)
	    ptr_particle->zetap[j] = ptr_buffer[i].zetap[j];
    
	}
      	
      
	free(ptr_buffer);
	ptr_buffer = NULL;
      
      }

    }
  
  
  
  
  
}



// void MAPPING_ID(const double A[], const int *npl, const double *dt, const double *a1, const double *a2, const double *a3, const double *b1, const double *b2, const double *b3) {
//    
//   volume_t *iter_volume;
//   
//   id_np_t *iter_id, *tmp_id;
//   
//   notionalParticle_t *ptr_particle;
//   
//   iter_volume = eulerian;
//   
//   volume_t svolume;
//     
//   while(iter_volume != NULL) {
//        
//     if(iter_volume->rvolume == NULL) {
//       
//       //MAPPING ID HERE **************************************
//       
//       HASH_ITER(hh, iter_volume->id_np, iter_id, tmp_id) {
//         
// 	      if( (notionalParticleArray[iter_particle->id].inew.igx != notionalParticleArray[iter_particle->id].index.igx) || (notionalParticleArray[iter_particle->id].inew.igy != notionalParticleArray[iter_particle->id].index.igy) || (notionalParticleArray[iter_particle->id].inew.igz != notionalParticleArray[iter_particle->id].index.igz) ) {
// 	
// // 	printf("\n\n\t\tIF");
// // 	gets(in);
// 	
// 	move_particle = malloc(sizeof(id_np_t));
// // 	printf("\n\n\t\tmove_particle = malloc(sizeof(id_np_t))");
// 	move_particle->id = iter_particle->id;
// // 	printf("\n\n\t\tmove_particle->id = iter_particle->id");
// 	move_volume = malloc(sizeof(volume_t));
// // 	printf("\n\n\t\tmove_volume = malloc(sizeof(volume_t))");
// // 	printf("\n\n\t\tmove_volume = malloc(sizeof(volume_t))");
// // 	move_volume = SEARCH_TABLE(0, notionalParticleArray[iter_particle->id].posix, notionalParticleArray[iter_particle->id].posiy, notionalParticleArray[iter_particle->id].posiz, A, levels);
// 	svolume = SEARCH_PARTICLE(A, notionalParticleArray[iter_particle->id].posix, notionalParticleArray[iter_particle->id].posiy, notionalParticleArray[iter_particle->id].posiz, ls);
// // 	printf("\n\n\t\tsvolume = SEARCH_PARTICLE()");
// 	move_volume = &svolume;
// // 	printf("\n\n\t\tmove_volume->svolume");
// 	HASH_ADD_INT(move_volume->id_np, id, move_particle);
// // 	printf("\n\n\t\tHASH_ADD_INT(move_volume->id_np, id, move_particle)");
// 	move_particle = NULL;
// // 	printf("\n\n\t\tmove_particle = NULL");
// 	free(move_particle);
// // 	printf("\n\n\t\tfree(move_particle)");
// 	
// 	move_volume = NULL;
// // 	printf("\n\n\t\tmove_volume = NULL");
// 	free(move_volume);
// // 	printf("\n\n\t\tfree(move_volume)");
// 	
// 	notionalParticleArray[iter_particle->id].index = notionalParticleArray[iter_particle->id].inew;
// // 	printf("\n\n\t\tnotionalParticleArray[iter_particle->id].index = notionalParticleArray[iter_particle->id].inew;");
// 
// // 	gets(in);
// 	
// 	tmp_particle = malloc(sizeof(id_np_t));
// // 	printf("\n\n\t\ttmp_particle = malloc(sizeof(id_np_t))");
// 	tmp_particle = iter_particle;
// // 	printf("\n\n\t\ttmp_particle = iter_particle");
// 	HASH_DEL(iter_volume->id_np, tmp_particle);
// // 	printf("\n\n\t\tHASH_DEL(iter_volume->id_np, tmp_particle)");
// 	tmp_particle = NULL;
// // 	printf("\n\n\t\ttmp_particle = NULL");
// 	free(tmp_particle);
// // 	printf("\n\n\t\tfree(tmp_particle)");
// 	
// // 	gets(in);
// 	
//       }
//       }
// 	
//       //***************************************************************
//       
//       if(iter_volume->hh.next != NULL)
// 	iter_volume = iter_volume->hh.next;  
//       
//       else if(iter_volume->hh.next == NULL && iter_volume->rl == 0) 
// 	 iter_volume = NULL;      
//       
//       else if(iter_volume->hh.next == NULL && iter_volume->rl != 0)  
// 	iter_volume = ptr_last[iter_volume->rl-1];
//       
//     }    
// 
//     else if(iter_volume->rvolume != NULL) {  
//       
//       ptr_last[iter_volume->rl] = iter_volume->hh.next;     
//       iter_volume = iter_volume->rvolume;      
//       
//     }    
//     
//   }
//     
//   for(i=0; i<levels; i++)
//     ptr_last[i] = NULL;
// 
// 
// }
//   
//   
//   
//   
//   
//   
//   
  
//   /*
// //   printf("\n\t\tMOVE_MICROMIXING");
//   
// //   int part_counter;
// //   
//   double delta;  
//   double zetap_S_mean;  
//   double weight_sum;
//   double omegam;
// //   
// //   double grad_gx, grad_gy, grad_gz;
// 
//   
//   long int idum3 = -3213853;
//   long int *ptr_idum3 = &idum3;
//   
//   volume_t *iter_volume, *tmp_volume, *move_volume;
//   
//   volume_t svolume;
//   
//   id_np_t *iter_id, *tmp_id, *move_id;
//     
//   notionalParticle_t *iter_particle, *tmp_particle;
//   
// //   printf("\n\n\t\tHASH_ITER(hh, eulerian, iter_volume, tmp_volume)");
// //   gets(in);
// 
//   HASH_ITER(hh, eulerian, iter_volume, tmp_volume) {
//     
//     part_counter = HASH_COUNT(iter_volume->id_np);
//     
//     delta = 0.0;  
//     zetap_S_mean = 0.0; 
//     weight_sum = 0.0;
//     omegam = 0.0;
//       
// //       printf("\n\n\t\tHASH_ITER(hh, iter_volume->id_np, iter_particle, tmp_particle)");
// //       gets(in);
//     HASH_ITER(hh, iter_volume->id_np, iter_particle, tmp_particle) {
// 	
// 	notionalParticleArray[iter_particle->id].weight = (rho_cte*dx*dy*dz)/part_counter; //Atualizar temporalmente em função de np_vol se alterar
// 	
// // 	printf("\n\n\t\t*** Volume %d %d %d - particle %d", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz, iter_particle->id);
// // 	gets(in);
// 	zetap_S_mean += (notionalParticleArray[iter_particle->id].zetap[0])*(notionalParticleArray[iter_particle->id].weight);
// // 	printf("\n\n\t\tzetap_S_mean = %f", zetap_S_mean);
// // 	gets(in);
// 	weight_sum += notionalParticleArray[iter_particle->id].weight;
// // 	printf("\n\n\t\tweight_sum = %f", weight_sum);
// // 	gets(in);
// // 	printf("\n\n\t\tSc = %f; ni_cte = %f; dx = %f; dy = %f; dz = %f", Sc, ni_cte, A[6], A[6+levels], A[6+2*levels]);
// // 	gets(in);
// 	
// 	grad_gx = (1/Sc)*(ni_cte)/A[6];
// 	grad_gy = (1/Sc)*(ni_cte)/A[6+levels];
// 	grad_gz = (1/Sc)*(ni_cte)/A[6+2*levels];
// // 	printf("\n\n\t\tgrad_gx = %f", grad_gx);
// // 	printf("\n\n\t\tgrad_gy = %f", grad_gy);
// // 	printf("\n\n\t\tgrad_gz = %f", grad_gz);
// // 
// // 	
// 	
// 	notionalParticleArray[iter_particle->id].posixo = notionalParticleArray[iter_particle->id].posix;
// 	notionalParticleArray[iter_particle->id].posiyo = notionalParticleArray[iter_particle->id].posiy;
// 	notionalParticleArray[iter_particle->id].posizo = notionalParticleArray[iter_particle->id].posiz;
// 	notionalParticleArray[iter_particle->id].posix = notionalParticleArray[iter_particle->id].posixo + U*dte + 0.0*dte*grad_gx + sqrt(dte*(two/Sc)*(ni_cte))*gasdev(ptr_idum3);
// 	notionalParticleArray[iter_particle->id].posiy = notionalParticleArray[iter_particle->id].posiyo + V*dte + 0.0*dte*grad_gy + sqrt(dte*(two/Sc)*(ni_cte))*gasdev(ptr_idum3);      
// 	notionalParticleArray[iter_particle->id].posiz = notionalParticleArray[iter_particle->id].posizo + W*dte + 0.0*dte*grad_gz + sqrt(dte*(two/Sc)*(ni_cte))*gasdev(ptr_idum3);
// 	
// // 	printf("\n\n\t\tposixo = %f", notionalParticleArray[iter_particle->id].posixo);
// // 	printf("\n\n\t\tposiyo = %f", notionalParticleArray[iter_particle->id].posiyo);
// // 	printf("\n\n\t\tposizo = %f", notionalParticleArray[iter_particle->id].posizo);
// // 	
// 	
// 	
// 	//Wall Boundary Conditions  
// 	
// 	if(notionalParticleArray[iter_particle->id].posix < a1)
// 	  notionalParticleArray[iter_particle->id].posix = a1 + (a1 - notionalParticleArray[iter_particle->id].posix);
// 	if(notionalParticleArray[iter_particle->id].posix > a2)
// 	  notionalParticleArray[iter_particle->id].posix = a2 - (notionalParticleArray[iter_particle->id].posix - a2);
//       
// 	if(notionalParticleArray[iter_particle->id].posiy < b1)
// 	  notionalParticleArray[iter_particle->id].posiy = b1 + (b1 - notionalParticleArray[iter_particle->id].posiy);
// 	if(notionalParticleArray[iter_particle->id].posiy > b2)
// 	  notionalParticleArray[iter_particle->id].posiy = b2 - (notionalParticleArray[iter_particle->id].posiy - b2);
//       
// 	if(notionalParticleArray[iter_particle->id].posiz < c1)
// 	  notionalParticleArray[iter_particle->id].posiz = c1 + (c1 - notionalParticleArray[iter_particle->id].posiz);
// 	if(notionalParticleArray[iter_particle->id].posiz > c2)
// 	  notionalParticleArray[iter_particle->id].posiz = c2 - (notionalParticleArray[iter_particle->id].posiz - c2);
//    /*       
// 	printf("\n\n\t\tposix = %f", notionalParticleArray[iter_particle->id].posix);
// 	printf("\n\n\t\tposiy = %f", notionalParticleArray[iter_particle->id].posiy);
// 	printf("\n\n\t\tposiz = %f", notionalParticleArray[iter_particle->id].posiz);
// 	gets(in);
// 	*/
// 	
// 	
// 	//Identifier (id_np_t) must be allocated in the correct volume (hash table) as a function of its own new position
// 	//Be carefull because the loop should visit the particle just once - so, the inew variable must be set here
//       
// 	svolume = SEARCH_PARTICLE(A, notionalParticleArray[iter_particle->id].posix, notionalParticleArray[iter_particle->id].posiy, notionalParticleArray[iter_particle->id].posiz, ls);
// 	notionalParticleArray[iter_particle->id].inew.igx = svolume.index.igx;
// 	notionalParticleArray[iter_particle->id].inew.igy = svolume.index.igy;
// 	notionalParticleArray[iter_particle->id].inew.igz = svolume.index.igz;
// 	
// // 	printf("\n\n\t\t*** Volume %d %d %d - particle %d", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz, iter_particle->id);
// // 	printf("\n\n\t\t*** New Volume %d %d %d - particle %d", notionalParticleArray[iter_particle->id].inew.igx, notionalParticleArray[iter_particle->id].inew.igy, notionalParticleArray[iter_particle->id].inew.igz, iter_particle->id);
// // 	gets(in);
// // 	
// 	
//       }
//       
//       
//       
//       zetap_S_mean = zetap_S_mean/weight_sum;
// //       printf("\n\n\t\tzetap_S_mean = %f", zetap_S_mean);
//        
// //       printf("\n\n\t\tparticles_counter = %d", particles_counter);
// //       gets(in);
//       HASH_ITER(hh, iter_volume->id_np, iter_particle, tmp_particle) {
// 	
// 	delta = pow(A[6]*A[6+levels]*A[6+2*levels], 1/3);
// 	// Turbulent frequency - Colluci 1998 and Fox 2005
// 	omegam = Comega*(ni_cte/Sc)/pow(delta, 2.0);
// 	// Micromixing model IEM
// 	notionalParticleArray[iter_particle->id].zetap[0] = notionalParticleArray[iter_particle->id].zetap[0] - omegam*dte*(notionalParticleArray[iter_particle->id].zetap[0] - zetap_S_mean);
//       
//       }
//           
//   }
//   
//   //Verify if the particle crosses the volume boundary using its inew (index_t) as parameter 
//   //and re-allocate its in the new hash table and delete it from the old hash table.
//   
// //   printf("\n\n\t\tCHANGING - HASH_ITER(hh, eulerian, iter_volume, tmp_volume)");
//   HASH_ITER(hh, eulerian, iter_volume, tmp_volume) {
//     
// 
// //     printf("\n\n\t\tfor(iter_particle=iter_volume->id_np; iter_particle!=NULL; iter_particle=iter_particle->hh.next)");
// 
//     for(iter_particle=iter_volume->id_np; iter_particle!=NULL; iter_particle=iter_particle->hh.next) {
//   
// //     printf("\n\n\t\t*** Volume %d %d %d - particle %d", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz, iter_particle->id);
//     
// //     gets(in);
//       
//       if( (notionalParticleArray[iter_particle->id].inew.igx != notionalParticleArray[iter_particle->id].index.igx) || (notionalParticleArray[iter_particle->id].inew.igy != notionalParticleArray[iter_particle->id].index.igy) || (notionalParticleArray[iter_particle->id].inew.igz != notionalParticleArray[iter_particle->id].index.igz) ) {
// 	
// // 	printf("\n\n\t\tIF");
// // 	gets(in);
// 	
// 	move_particle = malloc(sizeof(id_np_t));
// // 	printf("\n\n\t\tmove_particle = malloc(sizeof(id_np_t))");
// 	move_particle->id = iter_particle->id;
// // 	printf("\n\n\t\tmove_particle->id = iter_particle->id");
// 	move_volume = malloc(sizeof(volume_t));
// // 	printf("\n\n\t\tmove_volume = malloc(sizeof(volume_t))");
// // 	printf("\n\n\t\tmove_volume = malloc(sizeof(volume_t))");
// // 	move_volume = SEARCH_TABLE(0, notionalParticleArray[iter_particle->id].posix, notionalParticleArray[iter_particle->id].posiy, notionalParticleArray[iter_particle->id].posiz, A, levels);
// 	svolume = SEARCH_PARTICLE(A, notionalParticleArray[iter_particle->id].posix, notionalParticleArray[iter_particle->id].posiy, notionalParticleArray[iter_particle->id].posiz, ls);
// // 	printf("\n\n\t\tsvolume = SEARCH_PARTICLE()");
// 	move_volume = &svolume;
// // 	printf("\n\n\t\tmove_volume->svolume");
// 	HASH_ADD_INT(move_volume->id_np, id, move_particle);
// // 	printf("\n\n\t\tHASH_ADD_INT(move_volume->id_np, id, move_particle)");
// 	move_particle = NULL;
// // 	printf("\n\n\t\tmove_particle = NULL");
// 	free(move_particle);
// // 	printf("\n\n\t\tfree(move_particle)");
// 	
// 	move_volume = NULL;
// // 	printf("\n\n\t\tmove_volume = NULL");
// 	free(move_volume);
// // 	printf("\n\n\t\tfree(move_volume)");
// 	
// 	notionalParticleArray[iter_particle->id].index = notionalParticleArray[iter_particle->id].inew;
// // 	printf("\n\n\t\tnotionalParticleArray[iter_particle->id].index = notionalParticleArray[iter_particle->id].inew;");
// 
// // 	gets(in);
// 	
// 	tmp_particle = malloc(sizeof(id_np_t));
// // 	printf("\n\n\t\ttmp_particle = malloc(sizeof(id_np_t))");
// 	tmp_particle = iter_particle;
// // 	printf("\n\n\t\ttmp_particle = iter_particle");
// 	HASH_DEL(iter_volume->id_np, tmp_particle);
// // 	printf("\n\n\t\tHASH_DEL(iter_volume->id_np, tmp_particle)");
// 	tmp_particle = NULL;
// // 	printf("\n\n\t\ttmp_particle = NULL");
// 	free(tmp_particle);
// // 	printf("\n\n\t\tfree(tmp_particle)");
// 	
// // 	gets(in);
// 	
//       }
//            
//       
//     }    
//   }
//  
// }*/

  

  
// void HANDLE_PARTICLES(double stime, const double A[]) {
//      
//   printf("\n\n\t\t*** START HANDLE_PARTICLES() ***\n\n");
    
//   char in[1];
//   
//   int i = 0;
    
// //   int idum3 = -3213853;
// //   int *ptr_idum3 = &idum3;
  
  
// //   volume_t **ptr_last;
// //   ptr_last = malloc((levels-1)*sizeof(volume_t));  
// //   if(ptr_last == NULL) {
// //     printf("\nOut of memory\n");
// //     exit;
// //   }
    
// //   struct timeval tvBegin, tvEnd, tvDiff;
// //    
// //   gettimeofday(&tvBegin, NULL);


//   notionalParticle_t *notionalParticleArray;  
//   notionalParticleArray = malloc(((int)np_total+1)*sizeof(notionalParticle_t));
//   if(notionalParticleArray == NULL) {
//     printf("\nOut of memory\n");
//     exit;
//   }   
    
  /*INITIALIZE_PARTICLES_UNIFORMLY(A);
  printf("\nINITIALIZE_PARTICLES(A)\n");  
  gets(in); */ 

  //   INITIALIZE_PARTICLES(notionalParticleArray, A, ptr_last);

// //      PRINT_PARTICLES_ITERATOR(notionalParticleArray, A, ptr_last, stime);
  
//   INITIALIZE_ZETA(notionalParticleArray);
//   printf("\nINITIALIZE_ZETA(notionalParticleArray)\n");  
//   gets(in);  
    
//   PRINT_VTK(notionalParticleArray, 0.0);  
//   printf("\nPRINT_VTK(notionalParticleArray, totalParticles, 0.0)\n");  
//   gets(in);

  
  //PDF_PRINT(notionalParticleArray, 0.0);    
    
//     while(stime < FINAL_TIME) {
      
//       i++;      
//       stime += dte;
      //MOVE_MICROMIXING(notionalParticleArray, A, levels, ptr_idum3);
      
//       if(i % (int)(FINAL_TIME/(dte*nfiles)) == 0) {
	//PDF_PRINT(notionalParticleArray, stime);
// 	PRINT_VTK(notionalParticleArray, stime); 
//       }
      
//     }
    

  
// //   gettimeofday(&tvEnd, NULL);  
// //   timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
// //   printf("\n\n\t\tTime, in seconds, for PaSR with NPV (%i) = ", np_vol_ini);
// //   printf("%ld.%06ld\n\n", tvDiff.tv_sec, tvDiff.tv_usec);
 
//   printf("\n\n\t\t*** END ***\n\n");
//   
// }

/*
void INITIALIZE_PARTICLES(notionalParticle_t *notionalParticleArray, const double A[], volume_t **ptr_last) {

  int i, part_counter = 0;
  
  double x0, y0, z0;
  
  int idum3 = -3213853;
  int *ptr_idum3 = &idum3;
  
  volume_t *iter_volume;
  
  id_np_t *ptr_particle;
  
  iter_volume = eulerian;
  
  while(iter_volume != NULL) {
       
    if(iter_volume->rvolume == NULL) {
      
      //Insert your code here *************************
       
      x0 = (iter_volume->index.igx-1)*A[6+iter_volume->rl];
      y0 = (iter_volume->index.igy-1)*A[6+levels+iter_volume->rl];
      z0 = (iter_volume->index.igz-1)*A[6+2*levels+iter_volume->rl]; 
       
      for(i=1; i<=np_vol_ini; i++) {
	
	part_counter++;
	
	notionalParticleArray[i].index.igx = iter_volume->index.igx;
	notionalParticleArray[i].index.igy = iter_volume->index.igy;
	notionalParticleArray[i].index.igz = iter_volume->index.igz;
	
	notionalParticleArray[i].inew.igx = notionalParticleArray[i].index.igx;
	notionalParticleArray[i].inew.igy = notionalParticleArray[i].index.igy;
	notionalParticleArray[i].inew.igz = notionalParticleArray[i].index.igz;
            
	notionalParticleArray[i].posixo = x0 + ran2(ptr_idum3)*(A[6]+iter_volume->rl);
	notionalParticleArray[i].posiyo = y0 + ran2(ptr_idum3)*(A[6+levels]+iter_volume->rl);
	notionalParticleArray[i].posizo = z0 + ran2(ptr_idum3)*(A[6+2*levels+iter_volume->rl]);
	
        notionalParticleArray[i].posix = notionalParticleArray[i].posixo;
	notionalParticleArray[i].posiy = notionalParticleArray[i].posiyo;
	notionalParticleArray[i].posiz = notionalParticleArray[i].posizo;
	
	notionalParticleArray[i].treact = zeros;
	notionalParticleArray[i].dwall = zeros;
	notionalParticleArray[i].omegap = zeros;
	notionalParticleArray[i].rhop = rho_cte;
	notionalParticleArray[i].weight = rho_cte*(A[6]+iter_volume->rl)*(A[6+levels]+iter_volume->rl)*(A[6+2*levels+iter_volume->rl])/np_vol_ini; 
      
	ptr_particle = malloc(sizeof(id_np_t)); 
	ptr_particle->id = part_counter;      
	HASH_ADD_INT(iter_volume->id_np, id, ptr_particle);
	
      }
      //************************************************	
       
      if(iter_volume->hh.next != NULL)
	iter_volume = iter_volume->hh.next;  
      
      else if(iter_volume->hh.next == NULL && iter_volume->rl == 0) 
	 iter_volume = NULL;      
      
      else if(iter_volume->hh.next == NULL && iter_volume->rl != 0)  
	iter_volume = ptr_last[iter_volume->rl-1];
      
    }    

    else if(iter_volume->rvolume != NULL) {  
      
      ptr_last[iter_volume->rl] = iter_volume->hh.next;     
      iter_volume = iter_volume->rvolume;      
      
    }    
    
  }
    
  for(i=0; i<levels; i++)
    ptr_last[i] = NULL;
  


  
  
  
  
  
  
}*/





// void PRINT_PARTICLES_ITERATOR(notionalParticle_t *notionalParticleArray, const double A[], volume_t **ptr_last, double stime) {
// 
//   char in[1];
//   
//   int i;
//   
//   int part_counter;
//   
//   char file_name[50]; 
//   FILE *ptr_file; 
//   
//   sprintf(file_name, "Volumes_Content_%f.dat", stime);  
//   ptr_file = fopen(file_name, "w");  
//   fprintf(ptr_file, "Volume \t Level \n NP \n"); 
//   
//   volume_t *iter_volume;
//   
//   iter_volume = eulerian;
//   
//   while(iter_volume != NULL) {
//        
//     if(iter_volume->rvolume == NULL) {
//       
//       part_counter = HASH_COUNT(iter_volume->id_np);  
//       fprintf(ptr_file, "%i-%i-%i \t %i \t %i \n", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz, iter_volume->rl, part_counter);  
//        
//       if(iter_volume->hh.next != NULL)
// 	iter_volume = iter_volume->hh.next;  
//       
//       else if(iter_volume->hh.next == NULL && iter_volume->rl == 0) 
// 	 iter_volume = NULL;      
//       
//       else if(iter_volume->hh.next == NULL && iter_volume->rl != 0)  
// 	iter_volume = ptr_last[iter_volume->rl-1];
//       
//     }    
// 
//     else if(iter_volume->rvolume != NULL) {  
//       
//       ptr_last[iter_volume->rl] = iter_volume->hh.next;     
//       iter_volume = iter_volume->rvolume;      
//       
//     }    
//     
//   }
//     
//   for(i=0; i<levels; i++)
//     ptr_last[i] = NULL;
//   
//   fclose(ptr_file);
// 
// }





// 



// void MOVE_MICROMIXING(notionalParticle_t *notionalParticleArray, const double A[], const int ls, int *ptr_idum3) {
//   
//   
//   
// //   printf("\n\t\tMOVE_MICROMIXING");
//   
//   int part_counter;
//   
//   double delta;  
//   double zetap_S_mean;  
//   double weight_sum;
//   double omegam;
//   
//   double grad_gx, grad_gy, grad_gz;
//   
//   char in[1];
//   
//   volume_t *iter_volume, *tmp_volume, *move_volume;
//   
//   volume_t svolume;
//   
//   id_np_t *iter_particle, *tmp_particle, *move_particle;
//     
//   
// //   printf("\n\n\t\tHASH_ITER(hh, eulerian, iter_volume, tmp_volume)");
// //   gets(in);
// 
//   HASH_ITER(hh, eulerian, iter_volume, tmp_volume) {
//       
//       part_counter = HASH_COUNT(iter_volume->id_np);
//     
//       delta = zeros;  
//       zetap_S_mean = zeros; 
//       weight_sum = zeros;
//       omegam = zeros;
//       
// //       printf("\n\n\t\tHASH_ITER(hh, iter_volume->id_np, iter_particle, tmp_particle)");
// //       gets(in);
//       HASH_ITER(hh, iter_volume->id_np, iter_particle, tmp_particle) {
// 	
// 	notionalParticleArray[iter_particle->id].weight = (rho_cte*dx*dy*dz)/part_counter; //Atualizar temporalmente em função de np_vol se alterar
// 	
// // 	printf("\n\n\t\t*** Volume %d %d %d - particle %d", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz, iter_particle->id);
// // 	gets(in);
// 	zetap_S_mean += (notionalParticleArray[iter_particle->id].zetap[0])*(notionalParticleArray[iter_particle->id].weight);
// // 	printf("\n\n\t\tzetap_S_mean = %f", zetap_S_mean);
// // 	gets(in);
// 	weight_sum += notionalParticleArray[iter_particle->id].weight;
// // 	printf("\n\n\t\tweight_sum = %f", weight_sum);
// // 	gets(in);
// // 	printf("\n\n\t\tSc = %f; ni_cte = %f; dx = %f; dy = %f; dz = %f", Sc, ni_cte, A[6], A[6+levels], A[6+2*levels]);
// // 	gets(in);
// 	
// 	grad_gx = (1/Sc)*(ni_cte)/A[6];
// 	grad_gy = (1/Sc)*(ni_cte)/A[6+levels];
// 	grad_gz = (1/Sc)*(ni_cte)/A[6+2*levels];
// // 	printf("\n\n\t\tgrad_gx = %f", grad_gx);
// // 	printf("\n\n\t\tgrad_gy = %f", grad_gy);
// // 	printf("\n\n\t\tgrad_gz = %f", grad_gz);
// // 
// // 	
// 	
// 	notionalParticleArray[iter_particle->id].posixo = notionalParticleArray[iter_particle->id].posix;
// 	notionalParticleArray[iter_particle->id].posiyo = notionalParticleArray[iter_particle->id].posiy;
// 	notionalParticleArray[iter_particle->id].posizo = notionalParticleArray[iter_particle->id].posiz;
// 	notionalParticleArray[iter_particle->id].posix = notionalParticleArray[iter_particle->id].posixo + U*dte + 0.0*dte*grad_gx + sqrt(dte*(two/Sc)*(ni_cte))*gasdev(ptr_idum3);
// 	notionalParticleArray[iter_particle->id].posiy = notionalParticleArray[iter_particle->id].posiyo + V*dte + 0.0*dte*grad_gy + sqrt(dte*(two/Sc)*(ni_cte))*gasdev(ptr_idum3);      
// 	notionalParticleArray[iter_particle->id].posiz = notionalParticleArray[iter_particle->id].posizo + W*dte + 0.0*dte*grad_gz + sqrt(dte*(two/Sc)*(ni_cte))*gasdev(ptr_idum3);
// 	
// // 	printf("\n\n\t\tposixo = %f", notionalParticleArray[iter_particle->id].posixo);
// // 	printf("\n\n\t\tposiyo = %f", notionalParticleArray[iter_particle->id].posiyo);
// // 	printf("\n\n\t\tposizo = %f", notionalParticleArray[iter_particle->id].posizo);
// // 	
// 	
// 	
// 	//Wall Boundary Conditions  
// 	
// 	if(notionalParticleArray[iter_particle->id].posix < a1)
// 	  notionalParticleArray[iter_particle->id].posix = a1 + (a1 - notionalParticleArray[iter_particle->id].posix);
// 	if(notionalParticleArray[iter_particle->id].posix > a2)
// 	  notionalParticleArray[iter_particle->id].posix = a2 - (notionalParticleArray[iter_particle->id].posix - a2);
//       
// 	if(notionalParticleArray[iter_particle->id].posiy < b1)
// 	  notionalParticleArray[iter_particle->id].posiy = b1 + (b1 - notionalParticleArray[iter_particle->id].posiy);
// 	if(notionalParticleArray[iter_particle->id].posiy > b2)
// 	  notionalParticleArray[iter_particle->id].posiy = b2 - (notionalParticleArray[iter_particle->id].posiy - b2);
//       
// 	if(notionalParticleArray[iter_particle->id].posiz < c1)
// 	  notionalParticleArray[iter_particle->id].posiz = c1 + (c1 - notionalParticleArray[iter_particle->id].posiz);
// 	if(notionalParticleArray[iter_particle->id].posiz > c2)
// 	  notionalParticleArray[iter_particle->id].posiz = c2 - (notionalParticleArray[iter_particle->id].posiz - c2);
//    /*       
// 	printf("\n\n\t\tposix = %f", notionalParticleArray[iter_particle->id].posix);
// 	printf("\n\n\t\tposiy = %f", notionalParticleArray[iter_particle->id].posiy);
// 	printf("\n\n\t\tposiz = %f", notionalParticleArray[iter_particle->id].posiz);
// 	gets(in);
// 	*/
// 	
// 	
// 	//Identifier (id_np_t) must be allocated in the correct volume (hash table) as a function of its own new position
// 	//Be carefull because the loop should visit the particle just once - so, the inew variable must be set here
//       
// 	svolume = SEARCH_PARTICLE(A, notionalParticleArray[iter_particle->id].posix, notionalParticleArray[iter_particle->id].posiy, notionalParticleArray[iter_particle->id].posiz, ls);
// 	notionalParticleArray[iter_particle->id].inew.igx = svolume.index.igx;
// 	notionalParticleArray[iter_particle->id].inew.igy = svolume.index.igy;
// 	notionalParticleArray[iter_particle->id].inew.igz = svolume.index.igz;
// 	
// // 	printf("\n\n\t\t*** Volume %d %d %d - particle %d", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz, iter_particle->id);
// // 	printf("\n\n\t\t*** New Volume %d %d %d - particle %d", notionalParticleArray[iter_particle->id].inew.igx, notionalParticleArray[iter_particle->id].inew.igy, notionalParticleArray[iter_particle->id].inew.igz, iter_particle->id);
// // 	gets(in);
// // 	
// 	
//       }
//       
//       
//       
//       zetap_S_mean = zetap_S_mean/weight_sum;
// //       printf("\n\n\t\tzetap_S_mean = %f", zetap_S_mean);
//        
// //       printf("\n\n\t\tparticles_counter = %d", particles_counter);
// //       gets(in);
//       HASH_ITER(hh, iter_volume->id_np, iter_particle, tmp_particle) {
// 	
// 	delta = pow(A[6]*A[6+levels]*A[6+2*levels], 1/3);
// 	// Turbulent frequency - Colluci 1998 and Fox 2005
// 	omegam = Comega*(ni_cte/Sc)/pow(delta, 2.0);
// 	// Micromixing model IEM
// 	notionalParticleArray[iter_particle->id].zetap[0] = notionalParticleArray[iter_particle->id].zetap[0] - omegam*dte*(notionalParticleArray[iter_particle->id].zetap[0] - zetap_S_mean);
//       
//       }
//           
//   }
//   
//   //Verify if the particle crosses the volume boundary using its inew (index_t) as parameter 
//   //and re-allocate its in the new hash table and delete it from the old hash table.
//   
// //   printf("\n\n\t\tCHANGING - HASH_ITER(hh, eulerian, iter_volume, tmp_volume)");
//   HASH_ITER(hh, eulerian, iter_volume, tmp_volume) {
//     
// 
// //     printf("\n\n\t\tfor(iter_particle=iter_volume->id_np; iter_particle!=NULL; iter_particle=iter_particle->hh.next)");
// 
//     for(iter_particle=iter_volume->id_np; iter_particle!=NULL; iter_particle=iter_particle->hh.next) {
//   
// //     printf("\n\n\t\t*** Volume %d %d %d - particle %d", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz, iter_particle->id);
//     
// //     gets(in);
//       
//       if( (notionalParticleArray[iter_particle->id].inew.igx != notionalParticleArray[iter_particle->id].index.igx) || (notionalParticleArray[iter_particle->id].inew.igy != notionalParticleArray[iter_particle->id].index.igy) || (notionalParticleArray[iter_particle->id].inew.igz != notionalParticleArray[iter_particle->id].index.igz) ) {
// 	
// // 	printf("\n\n\t\tIF");
// // 	gets(in);
// 	
// 	move_particle = malloc(sizeof(id_np_t));
// // 	printf("\n\n\t\tmove_particle = malloc(sizeof(id_np_t))");
// 	move_particle->id = iter_particle->id;
// // 	printf("\n\n\t\tmove_particle->id = iter_particle->id");
// 	move_volume = malloc(sizeof(volume_t));
// // 	printf("\n\n\t\tmove_volume = malloc(sizeof(volume_t))");
// // 	printf("\n\n\t\tmove_volume = malloc(sizeof(volume_t))");
// // 	move_volume = SEARCH_TABLE(0, notionalParticleArray[iter_particle->id].posix, notionalParticleArray[iter_particle->id].posiy, notionalParticleArray[iter_particle->id].posiz, A, levels);
// 	svolume = SEARCH_PARTICLE(A, notionalParticleArray[iter_particle->id].posix, notionalParticleArray[iter_particle->id].posiy, notionalParticleArray[iter_particle->id].posiz, ls);
// // 	printf("\n\n\t\tsvolume = SEARCH_PARTICLE()");
// 	move_volume = &svolume;
// // 	printf("\n\n\t\tmove_volume->svolume");
// 	HASH_ADD_INT(move_volume->id_np, id, move_particle);
// // 	printf("\n\n\t\tHASH_ADD_INT(move_volume->id_np, id, move_particle)");
// 	move_particle = NULL;
// // 	printf("\n\n\t\tmove_particle = NULL");
// 	free(move_particle);
// // 	printf("\n\n\t\tfree(move_particle)");
// 	
// 	move_volume = NULL;
// // 	printf("\n\n\t\tmove_volume = NULL");
// 	free(move_volume);
// // 	printf("\n\n\t\tfree(move_volume)");
// 	
// 	notionalParticleArray[iter_particle->id].index = notionalParticleArray[iter_particle->id].inew;
// // 	printf("\n\n\t\tnotionalParticleArray[iter_particle->id].index = notionalParticleArray[iter_particle->id].inew;");
// 
// // 	gets(in);
// 	
// 	tmp_particle = malloc(sizeof(id_np_t));
// // 	printf("\n\n\t\ttmp_particle = malloc(sizeof(id_np_t))");
// 	tmp_particle = iter_particle;
// // 	printf("\n\n\t\ttmp_particle = iter_particle");
// 	HASH_DEL(iter_volume->id_np, tmp_particle);
// // 	printf("\n\n\t\tHASH_DEL(iter_volume->id_np, tmp_particle)");
// 	tmp_particle = NULL;
// // 	printf("\n\n\t\ttmp_particle = NULL");
// 	free(tmp_particle);
// // 	printf("\n\n\t\tfree(tmp_particle)");
// 	
// // 	gets(in);
// 	
//       }
//            
//       
//     }    
//   }
//  
// }


// void PDF_PRINT(notionalParticle_t *notionalParticleArray, double stime) {
//   
//   int i, j, k;
//   
//   int part_counter;
//   
//   double weight_sum;
//   
//   char in[1];
//     
//   volume_t *iter_volume, *tmp_volume;
//   
//   id_np_t *iter_particle, *tmp_particle;
//   
//   char file_name[100]; 
//   FILE *ptr_file; 
//   
//   char file_name2[100]; 
//   FILE *ptr_file2; 
//   
//   char file_name3[100]; 
//   FILE *ptr_file3; 
//   
//   double cdf[n_interv];    
//   double pdf[n_interv];
//   
//   double interv[n_interv+1];
//   double d_in = (zeta2_mix-zeta1_mix)/n_interv;
//     
//   interv[0] = zeta1_mix;
//   for(i=1; i<(n_interv+1); i++)     
//     interv[i] = interv[i-1] + d_in;
//   
//   sprintf(file_name2, "PV(%1.2f).dat", stime);  
//   ptr_file2 = fopen(file_name2, "w");  
//   fprintf(ptr_file2, "VOLUME \t NP \n");
//   
//   sprintf(file_name3, "WEIGHT(%1.2f).dat", stime);  
//   ptr_file3 = fopen(file_name3, "w");  
//   fprintf(ptr_file3, "VOLUME \t WEIGHT \n");
//     
//   HASH_ITER(hh, eulerian, iter_volume, tmp_volume) {     
// 
//     sprintf(file_name, "PDF_%i%i%i(%1.2f).dat", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz, stime);  
//     ptr_file = fopen(file_name, "w");  
//     fprintf(ptr_file, "NORM \t PDF \n"); 
//         
//     weight_sum = 0.0;
//     
// 
//     for(i=0; i<n_interv; i++) {
//       cdf[i] = 0.0;
//       pdf[i] = 0.0;
//     }
// 
//     HASH_ITER(hh, iter_volume->id_np, iter_particle, tmp_particle) {
//         
//       part_counter = HASH_COUNT(iter_volume->id_np);
// 
//       weight_sum += notionalParticleArray[iter_particle->id].weight;
//       
//       for(i=0; i<n_interv; i++) {
// 	if(notionalParticleArray[iter_particle->id].zetap[0] < interv[i+1]) {
// 	  cdf[i]++;  
// 	  break;
// 	}
//       }
//       
//     }
// 
//     
//     
//     for(i=0; i<n_interv; i++) {
//       pdf[i] = cdf[i]/part_counter;
//       fprintf(ptr_file, "%4.6f \t %4.6f \n", mid*(interv[i]+interv[i+1]), pdf[i]); 
//     }
//     
//     fprintf(ptr_file2, "%i-%i-%i \t %i \n", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz, part_counter);  
//     fprintf(ptr_file3, "%i-%i-%i \t %f \n", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz, weight_sum);   
//     
//     
//     fclose(ptr_file);  
//         
//   }
//   
//   fclose(ptr_file2);
//   fclose(ptr_file3);
//     
// }


  

// void data_distribution(char type, int counter[], int intervals, double delta_print) {
//   
//   int c;
//   
//   char file_name[30];
//   
//   sprintf(file_name, "ran2-%c-npv-%d.dat", type, np_vol_ini);
//   
//   FILE *file;
//   
//   file = fopen(file_name, "w");
//   
//   fprintf(file, "Mid Point \t Particles \n");
//   
//   for(c=0; c<intervals; c++)
//     fprintf(file, "%f \t %d \n", (c+0.5)*delta_print, counter[c]);
//   
//   fclose(file);
//   
// }


float ran1(int *idum) {
  
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (*idum <= 0 || !iy) { 								// Initialize.
  if (-(*idum) < 1) *idum=1; 								// Be sure to prevent idum = 0.
  else *idum = -(*idum);
  for (j=NTAB+7;j>=0;j--) { 								// Load the shuffle table (after 8 warm-ups).
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  if (j < NTAB) iv[j] = *idum;
  }
  iy=iv[0];
  }
  k=(*idum)/IQ; 									// Start here when not initializing.
  *idum=IA*(*idum-k*IQ)-IR*k; 								// Compute idum=(IA*idum) % IM without overif
  if(*idum < 0) *idum += IM; 								// flows by Schrage’s method.
  j=iy/NDIV_RAN1; 										// Will be in the range 0..NTAB-1.
  iy=iv[j]; 										// Output previously stored value and refill the
  iv[j] = *idum; 									// shuffle table.
  if ((temp=AM_RAN1*iy) > RNMX) return RNMX; 						// Because users don’t expect endpoint values.
  else return temp;
    
}


float ran2(long int *idum) {

  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (*idum <= 0) { 									//Initialize.
    if (-(*idum) < 1) *idum=1;								//Be sure to prevent idum = 0.
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { 								//Load the shuffle table (after 8 warm-ups).
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1; 									//Start here when not initializing.
  *idum=IA1*(*idum-k*IQ1)-k*IR1; 							//Compute idum=(IA1*idum) % IM1 without
  if (*idum < 0) *idum += IM1; 								//overflows by Schrage’s method.
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; 							//Compute idum2=(IA2*idum) % IM2 likewise.
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV; 										//Will be in the range 0..NTAB-1.
  iy=iv[j]-idum2; 									//Here idum is shuffled, idum and idum2 are
  iv[j] = *idum; 									//combined to generate output.
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;							// Because users don’t expect endpoint values.
  else return temp;

}  


float gasdev(int *idum) {

  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;
  if (*idum < 0) iset=0; 
  if (iset == 0) {
    do {
      v1=2.0*ran1(idum)-1.0; 
      v2=2.0*ran1(idum)-1.0; 
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0); 
    fac=sqrt(-2.0*log(rsq)/rsq);
 gset=v1*fac;
    iset=1;
    return v2*fac;
  } 
  else { 
    iset=0; 
    return gset;    
  } 
  
}

  



// void PRINT_GRID() {
//    
//   int i, j, k;
//   
//   double x, y, z;
//   
//   char string[10];
//   int file_number = 0;
//    
//   sprintf(string, "grid.vtk");
//     
//   FILE *pvtk;
//   pvtk = fopen(string, "w");
//   fprintf(pvtk, "# vtk DataFile Version 3.0\nStructured 3D Dataset\nASCII\nDATASET STRUCTURED_GRID\n");
//   fprintf(pvtk, "DIMENSIONS %d %d %d\n", nx+1, ny+1, nz+1);
//   fprintf(pvtk, "POINTS %d float\n", (nx+1)*(ny+1)*(nz+1));
//   
//   x = 0;
//   y = 0;
//   z = 0;
//   
//   for(i=1; i<=nx+1; i++)
//     for(j=1; j<=ny+1; j++)
//       for(k=1; k<=nz+1; k++) {
// 	x = (i-1)*dx;
// 	y = (j-1)*dy;
// 	z = (k-1)*dz;
// 	fprintf(pvtk, "%f %f %f\n", x, y, z);
//       }
//   
//   
//    
//   fclose(pvtk);
//      
// }


int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;

    return (diff<0);
}

