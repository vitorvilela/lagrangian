
//DEBUG
printf("\n*ptr_seed  = %i\n", *ptr_seed);
//DEBUG

//DEBUG
printf("\nTransport_Stoch_Particles()\n");
//DEBUG

//DEBUG
printf("\nct = %i\n", *ct_seed);
//DEBUG

//DEBUG
char in[1];
//DEBUG

//DEBUG
printf("\ntransf : %d\n", transf);
//DEBUG

//DEBUG
printf("\nptr_old_vol = Search_Particle(A, iter_particle->posixo, iter_particle->posiyo, iter_particle->posizo, *npl)\n");
printf("\nHASH_FIND_INT(ptr_old_vol->id_particle, id, rm_id)\n");
printf("\nHASH_DEL(ptr_old_vol->id_particle, rm_id)\n");
printf("\nfree(rm_id)\n");
printf("\nrm_id = NULL\n");
printf("\nptr_particle = iter_particle\n");
printf("\nHASH_DEL(lagrangian, ptr_particle)\n");
printf("\nfree(ptr_particle)\n");
printf("\nptr_particle = NULL\n");
//DEBUG

//DEBUG
printf("\nMPI_Recv(&transf_counter, 1, MPI_INT, 0, FLAG1, MPI_COMM_WORLD, &status)  -  my_rank %d\n", my_rank);
printf("\ntransf_counter : %d, my_rank %d\n", transf_counter, my_rank);
//DEBUG

//DEBUG
printf("\nadd_id = malloc(sizeof(id_t))\n");
printf("\n\ntake id (%d) from stack - proc(%d)", id, my_rank);
printf("\nHASH_ADD_INT(ptr_volume->id_particle, id, add_id)\n");
//DEBUG

//DEBUG
printf("\n\npush(%d) - proc(%d)", sid, my_rank);
printf("\n\nvol_proc(%d) = %f \n\n", my_rank, vol_proc); 
printf("\n\nnp_proc(%d) = %d de np_total (%d)\n\n", my_rank, np_proc, np_total); 
printf("\n\np(%d) - proc(%d)", p, my_rank);
printf("\n\npop(%d) - proc(%d)", sid, my_rank);
//DEBUG






//DEBUG

//DEBUG







