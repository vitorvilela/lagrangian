void Grid_Sweep(notionalParticle_t *notionalParticleArray, const double A[], volume_t **ptr_last, double stime) {

  char in[1];
  
  int i;
  
  int part_counter;
  
  /*
  char file_name[50]; 
  FILE *ptr_file; 
    sprintf(file_name, "Volumes_Content_%f.dat", stime);  
  ptr_file = fopen(file_name, "w");  
  fprintf(ptr_file, "Volume \t Level \n NP \n"); 
  */
  
  
  volume_t *iter_volume;
  
  iter_volume = eulerian;
  
  while(iter_volume != NULL) {
       
    if(iter_volume->rvolume == NULL) {
      
      //Aqui! o ponteiro *iter_volume está varrendo as células ltop local
      
      /*
      part_counter = HASH_COUNT(iter_volume->id_np);  
      fprintf(ptr_file, "%i-%i-%i \t %i \t %i \n", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz, iter_volume->rl, part_counter);  
      */
      
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
  
  //fclose(ptr_file);

}