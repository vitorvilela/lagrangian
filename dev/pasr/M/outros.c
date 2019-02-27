//   const int B[9] = {levels, 0, 1, 1, 1, nx, ny, nz, 1};
//   const double A[9] = {a1, a2, b1, b2, c1, c2, dx, dy, dz};
  
//   const int B[25] = {levels, 0, 1, 1, 1, 4, 1, 4, 1, 1, 3, 1, 3, 4, 1, 4, 17, 2, 7, 1, 7, 4, 1, 4, 33};
//   const double A[15] = {a1, a2, b1, b2, c1, c2, 0.25, 0.125, 0.0625, 0.25, 0.25, 0.25, 0.25, 0.125, 0.0625};
//   

volume_t* SEARCH_RECURSIVELY(volume_t *ptr_volume) {
  
  char in[1];
     
  
  if(ptr_volume->rvolume == NULL)
    return ptr_volume;
  else
    return SEARCH_RECURSIVELY(ptr_volume->rvolume);
  
}