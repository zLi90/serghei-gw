#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <pnetcdf.h>

static void handle_error(int status, int lineno) {
  
  fprintf(stderr, "Error at line %d: %s\n", lineno, ncmpi_strerror(status));
  MPI_Abort(MPI_COMM_WORLD, 1);
  
}

int main(int argc, char **argv) {

  int ncfile;
  int dimids[3];
  MPI_Offset st[3], ct[3];
  int ret; // error flag

  int var_ndims, var_natts;
  nc_type type;
  char filename[256], varname[NC_MAX_NAME+1];

  double *data = NULL;

  MPI_Init(&argc, &argv);

  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (argc > 2) {
    if (rank == 0) printf("Usage: %s filename\n", argv[0]);
    MPI_Finalize();
    exit(-1);
  }

  if (argc > 1) snprintf(filename, 256, "%s", argv[1]);
  else strcpy(filename, "testfile.nc");

  printf("Reading %s ...\n", filename);

  ret = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncfile);
  if (ret != NC_NOERR) handle_error(ret, __LINE__);
  printf("Opened %s ...\n", filename);
  printf("Looping through variables...\n");

  var_ndims = -1;
  var_natts = -1;
  strcpy(filename, "n/a");

  int ncols = -1;
  int nrows = -1;
  
  for (int i = 0; i < 9; i ++) {

    ret = ncmpi_inq_var(ncfile, i, varname, &type, &var_ndims, dimids, &var_natts);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);
    printf("%s | ndims: %d, natts: %d\n", varname, var_ndims, var_natts);
    
    if (!strcmp("single", varname)) {
      std::cout << "single : nothing to do." << std::endl;
    } else if (!strcmp("ncols", varname)) {
      std::cout << "ncols : set number of columns." << std::endl;
      st[0] = 0;
      ct[0] = var_ndims;
      ret = ncmpi_get_vara_int_all(ncfile, i, st, ct, &ncols);
      if (ret != NC_NOERR) handle_error(ret, __LINE__);
      std::cout << ncols << std::endl;
    } else if (!strcmp("nrows", varname)) {
      std::cout << "nrows : set number of rows." << std::endl;
      st[0] = 0;
      ct[0] = var_ndims;
      ret = ncmpi_get_vara_int_all(ncfile, i, st, ct, &nrows);
      if (ret != NC_NOERR) handle_error(ret, __LINE__);
      std::cout << nrows << std::endl;
    } else if (!strcmp("elevations", varname)) {
      std::cout << "elevations : set elevations." << std::endl;
      st[0] = 0;
      ct[0] = ncols * nrows;
      data = (double*) calloc(ncols*nrows, sizeof(double));
      ret = ncmpi_get_vara_double_all(ncfile, i, st, ct, data);
      if (ret != NC_NOERR) handle_error(ret, __LINE__);
      for (int k = 0; k < ncols*nrows; k ++) {
	printf("%d: %f\n", k, data[k]);
      }
    }
    
  }
  
  ret = ncmpi_close(ncfile);
  if (ret != NC_NOERR) handle_error(ret, __LINE__);
	
  MPI_Finalize();
  return 0;

}
