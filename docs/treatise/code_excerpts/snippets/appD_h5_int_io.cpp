// ===============================================================
// 抽取自仓库 [current]: src/io/disk_io_routines.cpp
// 行号区段：1286..1334
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void write_integer_to_h5(int integer, const char* int_name, hid_t file_id){
	hid_t   group_id;
	hid_t   dataset_id;
	herr_t  status;

	int     N_h5  [1];
	hsize_t dim_N [1] = {1};

	/* Open file and create/write content correponding to variable-name int_name */
	N_h5[0]     = integer;
	group_id    = H5Screate_simple(1, dim_N, NULL);
	dataset_id  = H5Dcreate(file_id, int_name, H5T_NATIVE_INT, group_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status      = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_h5); check_h5_write_call(status);
	
	/* Close dataset and group */
	status      = H5Dclose(dataset_id); check_h5_close_call(status);
	status      = H5Sclose(group_id); check_h5_close_call(status);
}
void read_integer_from_h5(int& integer, const char* int_name, const char* filename){

	hid_t  file_id;
	hid_t  dataset_id;
	herr_t status;
	file_id = H5Fopen(filename,
					  H5F_ACC_RDONLY,
					  H5P_DEFAULT);

	/* Open file and find content correponding to variable-name int_name */
	dataset_id = H5Dopen(file_id,
						 int_name,
						 H5P_DEFAULT);

	/* Read from file into N_h5 */
	int N_h5 [1];
	status = H5Dread (dataset_id,
					  H5T_NATIVE_INT,
					  H5S_ALL,
					  H5S_ALL,
					  H5P_DEFAULT,
					  N_h5);
	check_h5_read_call(status);

	/* Write value to input integer */
	integer = N_h5[0];

	/* Close file */
	status = H5Dclose(dataset_id); check_h5_close_call(status);
	status = H5Fclose(file_id); check_h5_close_call(status);
}
