
#ifndef POTENTIAL_MODEL_H
#define POTENTIAL_MODEL_H

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "type_defs.h"

class potential_model
{
public:
	potential_model();
	virtual ~potential_model() = default;

	// [EN] Modern factory returning an owning unique_ptr. Production code should
	// prefer create(); fetch_potential_ptr() below is retained as a thin
	// compatibility wrapper that releases the pointer for legacy raw-pointer
	// call sites. / [CN] 新工厂返回持有所有权的 unique_ptr，生产代码应优先使用；
	// fetch_potential_ptr() 作为兼容包装保留，向旧式裸指针调用点释放所有权。
	static std::unique_ptr<potential_model> create(run_params run_parameters);
	static potential_model *fetch_potential_ptr(run_params run_parameters);

	virtual void first_parameter_sampling(bool statement) = 0;
	virtual void update_parameters(double* parameters) = 0;
	virtual void setup_store_matrices(double* p_mesh, int Np, bool coupled, int &S, int &J, int &T, int &Tz) = 0;

	virtual void V(int i, int j, double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray) = 0;
};

#endif // POTENTIAL_MODEL_H

