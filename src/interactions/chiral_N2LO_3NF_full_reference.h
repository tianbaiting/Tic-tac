#ifndef CHIRAL_N2LO_3NF_FULL_REFERENCE_H
#define CHIRAL_N2LO_3NF_FULL_REFERENCE_H

#include <complex>
#include <string>

#include "three_nucleon_force_model.h"

// Correctness-first projection of the complete local N2LO spectator-1 3NF.
//
// This class evaluates the five-angle Golak projection directly in Tic-tac's
// Jj-coupled basis, with explicit spin/isospin product states.  It is intended
// as an executable reference for validating a factorised production PWD.  Its
// O(N_angle^5) cost makes it unsuitable for a full WPCD cache at realistic
// grids, so the factory exposes it under an explicitly honest reference name.
class chiral_N2LO_3NF_full_reference : public three_nucleon_force_model
{
public:
	chiral_N2LO_3NF_full_reference(double c_D, double c_E,
	                               double Lambda_3NF_MeV,
	                               double c1_gev, double c3_gev,
	                               double c4_gev, int angular_order = 4);

	bool enabled() const override;
	std::string name() const override {
		return "chiral_N2LO_full_5d_reference";
	}

	double lec_c1_gev() const override { return m_c1_gev; }
	double lec_c3_gev() const override { return m_c3_gev; }
	double lec_c4_gev() const override { return m_c4_gev; }
	int angular_order() const { return m_angular_order; }
	int angular_order_3nf() const override { return m_angular_order; }
	double axial_coupling_3nf() const override;
	double pion_decay_constant_mev_3nf() const override;
	double pion_mass_mev_3nf() const override;
	double chiral_scale_mev_3nf() const override { return 700.0; }
	double hbarc_mev_fm_3nf() const override;

	double W1_element(int alpha_r, int alpha_c,
	                  double p_r, double q_r,
	                  double p_c, double q_c,
	                  const pw_3N_statespace& pw_states) const override;

private:
	double m_c_D;
	double m_c_E;
	double m_Lambda;
	double m_c1_gev;
	double m_c3_gev;
	double m_c4_gev;
	double m_c1_fm;
	double m_c3_fm;
	double m_c4_fm;
	double m_fpi;
	double m_mpi;
	double m_lambda_chi;
	int m_angular_order;
};

#endif
