#ifndef CHIRAL_N2LO_3NF_FACTORIZED_H
#define CHIRAL_N2LO_3NF_FACTORIZED_H

#include "three_nucleon_force_model.h"

#include <string>

// Complete local N2LO spectator-1 force using the exact three-integral
// factorization of Hebeler et al., Phys. Rev. C 91, 044001 (2015), Eq. (6).
// Cartesian spin-momentum factors are expanded into a finite spherical-
// harmonic basis before the radial integrations.  No tensor or higher-partial-
// wave component is dropped.
class chiral_N2LO_3NF_factorized final : public three_nucleon_force_model
{
public:
	chiral_N2LO_3NF_factorized(double c_D, double c_E,
	                          double Lambda_3NF_MeV,
	                          double c1, double c3, double c4,
	                          int transfer_order);

	bool enabled() const override;
	std::string name() const override { return "chiral_N2LO_full_factorized"; }
	void update_parameters(const double* parameters) override;

	double W1_element(int alpha_r, int alpha_c,
	                  double p_r, double q_r,
	                  double p_c, double q_c,
	                  const pw_3N_statespace& pw_states) const override;

	double lec_c1_gev() const override { return m_c1_gev; }
	double lec_c3_gev() const override { return m_c3_gev; }
	double lec_c4_gev() const override { return m_c4_gev; }
	int angular_order_3nf() const override { return m_transfer_order; }
	double axial_coupling_3nf() const override;
	double pion_decay_constant_mev_3nf() const override;
	double pion_mass_mev_3nf() const override;
	double chiral_scale_mev_3nf() const override;
	double hbarc_mev_fm_3nf() const override;

private:
	double m_c_D;
	double m_c_E;
	double m_lambda;
	double m_c1_gev;
	double m_c3_gev;
	double m_c4_gev;
	int m_transfer_order;
};

#endif
