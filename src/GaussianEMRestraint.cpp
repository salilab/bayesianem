/**
 *  \file bayesianem/GaussianEMRestraint.cpp
 *  \brief Restraint two sets of gaussians (model and gmm derived from EM map)
 *
 *  Copyright 2007-2016 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/bayesianem/GaussianEMRestraint.h>
#include <IMP/math.h>
#include <IMP/atom/Atom.h>
#include <Eigen/LU>
#include <IMP/algebra/BoundingBoxD.h>
#include <IMP/algebra/vector_generators.h>
#include <IMP/isd/isd_config.h>
#include <IMP/isd/Scale.h>
#include <IMP/isd/em_utilities.h>
#include <limits>

IMPBAYESIANEM_BEGIN_NAMESPACE

GaussianEMRestraint::GaussianEMRestraint(Model *mdl, ParticleIndexes model_ps,
                                         ParticleIndexes density_ps,
                                         ParticleIndex global_sigma,
                                         Float model_cutoff_dist,
                                         Float cutoff_dist, Float slope,
                                         bool update_model, bool backbone_slope,
                                         std::string name)
    : Restraint(mdl, name), model_ps_(model_ps), density_ps_(density_ps),
      global_sigma_(global_sigma), slope_(slope), update_model_(update_model), cutoff_dist_(cutoff_dist) {

  msize_ = model_ps.size();
  dsize_ = density_ps.size();

  // check to make sure all particles are Gaussian and Mass
  for (int i = 0; i < msize_; i++) {
    IMP_USAGE_CHECK(core::Gaussian::get_is_setup(mdl, model_ps_[i]),
                    "Model particles must be Gaussian");
    IMP_USAGE_CHECK(atom::Mass::get_is_setup(mdl, model_ps_[i]),
                    "Model particles must have Mass");
  }
  for (int j = 0; j < dsize_; j++) {
    IMP_USAGE_CHECK(core::Gaussian::get_is_setup(mdl, density_ps_[j]),
                    "Density particles must be Gaussian");
    IMP_USAGE_CHECK(atom::Mass::get_is_setup(mdl, density_ps_[j]),
                    "Density particles must have Mass");
  }

	vec_score_dm_.reserve(dsize_);

	for (int j(0) ; j < dsize_ ; ++j){
		vec_score_dm_[j]=0.0;
	}

  compute_initial_scores();

  if (backbone_slope) {
    for (size_t nm = 0; nm < model_ps_.size(); nm++) {
      atom::AtomType a = atom::Atom(mdl, model_ps[nm]).get_atom_type();
      if (a == atom::AtomType("CA") || a == atom::AtomType("C") ||
          a == atom::AtomType("N")) {
        slope_ps_.push_back(model_ps_[nm]);
      }
    }
    std::cout << "limited slope to " << slope_ps_.size() << " ps out of "
              << model_ps_.size() << std::endl;
  } else
    slope_ps_ = model_ps_;
}

void GaussianEMRestraint::compute_initial_scores() {

  // precalculate DD score
  Eigen::Vector3d deriv;
  // this will be (3+N)/2*Log(2)-Log(Gamma(N/2)+Sum_i Log (overlap(d_i,D)) +
  // 0.5*Log(Pi)
  cached_score_term_ = 0.0;

	vec_score_dd_.reserve(dsize_);
  for (int i1 = 0; i1 < dsize_; i1++) {
		Float score_dd_i1 = 0.0;
    for (int i2 = 0; i2 < dsize_; i2++) {
      Float score = IMP::isd::score_gaussian_overlap(
          get_model(), ParticleIndexPair(density_ps_[i1], density_ps_[i2]),
          &deriv);
      score_dd_i1 += score;
    }
		vec_score_dd_.push_back(score_dd_i1);
  }

  for (int i1 = 0; i1 < dsize_; i1++) {
    cached_score_term_ += std::log(vec_score_dd_[i1]);
  }
  cached_score_term_ += (3.0 + dsize_) / 2.0 * std::log(2);
  cached_score_term_ += 0.5 * std::log(IMP::PI);
  cached_score_term_ -= std::lgamma(dsize_ / 2.0);

}

double
GaussianEMRestraint::unprotected_evaluate(DerivativeAccumulator *accum) const {
  Float scale = IMP::isd::Scale(get_model(), global_sigma_).get_scale();
  Eigen::Vector3d deriv;

  Float slope_score = 0.0;

  std::vector<Float> vec_score_dm_(dsize_, 0.0);

	double determinant;
	bool invertible;
	Eigen::Matrix3d inverse = Eigen::Matrix3d::Zero();
	for (int i(0) ; i < msize_ ; ++i){
		Float min_dist = std::numeric_limits<Float>::max();
		const ParticleIndex pi(model_ps_[i]);
		core::Gaussian g1(get_model(),pi);
		const Float m1 = atom::Mass(get_model(),pi).get_mass();
		Eigen::Matrix3d cov1 = g1.get_global_covariance();
		const Eigen::Vector3d v1(g1.get_coordinates().get_data());
		for (int j(0) ; j < dsize_ ; ++j){
			const int n = i*dsize_ + j;
			const ParticleIndex pj(density_ps_[j]);
			core::Gaussian g2(get_model(),pj);
			const Eigen::Vector3d v = Eigen::Vector3d(g2.get_coordinates().get_data()) - v1;
			const Float sd = v.norm();
			if (sd < min_dist) {
				min_dist = sd;
			}
			if (sd < cutoff_dist_){
				Eigen::Matrix3d covar = cov1 + g2.get_global_covariance();
				covar.computeInverseAndDetWithCheck(inverse,determinant,invertible);
				const Float mass12 = m1 * atom::Mass(get_model(),pj).get_mass();
				const Eigen::Vector3d tmp = inverse*v;
				const Float score = mass12 * 0.06349363593424097 / (std::sqrt(determinant)) * std::exp(-0.5*v.transpose()*tmp);
				vec_score_dm_[j] += score;
			}
		}
		slope_score += slope_ * min_dist;
	}


  double logterm = 0.0;
	for (int j(0) ; j < dsize_ ; ++j){
		const ParticleIndex pj(density_ps_[j]);
    const Float scoredd = vec_score_dd_[j];
    Float scoredm = vec_score_dm_[j];
    const Float lt =
        std::log(scale * (scoredm + 0.0000000001) / (scoredd + 0.0000000001));
    logterm += lt * lt;
  }
  const double log_score =
      dsize_ * 0.5 * std::log(logterm) + cached_score_term_ + slope_score;
  return log_score;
}


/* Return all particles whose attributes are read by the restraints. To
   do this, ask the pair score what particles it uses.*/
ModelObjectsTemp GaussianEMRestraint::do_get_inputs() const {
  ModelObjectsTemp ret;
  for (int i = 0; i < msize_; i++) {
    ret.push_back(get_model()->get_particle(model_ps_[i]));
  }
  for (int j = 0; j < dsize_; j++) {
    ret.push_back(get_model()->get_particle(density_ps_[j]));
  }
  ret.push_back(get_model()->get_particle(global_sigma_));
  return ret;
}

IMPBAYESIANEM_END_NAMESPACE
