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
#include <algorithm>

IMPBAYESIANEM_BEGIN_NAMESPACE

GaussianEMRestraint::GaussianEMRestraint(Model *mdl, ParticleIndexes model_ps,
                                         ParticleIndexes density_ps,
                                         ParticleIndex global_sigma,
                                         Float model_cutoff_dist,
                                         Float cutoff_dist, Float slope,
                                         bool update_model, bool backbone_slope,
                                         std::string name)
    : Restraint(mdl, name), model_ps_(model_ps), density_ps_(density_ps),
      global_sigma_(global_sigma), slope_(slope), update_model_(update_model), cutoff_dist_(cutoff_dist), serial_mover_(NULL) {

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

	vec_score_dm_.reserve(dsize_*msize_);
	for (int j(0) ; j < dsize_*msize_ ; ++j){
		vec_score_dm_.push_back(0.0);
	}

	slope_scores_.reserve(msize_);
	cached_model_indices_.reserve(msize_);
	for (int j(0) ; j < msize_ ; ++j){
		slope_scores_.push_back(0.0);
		cached_model_indices_.push_back(j);
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

	compute_individual_scores(cached_model_indices_);

}

Ints get_changed_indices(
	const IMP::core::SerialMover * sm,
	const std::vector<Ints> &indices,
	int num_ps){

	Ints ret;
	Ints changed_indices;

	const int i = sm->get_mover_index();
	return indices[i];
}


void
GaussianEMRestraint::compute_individual_scores(const Ints &changed_model_ps) const {

	double determinant;
	bool invertible;
	Eigen::Matrix3d inverse = Eigen::Matrix3d::Zero();
	for (int l(0) ; l < changed_model_ps.size() ; ++l){
		const int i=changed_model_ps[l];
		const ParticleIndex pi(model_ps_[i]);
		Float min_dist = std::numeric_limits<Float>::max();
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
			const Float xyzrdist = IMP::core::get_distance(IMP::core::XYZR(get_model(), pi), IMP::core::XYZR(get_model(), pj));
			if (xyzrdist < cutoff_dist_){
				Eigen::Matrix3d covar = cov1 + g2.get_global_covariance();
				covar.computeInverseAndDetWithCheck(inverse,determinant,invertible);
				const Float mass12 = m1 * atom::Mass(get_model(),pj).get_mass();
				const Eigen::Vector3d tmp = inverse*v;
				const Float score = mass12 * 0.06349363593424097 / (std::sqrt(determinant)) * std::exp(-0.5*v.transpose()*tmp);
				vec_score_dm_[n] = score;
			}
		}
		slope_scores_[i] = slope_ * min_dist;
	}
}

double
GaussianEMRestraint::unprotected_evaluate(DerivativeAccumulator *accum) const {
  Float scale = IMP::isd::Scale(get_model(), global_sigma_).get_scale();
  Eigen::Vector3d deriv;


	if(serial_mover_ && serial_mover_->get_mover_index()>=0){
		Ints changed_model_ps = get_changed_indices(serial_mover_, model_indexes_, model_ps_.size());
		compute_individual_scores(changed_model_ps);
	} else{
		compute_individual_scores(cached_model_indices_);
	}

	Float slope_score = 0.0;
	for (int i(0) ; i < msize_ ; ++i){
		slope_score+=slope_scores_[i];
	}

  double logterm = 0.0;
	for (int j(0) ; j < dsize_ ; ++j){
    const Float scoredd = vec_score_dd_[j];
    Float scoredm = 0.0;
		for (int i(0) ; i < msize_ ; ++i){
			const int n = i*dsize_ + j;
			scoredm+=vec_score_dm_[n];
		}
    const Float lt =
        std::log(scale * (scoredm + 0.0000000001) / (scoredd + 0.0000000001));
    logterm += lt * lt;
  }
  const double log_score =
      dsize_ * 0.5 * std::log(logterm) + cached_score_term_ + slope_score;
  return log_score;
}

Ints get_all_indexes(IMP::core::MonteCarloMover *m, const ParticleIndexes &ref){
	Ints ret;
	ModelObjectsTemp inputs = m->get_inputs();
	for(int i(0) ; i<inputs.size() ; ++i){
		ModelObject *o = inputs[i];
		Particle *p = dynamic_cast<Particle *>(o);
		if(p && IMP::core::Gaussian::get_is_setup(p)){
			for(int j(0) ; j<ref.size() ; ++j){
				if(p->get_index().get_index() == ref[j].get_index()){
					ret.push_back(j);
				}
			}
		}
		else if(p && IMP::core::RigidBody::get_is_setup(p)){
			ParticleIndexes pis = IMP::core::RigidBody(p).get_body_member_particle_indexes();
			for (int k(0) ; k<pis.size() ; ++k){
				for(int j(0) ; j<ref.size() ; ++j){
					if(pis[k].get_index() == ref[j].get_index()){
						ret.push_back(j);
					}
				}				
			}
		}
	}
	std::sort(ret.begin(), ret.end());
	return ret;
}

void GaussianEMRestraint::set_incremental(IMP::core::SerialMover *sm){
	serial_mover_ = sm;
	model_indexes_.clear();
	density_indexes_.clear();
	const IMP::core::MonteCarloMovers movers = sm->get_movers();
	for(int i(0) ; i<movers.size() ; ++i){
		IMP::core::MonteCarloMover *m=movers[i];
		Ints model_ps_in_mover = get_all_indexes(m, model_ps_);
		model_indexes_.push_back(model_ps_in_mover);
	}
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
