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
#include <boost/math/special_functions/gamma.hpp>

IMPBAYESIANEM_BEGIN_NAMESPACE

GaussianEMRestraint::GaussianEMRestraint(Model *mdl, ParticleIndexes model_ps,
                                         ParticleIndexes density_ps,
                                         ParticleIndex global_sigma,
                                         Float model_cutoff_dist,
                                         Float density_cutoff_dist, Float slope,
                                         bool update_model, bool backbone_slope,
                                         std::string name)
    : Restraint(mdl, name), model_cutoff_dist_(model_cutoff_dist),
      density_cutoff_dist_(density_cutoff_dist),
      model_ps_(model_ps), density_ps_(density_ps),
      global_sigma_(global_sigma), slope_(slope), update_model_(update_model) {

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

  // Set up md container
  md_container_ = new container::CloseBipartitePairContainer(
      new container::ListSingletonContainer(mdl, model_ps),
      new container::ListSingletonContainer(mdl, density_ps),
      density_cutoff_dist);

  mm_container_ = new container::ClosePairContainer(
      new container::ListSingletonContainer(mdl, model_ps), model_cutoff_dist);

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
  dd_score_ = 0.0;
  self_mm_score_ = 0.0;
  // this will be (3+N)/2*Log(2)-Log(Gamma(N/2)+Sum_i Log (overlap(d_i,D)) +
  // 0.5*Log(Pi)
  cached_score_term_ = 0.0;

  for (int i1 = 0; i1 < dsize_; i1++) {
    map_score_dd_[density_ps_[i1]] = 0.0;
    for (int i2 = 0; i2 < dsize_; i2++) {
      Float score = IMP::isd::score_gaussian_overlap(
          get_model(), ParticleIndexPair(density_ps_[i1], density_ps_[i2]),
          &deriv);
      dd_score_ += score;
      map_score_dd_[density_ps_[i1]] += score;
    }
  }

  for (std::map<ParticleIndex, Float>::const_iterator iter =
           map_score_dd_.begin();
       iter != map_score_dd_.end(); ++iter) {
    cached_score_term_ += std::log(iter->second);
  }
  cached_score_term_ += (3.0 + dsize_) / 2.0 * std::log(2.0);
  cached_score_term_ += 0.5 * std::log(IMP::PI);
  cached_score_term_ -= boost::math::lgamma(dsize_ / 2.0);

  // precalculate the self-mm score and initialize
  for (int i = 0; i < msize_; i++) {
    Float score = IMP::isd::score_gaussian_overlap(
        get_model(), ParticleIndexPair(model_ps_[i], model_ps_[i]), &deriv);
    self_mm_score_ += score;
  }
  // std::cout<<"init dd: "<<dd_score_<<" init mm "<<self_mm_score_<<std::endl;
}

double
GaussianEMRestraint::unprotected_evaluate(DerivativeAccumulator *accum) const {
  // score is the square difference between two GMMs
  Float scale = IMP::isd::Scale(get_model(), global_sigma_).get_scale();
  KahanAccumulation md_score, mm_score;
  mm_score = KahanSum(mm_score, self_mm_score_);
  Eigen::Vector3d deriv;
  // std::cout<<"\neval"<<std::endl;
  boost::unordered_map<ParticleIndex, KahanVectorAccumulation> derivs_mm,
      derivs_md, slope_md;
  // boost::unordered_map<ParticleIndexPair,Float> md_dists(msize_*dsize_);

  Float slope_score = 0.0;

  if (slope_ > 0.0) {
    for (ParticleIndexes::const_iterator mit = slope_ps_.begin();
         mit != slope_ps_.end(); ++mit) {
      Float min_dist = -1.0;
      for (ParticleIndexes::const_iterator dit = density_ps_.begin();
           dit != density_ps_.end(); ++dit) {
        Eigen::Vector3d v =
            Eigen::Vector3d(
                core::XYZ(get_model(), *mit).get_coordinates().get_data()) -
            Eigen::Vector3d(
                core::XYZ(get_model(), *dit).get_coordinates().get_data());
        Float sd = v.norm();
        // md_dists[ParticleIndexPair(*mit,*dit)]=dist;
        if (min_dist <= 0.0) {
          min_dist = sd;
        } else if (sd < min_dist) {
          min_dist = sd;
        }
        slope_md[*mit] = KahanVectorSum(slope_md[*mit], v * slope_ / sd);
      }
      slope_score += slope_ * min_dist;
    }
  }

  std::map<ParticleIndex,Float> map_score_dm;
  // std::cout<<"calculating MD score"<<std::endl;
  IMP_CONTAINER_FOREACH(container::CloseBipartitePairContainer, md_container_, {
    Float score = IMP::isd::score_gaussian_overlap(get_model(), _1, &deriv);
 
    if (map_score_dm.find(_1[1]) != map_score_dm.end()) {
      map_score_dm[_1[1]] += score;
    } else {
      map_score_dm[_1[1]] = score;
    }
 
    md_score = KahanSum(md_score, score);
    if (accum) {
      derivs_md[_1[0]] = KahanVectorSum(derivs_md[_1[0]], -deriv);
    }
  });

  double logterm = 0.0;
  for (std::map<ParticleIndex, Float>::const_iterator iter =
           map_score_dd_.begin();
       iter != map_score_dd_.end(); ++iter) {
    Float scoredd = iter->second;
    Float scoredm = 0.0;
    if (map_score_dm.find(iter->first) != map_score_dm.end()) {
      scoredm = map_score_dm[iter->first];
    }
    Float lt =
        std::log(scale * (scoredm + 0.0000000001) / (scoredd + 0.0000000001));
    logterm += lt * lt;
  }
  double log_score =
      dsize_ * 0.5 * std::log(logterm) + cached_score_term_ + slope_score;
  return log_score;
}





void GaussianEMRestraint::debug() {
  // score is the square difference between two GMMs
  Float scale = IMP::isd::Scale(get_model(), global_sigma_).get_scale();
  KahanAccumulation md_score, mm_score;
  mm_score = KahanSum(mm_score, self_mm_score_);
  Eigen::Vector3d deriv;
  // std::cout<<"\neval"<<std::endl;
  boost::unordered_map<ParticleIndex, KahanVectorAccumulation> derivs_mm,
      derivs_md, slope_md;
  // boost::unordered_map<ParticleIndexPair,Float> md_dists(msize_*dsize_);
  
  Float slope_score = 0.0;
  log2_debug_.clear();
  pis_debug_.clear();

  if (slope_ > 0.0) {
    for (ParticleIndexes::const_iterator mit = slope_ps_.begin();
         mit != slope_ps_.end(); ++mit) {
      Float min_dist = -1.0;
      for (ParticleIndexes::const_iterator dit = density_ps_.begin();
           dit != density_ps_.end(); ++dit) {
        Eigen::Vector3d v =
            Eigen::Vector3d(
                core::XYZ(get_model(), *mit).get_coordinates().get_data()) -
            Eigen::Vector3d(
                core::XYZ(get_model(), *dit).get_coordinates().get_data());
        Float sd = v.norm();
        // md_dists[ParticleIndexPair(*mit,*dit)]=dist;
        if (min_dist <= 0.0) {
          min_dist = sd;
        } else if (sd < min_dist) {
          min_dist = sd;
        }
        slope_md[*mit] = KahanVectorSum(slope_md[*mit], v * slope_ / sd);
      }
      slope_score += slope_ * min_dist;
    }
  }

  std::map<ParticleIndex,Float> map_score_dm;
  // std::cout<<"calculating MD score"<<std::endl;
  IMP_CONTAINER_FOREACH(container::CloseBipartitePairContainer, md_container_, {
    Float score = IMP::isd::score_gaussian_overlap(get_model(), _1, &deriv);
 
    if (map_score_dm.find(_1[1]) != map_score_dm.end()) {
      map_score_dm[_1[1]] += score;
    } else {
      map_score_dm[_1[1]] = score;
    }
 
    md_score = KahanSum(md_score, score);
  });

  double logterm = 0.0;
  for (std::map<ParticleIndex, Float>::const_iterator iter =
           map_score_dd_.begin();
       iter != map_score_dd_.end(); ++iter) {
    Float scoredd = iter->second;
    Float scoredm = 0.0;
    if (map_score_dm.find(iter->first) != map_score_dm.end()) {
      scoredm = map_score_dm[iter->first];
    }
    Float lt =
        std::log(scale * (scoredm + 0.0000000001) / (scoredd + 0.0000000001));
    logterm += lt * lt;
    pis_debug_.push_back(iter->first);
    log2_debug_.push_back(lt*lt);
  }
  double log_score =
      dsize_ * 0.5 * std::log(logterm) + cached_score_term_ + slope_score;
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
  ret.push_back(md_container_);
  ret.push_back(mm_container_);
  return ret;
}

RestraintInfo *GaussianEMRestraint::get_static_info() const {
  IMP_NEW(RestraintInfo, ri, ());
  ri->add_string("type", "IMP.bayesianem.GaussianEMRestraint");
  if (!density_fn_.empty()) {
    ri->add_filename("filename", density_fn_);
  }
  ri->add_float("model cutoff", model_cutoff_dist_);
  ri->add_float("density cutoff", density_cutoff_dist_);
  return ri.release();
}

RestraintInfo *GaussianEMRestraint::get_dynamic_info() const {
  IMP_NEW(RestraintInfo, ri, ());
  Float scale = IMP::isd::Scale(get_model(), global_sigma_).get_scale();
  ri->add_float("global sigma", scale);
  return ri.release();
}

IMPBAYESIANEM_END_NAMESPACE
