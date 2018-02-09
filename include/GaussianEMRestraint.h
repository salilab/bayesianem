/**
 *  \file IMP/bayesianem/GaussianEMRestraint.h
 *  \brief Restraint two sets of gaussians (model and gmm derived from EM map)
 *
 *  Copyright 2007-2016 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPBAYESIANEM_GAUSSIAN_EM_RESTRAINT_H
#define IMPBAYESIANEM_GAUSSIAN_EM_RESTRAINT_H

#include "bayesianem_config.h"
#include <IMP/PairContainer.h>
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/container_macros.h>
#include <IMP/core/XYZ.h>
#include <IMP/core/Gaussian.h>
#include <IMP/algebra/Gaussian3D.h>
#include <IMP/container/CloseBipartitePairContainer.h>
#include <IMP/container/ClosePairContainer.h>
#include <IMP/atom/Mass.h>
#include <math.h>
#include <Eigen/Dense>
#include <boost/unordered_map.hpp>

IMPBAYESIANEM_BEGIN_NAMESPACE

#if !defined(SWIG) && !defined(IMP_DOXYGEN)
struct KahanAccumulation{
double sum;
double correction;
KahanAccumulation():
  sum(0.0),correction(0.0)
  {}
};
struct KahanVectorAccumulation{
  Eigen::Vector3d sum;
  Eigen::Vector3d correction;
KahanVectorAccumulation():
  sum(Eigen::Vector3d(0,0,0)),
    correction(Eigen::Vector3d(0,0,0))
  {}
};
inline KahanAccumulation KahanSum(KahanAccumulation accumulation, double value){
  KahanAccumulation result;
  double y = value - accumulation.correction;
  double t = accumulation.sum + y;
  result.correction = (t - accumulation.sum) - y;
  result.sum = t;
  return result;
}
inline KahanVectorAccumulation
KahanVectorSum(KahanVectorAccumulation accumulation, Eigen::Vector3d value){
  KahanVectorAccumulation result;
  Eigen::Vector3d y = value - accumulation.correction;
  Eigen::Vector3d t = accumulation.sum + y;
  result.correction = (t - accumulation.sum) - y;
  result.sum = t;
  return result;
}
#endif


//! Creates a restraint between two Gaussian Mixture Models, "model" and "density"
//
/** This restrains two sets of GMMs with a density overlap function.
    The function is correlation of the two GMMs \f$f_M\f$ and \f$f_D\f$:
    \f[
       \frac{2\int{f_M(x)f_D(x)dx}}{\int{f_M^2(x)+f_D^2(x)dx}}
    \f]
    Where the integral is the "overlap function" given by:
    \f[
    ov(f_M,f_D) = \sum_{i=1}^{N_M} \sum_{j=1}^{N_D} \frac{1}{(2 \pi)^{3/2}|\Sigma_{Mi}+\Sigma_{Dj}|^{1/2}}\exp\left [-\frac{1}{2}(\boldsymbol\mu_{Mi} - \boldsymbol\mu_{Dj})^\top (\Sigma_{Mi}+\Sigma_{Dj})^{-1} (\boldsymbol\mu_{Mi} - \boldsymbol \mu_{Dj})\right ]
    \f]
    \note Source: Greenberg, Pellarin, Sali. In preparation.
 */
class IMPBAYESIANEMEXPORT GaussianEMRestraint : public Restraint
{

 public:
  //! Setup the GaussianEMRestraint
  /**
     \param[in] mdl the Model object to operate on
     \param[in] model_ps particles for the model GMM
     \param[in] density_ps particles for the density GMM
     \param[in] global_sigma Particle to modulate the uncertainty
     \param[in] model_cutoff_dist Cutoff for the model-model interactions
     \param[in] density_cutoff_dist Cutoff for model-density interactions
     \param[in] slope Gentle term to move all particles to the density
     \param[in] update_model (DEPRECATED) update model each cycle
     \param[in] backbone_slope Limit the slope only to backbone particles
     \param[in] name Name of this restraint
     \note the model and density particles must be set up as Gaussian
  */
  GaussianEMRestraint(Model *mdl,
                      ParticleIndexes model_ps, ParticleIndexes density_ps,
                      ParticleIndex global_sigma,
                      Float model_cutoff_dist,Float density_cutoff_dist,Float slope,
                      bool update_model=true, bool backbone_slope=false,
                      std::string name="GaussianEMRestraint%1%");

  //! Returns exp(score)
  double get_probability() const {
    return exp(-unprotected_evaluate(NULL));
  }

  //! Pre-calculate the density-density and model-model scores
  /** This is automatically called by the constructor.
      You only need to call it manually if you change Gaussian variances
  */
  void compute_initial_scores();
  
  ParticleIndexes const get_indexes(){ 
                         return pis_debug_;
			 };
  Floats const get_log2(){
                         return log2_debug_;
			 };
			 
  void debug();
  //! Set restraint slope
  void set_slope(Float s){slope_=s;}

  //! Get restraint slope
  Float get_slope(){return slope_;}

  virtual double
    unprotected_evaluate(IMP::DerivativeAccumulator *accum)
    const IMP_OVERRIDE;
  virtual IMP::ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE;
  void show(std::ostream &out) const { out << "GEM restraint"; }
  IMP_OBJECT_METHODS(GaussianEMRestraint);

 private:
  ParticleIndexes model_ps_;
  ParticleIndexes density_ps_;
  ParticleIndex global_sigma_;
  ParticleIndexes pis_debug_;
  Floats log2_debug_;
  Float slope_;
  bool update_model_;
  int msize_,dsize_;
  Float normalization_;
  Float dd_score_;
  Float self_mm_score_;
  PointerMember<container::CloseBipartitePairContainer> md_container_;
  Pointer<container::ClosePairContainer> mm_container_;
  ParticleIndexes slope_ps_; //experiment
  std::map<ParticleIndex,Float> map_score_dd_;
  Float cached_score_term_; 

  //variables needed to tabulate the exponential
  Floats exp_grid_;
  double invdx_;
  double argmax_;


};

IMPBAYESIANEM_END_NAMESPACE

#endif  /* IMPBAYESIANEM_GAUSSIAN_EM_RESTRAINT_H */
