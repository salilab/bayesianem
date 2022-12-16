/**
 *  \file IMP/bayesianem/utilities.h
 *  \brief Useful utilities
 *
 *  Copyright 2007-2016 IMP Inventors. All rights reserved.
 */

#ifndef IMPBAYESIANEM_UTILITIES_H
#define IMPBAYESIANEM_UTILITIES_H

#include <IMP/bayesianem/bayesianem_config.h>
//#include <IMP/algebra/algebra_config.h>
#include "IMP/algebra/Gaussian3D.h"
#include <IMP/algebra/grid_utility.h>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include <IMP/em.h>
#include <IMP/isd.h>

#include <limits>
#include <vector>

IMPBAYESIANEM_BEGIN_NAMESPACE

typedef IMP::algebra::DenseGrid3D<double> DensityGrid;

double get_rmsd_excluding_higher_than( const IMP::core::XYZs &m1,
                                              const IMP::core::XYZs &m2,
                                              const double t) {
  double rmsd = 0.0;
  int N(0);
  for (int i=0; i<m1.size(); ++i) {
    double tmprmsd = get_squared_distance(m1[i].get_coordinates(), m2[i].get_coordinates());
    if(tmprmsd<t*t){
        rmsd+=tmprmsd;
        ++N;
    }
  }
  return std::sqrt(rmsd / N);
  //return Floats({std::sqrt(rmsd / N),((double)N/m1.size())});
}

double get_percentage_closer_than(const IMP::core::XYZs &m1,
                                  const IMP::core::XYZs &m2,
                                  const double t) {
  int N(0);
  for (int i=0; i<m1.size(); ++i) {
    double tmprmsd = get_squared_distance(m1[i].get_coordinates(), m2[i].get_coordinates());
    if(tmprmsd<t*t){
        ++N;
    }
  }
  return ((double) N)/m1.size();
  //return Floats({std::sqrt(rmsd / N),((double)N/m1.size())});
}

double get_rmsd_of_best_population(const IMP::core::XYZs &m1,
                                   const IMP::core::XYZs &m2,
                                   const double percentage){
	std::vector<double> sq_distances(m1.size());
	for (int i=0; i<m1.size(); ++i) {
		sq_distances[i] = get_squared_distance(m1[i].get_coordinates(), m2[i].get_coordinates());
	}
	std::sort(sq_distances.begin(), sq_distances.end());
	double sd=0.0;
	int N=0;
	for (int i=0; ((double) i) < (percentage*m1.size()); ++i) {
		sd+=sq_distances[i];
		++N;
	}
	return std::sqrt(sd/N);
}

double get_rmsd_of_best_population(const IMP::atom::Selection &s1,
                                   const IMP::atom::Selection &s2,
                                   const double percentage){

	IMP::Particles ps1=s1.get_selected_particles();
	IMP::Particles ps2=s2.get_selected_particles();

	IMP::core::XYZs ds1;
	IMP::core::XYZs ds2;

	for (int i=0; i<ps1.size(); ++i){ds1.push_back(IMP::core::XYZ(ps1[i]));}
	for (int i=0; i<ps2.size(); ++i){ds2.push_back(IMP::core::XYZ(ps2[i]));}

	return get_rmsd_of_best_population(ds1,ds2,percentage);
}

IMP::algebra::Transformation3D get_transformation_aligning_first_to_second(const IMP::atom::Selection &s1,
                                   const IMP::atom::Selection &s2){

	IMP::Particles ps1=s1.get_selected_particles();
	IMP::Particles ps2=s2.get_selected_particles();

	IMP::core::XYZs ds1;
	IMP::core::XYZs ds2;

	for (int i=0; i<ps1.size(); ++i){ds1.push_back(IMP::core::XYZ(ps1[i]));}
	for (int i=0; i<ps2.size(); ++i){ds2.push_back(IMP::core::XYZ(ps2[i]));}

	return IMP::algebra::get_transformation_aligning_first_to_second(ds1,ds2);
}



double get_rmsd_of_best_population(const IMP::algebra::Vector3Ds &m1,
                                   const IMP::algebra::Vector3Ds &m2,
                                   const double percentage){
	std::vector<double> sq_distances(m1.size());
	for (int i=0; i<m1.size(); ++i) {
		sq_distances[i] = get_squared_distance(m1[i], m2[i]);
	}
	std::sort(sq_distances.begin(), sq_distances.end());
	double sd=0.0;
	int N=0;
	for (int i=0; ((double) i) < (percentage*m1.size()); ++i) {
		sd+=sq_distances[i];
		++N;
	}
	return std::sqrt(sd/N);
}

double gem_score_cc(Particles model_ps, Particles density_ps){
	Eigen::Vector3d deriv;
	double mm_score(0.0);
	double md_score(0.0);
	double dd_score(0.0);

	int nm = model_ps.size();
	int nd = density_ps.size();
	IMP::Model *mdl = model_ps[0]->get_model();
	
	for(int nd1=0 ; nd1<nd ; ++nd1){
		for(int nd2=0 ; nd2<nd ; ++nd2){
			dd_score += IMP::isd::score_gaussian_overlap(mdl, ParticleIndexPair(density_ps[nd1]->get_index(), density_ps[nd2]->get_index()), &deriv);
		}
	}
	for(int nm1=0 ; nm1<nm ; ++nm1){
		for(int nm2=0 ; nm2<nm ; ++nm2){
			mm_score += IMP::isd::score_gaussian_overlap(mdl, ParticleIndexPair(model_ps[nm1]->get_index(), model_ps[nm2]->get_index()), &deriv);
		}
		for(int nd1=0 ; nd1<nd ; ++nd1){
			md_score += IMP::isd::score_gaussian_overlap(mdl, ParticleIndexPair(model_ps[nm1]->get_index(), density_ps[nd1]->get_index()), &deriv);
		}
	}
	double cc = md_score/std::sqrt(mm_score*dd_score);
	return cc;
}


namespace {
double get_gaussian_eval_prefactor(double determinant) {
  return 1.0 / pow(2 * algebra::PI, 2.0 / 3.0) / std::sqrt(determinant);
}
Eigen::Vector3d get_vec_from_center(IMP::algebra::GridIndex3D i,
                                        DensityGrid const &g,
                                        Eigen::Vector3d const &center) {
  IMP::algebra::Vector3D aloc = g.get_center(i);
  Eigen::Vector3d loc(aloc[0], aloc[1], aloc[2]);
  Eigen::Vector3d r(loc - center);
  return r;
}
void update_value(DensityGrid *g,
                  DensityGrid::Index i, Eigen::Vector3d r,
                  Eigen::Matrix3d inverse, double pre, Float weight) {
  double d(r.transpose() * (inverse * r));
  double score(pre * weight * std::exp(-0.5 * (d)));
  (*g)[i] += score;
}
}

DensityGrid get_grid(IMP::em::DensityMap *in) {
  IMP_FUNCTION_LOG;
  IMP_CHECK_OBJECT(in);
  DensityGrid ret(in->get_header()->get_spacing(), get_bounding_box(in));
  IMP_USAGE_CHECK(ret.get_number_of_voxels(0) ==
                      static_cast<unsigned int>(in->get_header()->get_nx()),
                  "X voxels don't match");
  IMP_USAGE_CHECK(ret.get_number_of_voxels(1) ==
                      static_cast<unsigned int>(in->get_header()->get_ny()),
                  "Y voxels don't match");
  IMP_USAGE_CHECK(ret.get_number_of_voxels(2) ==
                      static_cast<unsigned int>(in->get_header()->get_nz()),
                  "Z voxels don't match");
  for (unsigned int i = 0; i < ret.get_number_of_voxels(0); ++i) {
    for (unsigned int j = 0; j < ret.get_number_of_voxels(1); ++j) {
      for (unsigned int k = 0; k < ret.get_number_of_voxels(2); ++k) {
        DensityGrid::ExtendedIndex ei(i, j, k);
        DensityGrid::Index gi = ret.get_index(ei);
        long vi = in->get_voxel_by_location(i, j, k);
        ret[gi] = in->get_value(vi);
      }
    }
  }
  return ret;
}

IMP::em::DensityMap *get_masked_map(const IMP::algebra::Gaussian3Ds &gmm,
                                    const Floats &weights,
                                    IMP::em::DensityMap *densitymap,
                                    double threshold) {
  DensityGrid mask = IMP::bayesianem::get_grid(densitymap);

  for (const DensityGrid::Index & i : mask.get_all_indexes()) {
    mask[i] = 0.0;
  }

  for (unsigned int ng = 0; ng < gmm.size(); ng++) {
    Eigen::Matrix3d covar = get_covariance(gmm[ng]);
    Eigen::Matrix3d inverse = Eigen::Matrix3d::Zero(3, 3);

    double determinant;
    bool invertible;
    covar.computeInverseAndDetWithCheck(inverse, determinant, invertible);
    IMP_INTERNAL_CHECK((!invertible || determinant < 0),
                       "\n\n\n->>>>not proper matrix!!\n\n\n");
    double pre(get_gaussian_eval_prefactor(determinant));
    Eigen::Vector3d evals = covar.eigenvalues().real();
    double maxeval = sqrt(evals.maxCoeff());
    double cutoff = 2.5 * maxeval;
    double cutoff2 = cutoff * cutoff;
    IMP::algebra::Vector3D c = gmm[ng].get_center();
    IMP::algebra::Vector3D lower =
        c - IMP::algebra::Vector3D(cutoff, cutoff, cutoff);
    IMP::algebra::Vector3D upper =
        c + IMP::algebra::Vector3D(cutoff, cutoff, cutoff);
    IMP::algebra::GridIndex3D lowerindex = mask.get_nearest_index(lower);
    IMP::algebra::GridIndex3D upperindex = mask.get_nearest_index(upper);
    Eigen::Vector3d center(c.get_data());
    IMP_INTERNAL_CHECK(invertible, "matrix wasn't invertible! uh oh!");
    IMP_GRID3D_FOREACH_SMALLER_EXTENDED_INDEX_RANGE(mask, upperindex,
                                                    lowerindex, upperindex, {
      IMP::algebra::GridIndex3D i(voxel_index[0], voxel_index[1],
                                  voxel_index[2]);
      Eigen::Vector3d r(get_vec_from_center(i, mask, center));
      if (r.squaredNorm() < cutoff2) {
        update_value(&mask, i, r, inverse, pre, weights[ng]);
      }
    })
  }
  Pointer<IMP::em::DensityMap> maskmap = IMP::em::create_density_map(mask);
  Pointer<IMP::em::DensityMap> ret = IMP::em::multiply(maskmap, densitymap);
  return ret.release();
}

IMP::em::DensityMap *get_sub_map(const IMP::em::DensityMap *dm,
                                 const IMP::em::DensityMap *sub_gmm,
                                 const IMP::em::DensityMap *gmm) {
  const IMP::em::DensityHeader *header = sub_gmm->get_header();
  Pointer<IMP::em::DensityMap> m_map(IMP::em::create_density_map(dm));
  double *data1 = sub_gmm->get_data();
  double *data2 = gmm->get_data();
  double *new_data = m_map->get_data();
  for (long i = 0; i < header->get_number_of_voxels(); i++) {
    if (data2[i] != 0.0) {
      double const w(data1[i] / data2[i]);
      new_data[i] *= w;
    } else {
      new_data[i] = 0.0;
    }
  }
  return m_map.release();
}

double get_overlap_fast(const IMP::algebra::Gaussian3Ds &gmm,
                        const Floats &weights, IMP::em::DensityMap *dm,
                        double factor = 2.5) {
  DensityGrid rasterized = IMP::algebra::get_rasterized_fast(
      gmm, weights, dm->get_header()->get_spacing(), get_bounding_box(dm),
      factor);
  DensityGrid orig = IMP::bayesianem::get_grid(dm);
  /*
  double m1(0.0);
  double m2(0.0);
  for (const IMP::algebra::DensityGrid::Index & i :
  orig.get_all_indexes()) {
          IMP::algebra::Vector3D position(orig.get_center(i));
          IMP::algebra::DensityGrid::Index
  j(rasterized.get_nearest_index(position));
          m1+=rasterized[j];
          m2+=orig[i];
  }
  double scale = m2/m1;
  */
  double score(0.0);
  int norm(0);
  for (const DensityGrid::Index & i : orig.get_all_indexes()) {
    IMP::algebra::Vector3D position(orig.get_center(i));
    DensityGrid::Index j(rasterized.get_nearest_index(position));
    // rasterized[j]*=scale;
    double x = rasterized[j] - orig[i];
    score += x * x;
    ++norm;
  }
  return score / norm;
}

IMP::algebra::Vector3Ds get_overlap_binned(const IMP::algebra::Gaussian3Ds &gmm,
                        const Floats &weights, IMP::em::DensityMap *dm,
                        double factor = 2.5, int Nbin=100) {
  
	Floats score(Nbin,0.);
	Floats densities(Nbin,0.);
	Ints counts(Nbin,0);

	double max(0.0);
	double min(std::numeric_limits<double>::max());
	for (long i=0 ; i<dm->get_number_of_voxels() ; ++i){
		double const rho = dm->get_value(i);
		if (rho>0){
			//std::cerr << "rho " << i << " " << rho << "\n";
		}
		if (rho>max){
			max=rho;
		}
		if ((rho>0.0) && (rho<min)){
			min=rho;
		}
	}
	double const dx( (max-min)/(Nbin) );
	//std::cerr << "min " << min << "\nmax " << max << "\ndx " << dx << "\n";

	for(int i(0) ; i<Nbin ; ++i){
		densities[i]=min+dx*i;
	}

  DensityGrid rasterized = IMP::algebra::get_rasterized_fast(
      gmm, weights, dm->get_header()->get_spacing(), get_bounding_box(dm),
      factor);
	IMP::em::DensityMap* rasterized_map = IMP::em::create_density_map(rasterized);
	for (long i=0 ; i<dm->get_number_of_voxels() ; ++i){
		double const rho = dm->get_value(i);
		if(rho>=min){
			size_t const index = floor((rho-min)/dx);
			//IMP::algebra::Vector3D const position = dm->get_location_by_voxel(i);
			double const gmm_val = rasterized_map->get_value(i);
			double const x = gmm_val - rho;
			score[index] += x * x;
			++counts[index];
		}
	}
	IMP::algebra::Vector3Ds ret(Nbin);
	for(int  i(0) ; i<Nbin ; ++i){
		if(counts[i]>0){
			score[i]/=counts[i];
		}
		ret[i][0]=densities[i];
		ret[i][1]=score[i];
		ret[i][2]=counts[i];
	}
	return ret;
}

IMP::algebra::Rotation3D get_rotation_matrix(const IMP::algebra::Vector3D& x,
                                          const IMP::algebra::Vector3D& y) {
    IMP::algebra::Vector3D z = IMP::algebra::get_vector_product(x, y);
    return IMP::algebra::get_rotation_from_matrix(x[0], x[1], x[2], y[0], y[1], y[2],
                                             z[0], z[1], z[2]);
}

float sgn(double v) {
    return ( v < 0.0 ) ? -1.0 : ( ( v > 0.0 ) ? 1.0 : 0 );
} 

IMP::algebra::PrincipalComponentAnalysis NormalizePCA(const IMP::algebra::PrincipalComponentAnalysis& pca, 
   const IMP::Particles& ps)
{
                   // prepare vector for new pcas
		   IMP::algebra::Vector3Ds newpcs; 
                   // cycle on pca
		   for (unsigned int i = 0; i < 3; ++i) {
		       // get pca
		       IMP::algebra::Vector3D x = pca.get_principal_component(i);
		       // calculate the sign of the sum of the signed inner products
		       double ips = 0.0;
		       // cycle on all particles
		       for (unsigned int j = 0 ; j < ps.size() ; ++j ) {
			   // scalar product
			   IMP::algebra::Vector3D y = IMP::core::XYZ(ps[j]).get_coordinates()-pca.get_centroid();
			   double ip = x.get_scalar_product(y); 
			   // increment ips
			   ips += sgn(ip) * ip * ip;
		       }
		       // check sign of ips
		       if ( ips < 0. )  x *=  -1.0;
		       // store pca
		       newpcs.push_back(x);
		   }
		   // create new pcas
		   IMP::algebra::PrincipalComponentAnalysis newpcas =
		   IMP::algebra::PrincipalComponentAnalysis(newpcs, pca.get_principal_values(), pca.get_centroid());
		   
		   return newpcas;        
}


// MB old 

IMP::algebra::Transformation3Ds PCAalign(const IMP::algebra::PrincipalComponentAnalysis& pca1, 
                                         const IMP::algebra::PrincipalComponentAnalysis& pca2)
     {
  algebra::Transformation3Ds transforms;

  IMP::algebra::Vector3D x = pca2.get_principal_component(0);
  IMP::algebra::Vector3D y = pca2.get_principal_component(1);
  IMP::algebra::Rotation3D r = IMP::bayesianem::get_rotation_matrix(x, y);
  IMP::algebra::Transformation3D map_trans(
      IMP::algebra::Transformation3D(r, -(r * pca2.get_centroid())));
  IMP::algebra::Transformation3D inverse_map_trans = map_trans.get_inverse();

  // align the principal components by enumeration 6 xy choices

  IMP::algebra::Vector3D xi = pca1.get_principal_component(0);
  IMP::algebra::Vector3D yi = pca1.get_principal_component(1);
  IMP::algebra::Rotation3Ds rotations(1);
  rotations[0] = IMP::bayesianem::get_rotation_matrix(xi, yi);
  for (unsigned int k = 0; k < rotations.size(); k++) {
    IMP::algebra::Transformation3D points_trans =
        IMP::algebra::Transformation3D(
            rotations[k], -(rotations[k] * pca1.get_centroid()));
    IMP::algebra::Transformation3D ps2dens =
        inverse_map_trans * points_trans;
    transforms.push_back(ps2dens);
  }
  return transforms;
}


/*
IMP::algebra::Transformation3D PCAalign(const IMP::algebra::PrincipalComponentAnalysis& pca1, 
                                        const IMP::algebra::PrincipalComponentAnalysis& pca2)
{
      // get the two sets of pcas
      IMP::algebra::Vector3Ds pca1v = pca1.get_principal_components();
      IMP::algebra::Vector3Ds pca2v = pca2.get_principal_components();
      
      // calculate transformations
      for (int i = 0; i < 3; i++) {      
          pca1v[i]+=pca1.get_centroid();
          pca2v[i]+=pca2.get_centroid();	
      }  
      pca1v.push_back(pca1.get_centroid());
      pca2v.push_back(pca2.get_centroid());
      algebra::Transformation3D t1 = algebra::get_transformation_aligning_first_to_second(pca1v, pca2v);
      
      return t1;
}
*/



IMPBAYESIANEM_END_NAMESPACE

#endif /* IMPBAYESIANEM_UTILITIES_H */
