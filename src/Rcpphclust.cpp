#include <algorithm>
#include <functional>
#include <stack>
#include <stdexcept>
#include <omp.h>

#include <RcppsvaHClust.h>

// Rcpp extension function
namespace Rcpp {

template <> Eigen::RowMajorNumericMatrix as(SEXP x) {
  return Eigen::RowMajorNumericMatrix(as<Eigen::MapNumericMatrix>(x));
}

template <> RcppsvaHClust::LinkageKinds as(SEXP x){
  switch (as<int>(x)) {
  default: throw not_compatible("Linkage method invalid or not yet supported"); 
  case 1: return RcppsvaHClust::WARD;
  case 2: return RcppsvaHClust::AVERAGE;
  case 3: return RcppsvaHClust::SINGLE;
  case 4: return RcppsvaHClust::COMPLETE;
  }
}

template <> RcppsvaHClust::DistanceKinds as(SEXP x){
  switch (as<int>(x)) {
  default: throw not_compatible("Distance method invalid or not yet supported"); 
  case 1: return RcppsvaHClust::EUCLIDEAN;
  case 2: return RcppsvaHClust::MANHATTAN;
  case 3: return RcppsvaHClust::MAXIMUM;
  case 4: return RcppsvaHClust::MINKOWSKI;
  case 5: return RcppsvaHClust::SCOSINE;
  case 6: return RcppsvaHClust::UCOSINE;
  }
}

template <> SEXP wrap( const RcppsvaHClust::Hclust& hclust ) {
  return List::create( _["merge"] = hclust.merge, _["height"] = hclust.height, _["order"] = hclust.order ); 
}

} // Rcpp

RcppExport SEXP Rcppsva_hclust_from_data(SEXP data, SEXP link, SEXP dist, SEXP minkowski) {
  BEGIN_RCPP
  
  Eigen::RowMajorNumericMatrix data_e(Rcpp::as<Eigen::RowMajorNumericMatrix>(data));
  
  RcppsvaHClust::LinkageKinds  lk = Rcpp::as<RcppsvaHClust::LinkageKinds>(link);
  RcppsvaHClust::DistanceKinds dk = Rcpp::as<RcppsvaHClust::DistanceKinds>(dist);
  
  switch (lk) {
  default: 
    throw std::invalid_argument("Linkage or distance method not yet supported");
  case RcppsvaHClust::WARD: {
    typedef RcppsvaHClust::NumericCluster::center cluster_type;
    
    RcppsvaHClust::ClusterVector<cluster_type> clusters(data_e.rows());	
    RcppsvaHClust::init_clusters_from_rows(data_e, clusters);
    
    RcppsvaHClust::cluster_via_rnn( RcppsvaHClust::wards_link<cluster_type>(), clusters );
    
    return Rcpp::wrap(clusters);	
  }
  case RcppsvaHClust::AVERAGE: {
    typedef RcppsvaHClust::NumericCluster::obs cluster_type;
    
    RcppsvaHClust::ClusterVector<cluster_type> clusters(data_e.rows());
    RcppsvaHClust::init_clusters_from_rows(data_e, clusters);
    
    RcppsvaHClust::cluster_via_rnn( RcppsvaHClust::average_link<cluster_type>( RcppsvaHClust::stored_data_rows(data_e, dk, Rcpp::as<double>(minkowski)) ), clusters );
    
    return Rcpp::wrap(clusters);
  }
  case RcppsvaHClust::SINGLE: {
    typedef RcppsvaHClust::NumericCluster::plain cluster_type;
    
    RcppsvaHClust::ClusterVector<cluster_type> clusters(data_e.rows());
    RcppsvaHClust::init_clusters_from_rows(data_e, clusters);
    
    RcppsvaHClust::cluster_via_slink( RcppsvaHClust::stored_data_rows(data_e, dk, Rcpp::as<double>(minkowski)), clusters );
    
    return Rcpp::wrap(clusters);
  }
  case RcppsvaHClust::COMPLETE: {
    typedef RcppsvaHClust::NumericCluster::obs cluster_type;
    
    RcppsvaHClust::ClusterVector<cluster_type> clusters(data_e.rows());
    RcppsvaHClust::init_clusters_from_rows(data_e, clusters);
    
    RcppsvaHClust::cluster_via_rnn( RcppsvaHClust::complete_link<cluster_type>( RcppsvaHClust::stored_data_rows(data_e, dk, Rcpp::as<double>(minkowski)) ), clusters );
    
    return Rcpp::wrap(clusters);
  }
    
  }
  
  END_RCPP
}
