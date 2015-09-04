#ifndef RCPPSVAHCLUST_H
#define RCPPSVAHCLUST_H

#include <RcppCommon.h>
#include <RcppArmadillo.h>

#define EIGEN_MATRIXBASE_PLUGIN <RcppsvaEigenMatrixPlugin.h>
#define EIGEN_ARRAYBASE_PLUGIN <RcppsvaEigenArrayPlugin.h>
#include <RcppEigenForward.h>
#include <RcppsvaForward.h>

#include <RcppEigenWrap.h>
#include <RcppsvaEigenSugar.h>

#include <RcppsvaHClust/cluster.h>
#include <RcppsvaHClust/algorithm.h>
#include <RcppsvaHClust/method.h>
#include <RcppsvaHClust/hclust.h>

#endif
