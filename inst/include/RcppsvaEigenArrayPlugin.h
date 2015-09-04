#ifndef RCPPSVAEIGENARRAYPLUGIN_H
#define RCPPSVAEIGENARRAYPLUGIN_H

RealScalar squaredNorm() const {
  return square().sum();
}

RealScalar norm() const {
  return sqrt(square().sum());
}

#endif
