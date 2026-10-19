#ifndef BRIDGE_REAL_H
#define BRIDGE_REAL_H

static inline double to_double(double val) {
  return val;
}

static inline double to_double(float val) {
  double tmp=val;
  return tmp;
}

static inline float to_float(double val) {
  float tmp=val;
  return val;
}

static inline float to_float(float val) {
  return val;
}

#endif
