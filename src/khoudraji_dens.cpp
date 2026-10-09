// Khoudraji-Archimedean copula densities (d = 2 or 3)
// Families: 1=Clayton, 2=Gumbel, 3=Frank, 4=Joe, 5=AMH
#include <Rcpp.h>
#include <cmath>
#include <limits>
using namespace Rcpp;

static const double EPS = 1e-10;
static const double DMIN = 1e-300;
static const double NINF = -std::numeric_limits<double>::infinity();

inline double clip01(double u) {
  if (u < EPS) return EPS;
  if (u > 1.0 - EPS) return 1.0 - EPS;
  return u;
}

inline int fam_code(const std::string& f) {
  if (f == "CLA" || f == "C") return 1;
  if (f == "GUM" || f == "G") return 2;
  if (f == "FRA" || f == "F") return 3;
  if (f == "JOE" || f == "J") return 4;
  if (f == "AMH" || f == "A") return 5;
  stop("Unknown family");
  return 0;
}

inline double phi1(double t, double th, int fam) {
  t = clip01(t);
  switch (fam) {
  case 1: return -th * std::pow(t, -th - 1.0);
  case 2: {
    double L = -std::log(t);
    return -th / t * std::pow(L, th - 1.0);
  }
  case 3: {
    if (std::fabs(th) < 1e-8) return -1.0 / t;
    double z = std::exp(-th * t);
    return th * z / (z - 1.0);
  }
  case 4: {
    double w = 1.0 - t;
    double D = 1.0 - std::pow(w, th);
    if (D < EPS) D = EPS;
    return -th * std::pow(w, th - 1.0) / D;
  }
  case 5: {
    double den = t * (1.0 - th * (1.0 - t));
    return (th - 1.0) / den;
  }
  }
  return 0.0;
}

inline double phi2(double t, double th, int fam) {
  t = clip01(t);
  switch (fam) {
  case 1: return th * (th + 1.0) * std::pow(t, -th - 2.0);
  case 2: {
    double L = -std::log(t);
    return th / (t * t) * std::pow(L, th - 1.0) +
           th * (th - 1.0) / (t * t) * std::pow(L, th - 2.0);
  }
  case 3: {
    if (std::fabs(th) < 1e-8) return 1.0 / (t * t);
    double z = std::exp(-th * t);
    double zm = z - 1.0;
    return th * th * z / (zm * zm);
  }
  case 4: {
    double w = 1.0 - t;
    double D = 1.0 - std::pow(w, th);
    if (D < EPS) D = EPS;
    return (th * (th - 1.0) * std::pow(w, th - 2.0) +
            th * std::pow(w, 2.0 * th - 2.0)) / (D * D);
  }
  case 5: {
    double a = t * (1.0 - th * (1.0 - t));
    return (1.0 - th) * (1.0 - th + 2.0 * th * t) / (a * a);
  }
  }
  return 0.0;
}

inline double phi3(double t, double th, int fam) {
  t = clip01(t);
  switch (fam) {
  case 1: return -th * (th + 1.0) * (th + 2.0) * std::pow(t, -th - 3.0);
  case 2: {
    double L = -std::log(t);
    double t3 = t * t * t;
    return -th * (th - 1.0) * (th - 2.0) * std::pow(L, th - 3.0) / t3 -
           3.0 * th * (th - 1.0) * std::pow(L, th - 2.0) / t3 -
           2.0 * th * std::pow(L, th - 1.0) / t3;
  }
  case 3: {
    if (std::fabs(th) < 1e-8) return -2.0 / (t * t * t);
    double z = std::exp(-th * t);
    double zm = z - 1.0;
    return th * th * th * z * (z + 1.0) / (zm * zm * zm);
  }
  case 4: {
    double w = 1.0 - t;
    double D = 1.0 - std::pow(w, th);
    if (D < EPS) D = EPS;
    double Np = -th * (th - 1.0) * (th - 2.0) * std::pow(w, th - 3.0) -
                th * (2.0 * th - 2.0) * std::pow(w, 2.0 * th - 3.0);
    double dDdt = th * std::pow(w, th - 1.0);
    double N = th * (th - 1.0) * std::pow(w, th - 2.0) +
               th * std::pow(w, 2.0 * th - 2.0);
    return (Np * D - 2.0 * N * dDdt) / (D * D * D);
  }
  case 5: {
    double h = 1e-6;
    return (phi2(t + h, th, fam) - phi2(t - h, th, fam)) / (2.0 * h);
  }
  }
  return 0.0;
}

inline double Carch(const double* u, int d, double th, int fam) {
  switch (fam) {
  case 1: {
    double s = -(d - 1);
    for (int j = 0; j < d; ++j) s += std::pow(clip01(u[j]), -th);
    if (s < EPS) s = EPS;
    return std::pow(s, -1.0 / th);
  }
  case 2: {
    double s = 0.0;
    for (int j = 0; j < d; ++j) s += std::pow(-std::log(clip01(u[j])), th);
    return std::exp(-std::pow(s, 1.0 / th));
  }
  case 3: {
    if (std::fabs(th) < 1e-8) {
      double p = 1.0;
      for (int j = 0; j < d; ++j) p *= clip01(u[j]);
      return p;
    }
    double num = 1.0;
    for (int j = 0; j < d; ++j) num *= (std::exp(-th * clip01(u[j])) - 1.0);
    double den = std::pow(std::exp(-th) - 1.0, d - 1.0);
    double arg = 1.0 + num / den;
    if (arg < EPS) arg = EPS;
    return -std::log(arg) / th;
  }
  case 4: {
    double p = 1.0;
    for (int j = 0; j < d; ++j) p *= (1.0 - std::pow(1.0 - clip01(u[j]), th));
    double s = 1.0 - p;
    if (s < 0.0) s = 0.0;
    return 1.0 - std::pow(s, 1.0 / th);
  }
  case 5: {
    if (d == 2) {
      double u1 = clip01(u[0]), u2 = clip01(u[1]);
      return u1 * u2 / (1.0 - th * (1.0 - u1) * (1.0 - u2));
    }
    double s = 0.0;
    for (int j = 0; j < d; ++j) {
      double uj = clip01(u[j]);
      s += std::log((1.0 - th * (1.0 - uj)) / uj);
    }
    double den = std::exp(s) - th;
    if (den < EPS) den = EPS;
    return (1.0 - th) / den;
  }
  }
  return 0.0;
}

inline double arch_dens_row(const double* u, int d, double th, int fam) {
  if ((fam == 3 && std::fabs(th) < 1e-8) ||
      (fam == 1 && std::fabs(th) < 1e-8) ||
      (fam == 2 && std::fabs(th - 1.0) < 1e-10))
    return 1.0;
  double Cv = Carch(u, d, th, fam);
  Cv = clip01(Cv);
  double p1 = phi1(Cv, th, fam);
  double p2 = phi2(Cv, th, fam);
  double d1prod = 1.0;
  for (int j = 0; j < d; ++j) d1prod *= phi1(u[j], th, fam);
  if (d == 2) {
    return -d1prod * p2 / (p1 * p1 * p1);
  }
  double p3 = phi3(Cv, th, fam);
  double g3 = -(p3 * p1 - 3.0 * p2 * p2) / (p1 * p1 * p1 * p1 * p1);
  return g3 * d1prod;
}

inline double khoud2_row(double u, double v, double th,
                         double a, double b, int fam) {
  u = clip01(u); v = clip01(v);
  double x = clip01(std::pow(u, a));
  double y = clip01(std::pow(v, b));
  double uv[2] = {x, y};
  double Cv = clip01(Carch(uv, 2, th, fam));
  double cpsi = arch_dens_row(uv, 2, th, fam);
  if (cpsi < 0.0) cpsi = 0.0;
  double C1 = phi1(x, th, fam) / phi1(Cv, th, fam);
  double C2 = phi1(y, th, fam) / phi1(Cv, th, fam);
  double t1 = a * b * cpsi;
  // The power factors are in the original coordinates.  Using x^{-a}
  // and y^{-b} here would incorrectly produce u^{-a^2} and v^{-b^2}.
  double t2 = a * (1.0 - b) * std::pow(v, -b) * C1;
  double t3 = b * (1.0 - a) * std::pow(u, -a) * C2;
  double t4 = (1.0 - a) * (1.0 - b) *
              std::pow(u, -a) * std::pow(v, -b) * Cv;
  double dens = t1 + t2 + t3 + t4;
  if (!(dens > 0.0) || !std::isfinite(dens)) return 0.0;
  return dens;
}

inline double logadd(double x, double y) {
  if (x == NINF) return y;
  if (y == NINF) return x;
  double m = x > y ? x : y;
  return m + std::log(std::exp(x - m) + std::exp(y - m));
}

inline double log_khoud2_row(double u, double v, double th,
                             double a, double b, int fam) {
  u = clip01(u); v = clip01(v);
  double x = clip01(std::pow(u, a));
  double y = clip01(std::pow(v, b));
  double uv[2] = {x, y};
  double Cv = clip01(Carch(uv, 2, th, fam));
  double cpsi = arch_dens_row(uv, 2, th, fam);
  double C1 = phi1(x, th, fam) / phi1(Cv, th, fam);
  double C2 = phi1(y, th, fam) / phi1(Cv, th, fam);
  double out = NINF;
  if (a > 0.0 && b > 0.0 && cpsi > 0.0)
    out = logadd(out, std::log(a) + std::log(b) + std::log(cpsi));
  if (a > 0.0 && b < 1.0 && C1 > 0.0)
    out = logadd(out, std::log(a) + std::log1p(-b) -
                 b * std::log(v) + std::log(C1));
  if (b > 0.0 && a < 1.0 && C2 > 0.0)
    out = logadd(out, std::log(b) + std::log1p(-a) -
                 a * std::log(u) + std::log(C2));
  if (a < 1.0 && b < 1.0 && Cv > 0.0)
    out = logadd(out, std::log1p(-a) + std::log1p(-b) -
                 a * std::log(u) - b * std::log(v) + std::log(Cv));
  return std::isfinite(out) ? out : std::log(DMIN);
}

inline double khoud3_row(const double* u_in, double th,
                         const double* alpha, int fam) {
  double u[3], x[3];
  for (int j = 0; j < 3; ++j) {
    u[j] = clip01(u_in[j]);
    x[j] = clip01(std::pow(u[j], alpha[j]));
  }
  double Cv = clip01(Carch(x, 3, th, fam));
  double p1C = phi1(Cv, th, fam);
  double p2C = phi2(Cv, th, fam);
  double p3C = phi3(Cv, th, fam);
  double g0 = Cv;
  double g1 = 1.0 / p1C;
  double g2 = -p2C / (p1C * p1C * p1C);
  double g3 = -(p3C * p1C - 3.0 * p2C * p2C) /
              (p1C * p1C * p1C * p1C * p1C);
  double phi1x[3];
  for (int j = 0; j < 3; ++j) phi1x[j] = phi1(x[j], th, fam);

  double dens = 0.0;
  for (int mask = 0; mask < 8; ++mask) {
    int inS[3];
    int k = 0;
    for (int j = 0; j < 3; ++j) {
      inS[j] = (mask >> j) & 1;
      k += inS[j];
    }
    double AS;
    if (k == 0) {
      AS = g0;
    } else if (k == 1) {
      int j = inS[0] ? 0 : (inS[1] ? 1 : 2);
      AS = g1 * phi1x[j];
    } else if (k == 2) {
      int j1 = -1, j2 = -1;
      for (int j = 0; j < 3; ++j) if (inS[j]) {
        if (j1 < 0) j1 = j; else j2 = j;
      }
      AS = g2 * phi1x[j1] * phi1x[j2];
    } else {
      AS = g3 * phi1x[0] * phi1x[1] * phi1x[2];
    }
    double fac = 1.0;
    for (int j = 0; j < 3; ++j) {
      if (inS[j]) fac *= alpha[j];
      else {
        fac *= (1.0 - alpha[j]);
        AS *= std::pow(u[j], -alpha[j]);
      }
    }
    dens += fac * AS;
  }
  if (!(dens > 0.0) || !std::isfinite(dens)) return 0.0;
  return dens;
}

inline double log_khoud3_row(const double* u_in, double th,
                             const double* alpha, int fam) {
  double u[3], x[3];
  for (int j = 0; j < 3; ++j) {
    u[j] = clip01(u_in[j]);
    x[j] = clip01(std::pow(u[j], alpha[j]));
  }
  double Cv = clip01(Carch(x, 3, th, fam));
  double p1C = phi1(Cv, th, fam);
  double p2C = phi2(Cv, th, fam);
  double p3C = phi3(Cv, th, fam);
  double g[4] = {
    Cv,
    1.0 / p1C,
    -p2C / (p1C * p1C * p1C),
    -(p3C * p1C - 3.0 * p2C * p2C) /
      (p1C * p1C * p1C * p1C * p1C)
  };
  double phi1x[3];
  for (int j = 0; j < 3; ++j) phi1x[j] = phi1(x[j], th, fam);

  double out = NINF;
  for (int mask = 0; mask < 8; ++mask) {
    int k = 0;
    double partial = 1.0;
    double logfac = 0.0;
    bool valid = true;
    for (int j = 0; j < 3; ++j) {
      bool inS = ((mask >> j) & 1) != 0;
      if (inS) {
        ++k;
        partial *= phi1x[j];
        if (alpha[j] <= 0.0) valid = false;
        else logfac += std::log(alpha[j]);
      } else {
        if (alpha[j] >= 1.0) valid = false;
        else logfac += std::log1p(-alpha[j]) -
                       alpha[j] * std::log(u[j]);
      }
    }
    partial *= g[k];
    if (valid && partial > 0.0 && std::isfinite(partial))
      out = logadd(out, logfac + std::log(partial));
  }
  return std::isfinite(out) ? out : std::log(DMIN);
}

// [[Rcpp::export]]
NumericVector d_khoud_cpp(NumericMatrix U, double th, NumericVector alpha,
                          std::string family) {
  int n = U.nrow();
  int d = U.ncol();
  int fam = fam_code(family);
  NumericVector out(n);
  if (d == 2) {
    double a = alpha[0], b = alpha[1];
    for (int i = 0; i < n; ++i)
      out[i] = khoud2_row(U(i, 0), U(i, 1), th, a, b, fam);
  } else if (d == 3) {
    double al[3] = {alpha[0], alpha[1], alpha[2]};
    double u[3];
    for (int i = 0; i < n; ++i) {
      u[0] = U(i, 0); u[1] = U(i, 1); u[2] = U(i, 2);
      out[i] = khoud3_row(u, th, al, fam);
    }
  } else {
    stop("only d=2 or d=3 implemented in C++");
  }
  return out;
}

// [[Rcpp::export]]
NumericVector log_d_khoud_cpp(NumericMatrix U, double th, NumericVector alpha,
                              std::string family) {
  int n = U.nrow();
  int d = U.ncol();
  int fam = fam_code(family);
  NumericVector out(n);
  if (d == 2) {
    for (int i = 0; i < n; ++i)
      out[i] = log_khoud2_row(U(i, 0), U(i, 1), th,
                              alpha[0], alpha[1], fam);
  } else if (d == 3) {
    double al[3] = {alpha[0], alpha[1], alpha[2]};
    double u[3];
    for (int i = 0; i < n; ++i) {
      u[0] = U(i, 0); u[1] = U(i, 1); u[2] = U(i, 2);
      out[i] = log_khoud3_row(u, th, al, fam);
    }
  } else {
    stop("only d=2 or d=3 implemented in C++");
  }
  return out;
}

// [[Rcpp::export]]
NumericVector arch_density_cpp(NumericMatrix U, double th, std::string family) {
  int n = U.nrow();
  int d = U.ncol();
  int fam = fam_code(family);
  NumericVector out(n);
  std::vector<double> u(d);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < d; ++j) u[j] = U(i, j);
    double dens = arch_dens_row(u.data(), d, th, fam);
    out[i] = dens > 0.0 ? dens : 0.0;
  }
  return out;
}
