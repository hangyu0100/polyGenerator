#ifndef _GENPOLY_GENPOLY
#define _GENPOLY_GENPOLY

#include <eigen3/Eigen/Dense>
#include <fstream>
#include <random>

typedef Eigen::MatrixXd Matrix ;
typedef Eigen::VectorXd ColVector ;
#define BOUND 100

class GenPoly {
public:
  GenPoly ( char* para[] ) ;
  ~GenPoly () ;
  // use it after dividing coefficients by gcd
  void WriteInfo () ;
  bool CreateIrreCons () ;
  void CreateReCons () ;
  std::vector<int> MixConstraints () ;
  void WritePoly (const std::vector<int>& consOrder, int polyIdx) ;
  int get_poly_num () {
    return _poly_num ;
  }
private:
  bool Minimize () ;
  int _poly_num ;
  int _total_cons_num ;
  int _cons_num ;
  int _re_cons_num ;
  int _vari_num ;
  int _zero_num ;
  int _distance ;
  bool _ask_closed ;
  std::ofstream _ofs1 ; 
  Matrix _constraints ;
  ColVector _constant ; 
  Matrix _reconstraints ;
  ColVector _reconstant ; 
  ColVector _center ;
} ;

#endif
