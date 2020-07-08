#include <iostream>
#include <pplp_polyhedron.h>
#include <xmmintrin.h>
//#include <raytracing/polyhedron.h>
#include "genPoly.h"

GenPoly::GenPoly ( char* para[] ) {
  std::string path(para[1]) ;
  _poly_num = std::stoi( para[2] ) ;
  _total_cons_num = std::stoi( para[3] ) ;
  _re_cons_num = std::stoi( para[4] ) ;   
  _vari_num = std::stoi( para[5] ) ;
  _zero_num = std::stoi( para[6] ) ;
  _distance = std::stoi( para[7] ) ;
  _ask_closed = std::stoi( para[8] ) ;
  _cons_num = _total_cons_num - _re_cons_num ;
  if(_zero_num >= _vari_num) {
    std::cerr << "Error: the number of zeros is larger than the number of variables." << std::endl ;
    std::terminate() ;
  }
  _ofs1.open(path, std::ofstream::out | std::ofstream::app) ;
  if ( ! _ofs1.is_open() ) {
    std::cerr << "Cannot open file." << std::endl ;
    std::terminate() ; 
  }
  _constraints.resize(_cons_num, _vari_num) ;
  _constant.resize(_cons_num) ;
  _reconstraints.resize(_re_cons_num, _vari_num) ;
  _reconstant.resize(_re_cons_num) ;
  _center.resize(_vari_num) ;
  for(int i = 0; i < _vari_num; ++ i) {
    _center(i) = 0.0 ;
  }
}

GenPoly::~GenPoly() {
  _ofs1.close() ;
}

void GenPoly::WriteInfo () {
  std::string line ;
  line = "constraints = " + std::to_string(_total_cons_num) ;
  _ofs1 << line << std::endl ;
  line = "redundancy = " + std::to_string(_re_cons_num) ;
  _ofs1 << line << std::endl ;
  line = "variables = " + std::to_string(_vari_num) ;
  _ofs1 << line << std::endl ; 
  line = "zeros = " + std::to_string(_zero_num) ;
  _ofs1 << line << std::endl ;  
  _ofs1 << std:: endl ;
  _ofs1 << "Begin" << std::endl << std:: endl ;
}
// Generate the poly outside firstly, the distance is 5
bool GenPoly::CreateIrreCons () {
  std::random_device rd ;
  std::mt19937_64 gen( rd() ) ;
  // create random direction and distance
  // the constraints are Cx <= b
  std::uniform_int_distribution<> disC(-50, 50) ;
  std::uniform_int_distribution<> disIdx(0, _vari_num-1) ;
  for(int i = 0; i < _cons_num; ++ i) {
    int j = 0 ;
    // generate C
    while(j < _vari_num) {
      int coef = disC(gen) ;
      // do not generate zero now
      if (coef == 0) continue ;
      else {
        _constraints(i, j) = coef ;
        ++ j ;
      }
    }
  } 
  // set zeros according to the density
  if (_zero_num != 0) {
    for (int i = 0; i < _cons_num; ++ i) {
      int z = 0 ;
      std::vector<bool> isZero(_vari_num, false) ;
      while (z < _zero_num) { 
        int zeroIdx = disIdx(gen) ;  
        if (isZero[zeroIdx] == true) continue ;
        else {
          _constraints(i, zeroIdx) = 0 ;
          isZero[zeroIdx] = true ;
          ++ z ; 
        }
      }  
    }
  }
  for (int i = 0; i < _cons_num; ++ i) {
    // b is fixed, now assume Cx <= b
    float curr = _constraints.row(i).norm() * _distance +
        _constraints.row(i) * _center ;
    // get integer
    // the _distance of two poly should be larger than 1
    _constant(i) = std::floor(curr) ;
  }
  for (int i = 0; i < _cons_num; ++ i) {
    // test with the point
    double temp = _center.dot( _constraints.row(i) ) ;  
    if(temp - _constant(i) > 0) {
      double curr = _constraints.row(i) * _center -
          _constraints.row(i).norm() * _distance ;
      _constant(i) = std::ceil(curr) ;
    }
  }
  return Minimize() ;    
}

bool GenPoly::Minimize () {
  // minimize the polyhedron
  std::random_device rd ;
  std::mt19937_64 gen( rd() ) ;
  std::uniform_int_distribution<> disC(-50, 50) ;
  std::uniform_int_distribution<> disB(0, 20) ;
  std::uniform_int_distribution<> disIdx(0, _vari_num-1) ;
  // set a bound for it
  int count = 0 ; 
  PPLP::Polyhedron poly(_cons_num, _vari_num) ;
  //Polyhedron poly(_cons_num, _vari_num) ;
  for (int i = 0; i < _cons_num; ++ i) {
    poly.SetConstraint( i, _constraints.row(i) ) ;
    poly.SetConstant( i, _constant(i) ) ;
  }
  while (count < BOUND) {
    ++ count ;
    if (poly.Minimize() == false) {
      return false ;
    }
    std::vector<int> reIdx = poly.GetInactiveIdx() ; 
    if ( (int)reIdx.size() == 0 ) break ;
    for(int i = 0; i < (int)reIdx.size(); ++ i) {
      int j = 0 ;
      // generate C
      while(j < _vari_num) {
        int coef = disC(gen) ;
        // do not generate zero now
        if (coef == 0) continue ;
        else {
          poly.SetCoef(reIdx[i], j, coef) ;
          ++ j ;
        }
      }
    } 
    // set zeros according to the density
    if (_zero_num != 0) {
      for (int i = 0; i < (int)reIdx.size(); ++ i) {
        int z = 0 ;
        std::vector<bool> isZero(_vari_num, false) ;
        while (z < _zero_num) { 
          int zeroIdx = disIdx(gen) ;  
          if (isZero[zeroIdx] == true) continue ;
          else {
            poly.SetCoef(reIdx[i], zeroIdx, 0) ;
            isZero[zeroIdx] = true ;
            ++ z ; 
          }
        }  
      }
    }
    for (int i = 0; i < (int)reIdx.size(); ++i) {
      // b is fixed
      int currIdx = reIdx[i] ;
      double curr = poly.get_coefficients().row(currIdx).norm() * _distance +
        poly.get_coefficients().row(currIdx) * _center ;
        poly.SetConstant( currIdx, std::floor(curr) ) ;
    }
    for (int i = 0; i < (int)reIdx.size(); ++ i) {
      int currIdx = reIdx[i] ;
      // test with the point
      double temp = _center.dot( poly.get_coefficients().row( reIdx[i] ) ) ;  
      if(temp - poly.GetConstant( reIdx[i] ) > 0) {
        double curr = poly.get_coefficients().row(currIdx) * _center -
            poly.get_coefficients().row(currIdx).norm() * _distance ;
        poly.SetConstant( currIdx, std::ceil(curr) ) ;
      }
    }
    poly.Init() ; 
  }   
  bool result = true ;
  // Cannot generate a polyhedron with current parameters, ignore it
  if (count == BOUND) {
    result = false ;
  }
  // get the new constraints and constants
  _constraints = poly.get_coefficients() ;
  _constant = poly.get_constants() ;

  if (_ask_closed == true && poly.IsOpen() == true) {
    result = false ; 
  }
   
  return result ;
}


void GenPoly::CreateReCons () {
  std::random_device rd ;
  std::mt19937_64 gen( rd() ) ;
  // create the redundant constraints
  std::uniform_int_distribution<> disLamda(1, 5) ;
  std::uniform_int_distribution<> disCons(0, _cons_num-1) ;
  std::uniform_int_distribution<> disRand(1, 50) ;
  Matrix constraints2(_re_cons_num, _vari_num) ;
  ColVector constant2(_re_cons_num) ;
  for(int i = 0; i < _re_cons_num; ++ i) {
    int lamda1 = disLamda(gen) ;
    int lamda2 = disLamda(gen) ;
    int lamda3 = disLamda(gen) ;
    int idx1 = disCons(gen) ;
    int idx2 = disCons(gen) ;
    int idx3 = disCons(gen) ;
    int randConstant = disRand(gen) ;
    _reconstraints.row(i) = lamda1 * _constraints.row(idx1) 
      + lamda2 * _constraints.row(idx2) + lamda3 * _constraints.row(idx3) ; 
    _reconstant(i) = lamda1 * _constant(idx1) 
      + lamda2 * _constant(idx2) + lamda3 * _constant(idx3) + randConstant ;
  }
}

std::vector<int> GenPoly::MixConstraints () {
  // mix the redundant and irredundant constraints 
  std::random_device rd ;
  std::mt19937_64 gen( rd() ) ;
  std::vector<int> consOrder ;
  for (int i = 0; i < _cons_num; ++ i) {
    consOrder.push_back(i) ;
  } 
  std::vector<int>::iterator it ; 
  for (int i = 0; i < _re_cons_num; ++ i) {
    it = consOrder.begin() ;
    std::uniform_int_distribution<> disOrdIdx(0, _cons_num + i) ;
    int currOrdIdx = disOrdIdx(gen) ;
    if (currOrdIdx == _cons_num + i) {
      consOrder.push_back(_cons_num + i) ;
    }
    else {
      consOrder.insert(it + currOrdIdx, _cons_num + i) ;
    } 
  }
  return consOrder ;
}


void GenPoly::WritePoly (const std::vector<int>& consOrder, int polyIdx) {
  // write into files
  // write polyhedra name
  std::string name = "P_" + std::to_string(_total_cons_num) + "_" 
      + std::to_string(_re_cons_num) 
      + "_" + std::to_string(_vari_num) + "_" + std::to_string(_zero_num) 
      + "_0_" + std::to_string(polyIdx) ;
  _ofs1 << name << std::endl ;
  for (int i = 0; i < _total_cons_num; ++ i) {
    int currIdx = consOrder[i] ;
    if (currIdx < _cons_num) {
      for (int j = 0; j < _vari_num; ++ j) {
        _ofs1 << _constraints(currIdx, j) << " " ;
      }
      _ofs1 << "<=" << " " << _constant(currIdx) << std::endl ;
    }
    else {
      for (int j = 0; j < _vari_num; ++ j) {
        _ofs1 << _reconstraints(currIdx - _cons_num, j) << " " ;
      }
      // TODO consider different operators later
      _ofs1 << "<=" << " " << _reconstant(currIdx - _cons_num) << std::endl ;
    }
  }
  _ofs1 << std::endl ;
}
