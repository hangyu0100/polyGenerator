#include <iostream>
#include <string>
#include <vector>
#include <exception>
#include <cstdio>
#ifndef _GENERATOR
#define _GENERATOR
#include "genPoly.h"
#endif

#define POLYBOUND 100

/*******************************************************************************
 * This function generate polyhedra for test
 * TODO consider the operator: <, <=, or =
*******************************************************************************/
int main (int argc, char* argv[]) {
  if(argc != 9) {
    std::cout << "usage: " << argv[0] << " <filename> <number_of_polyhedra> " 
      <<  "<total_number_of_constraints> <number_of_redundant_constraints> "
      << "<number_of_variables> <density (number of zeros) > <distance> " 
      << "<ask_closed_poly>"<< std::endl ; 
    return -1 ;
  }     
  GenPoly generator(argv) ; 
  generator.WriteInfo() ; 
  int count = 0, k = 0 ;
  while (k < generator.get_poly_num()) {
    if (count == POLYBOUND) {
      std::remove(argv[1]) ; 
      std::cerr << "Cannot generate polyhedra P_" << argv[3] << "_" 
        << argv[4] << "_" << argv[5] << "_" << argv[6] << std::endl ; 
      break ;
    }
    if (generator.CreateIrreCons() == false) {
      ++ count ;
      continue ;
    }
    generator.CreateReCons() ;
    std::vector<int> order = generator.MixConstraints() ;
    generator.WritePoly(order, k) ; 
    ++ k ;
  }
  return 0 ;
}
