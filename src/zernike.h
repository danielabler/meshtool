#ifndef ZERNIKE_H_
#define ZERNIKE_H_

#include <unistd.h>
#include <iostream>
#include <cmath>
#include <map>


namespace zernike {

  typedef std::map<int, double> j_coeff_map_type;

  struct zernike_surface {
    j_coeff_map_type zernike_coefficients;
    double pupil_radius ;
  };

  struct zernike_surface_return {
    double z_at_0_0;
    double inout;
    double z_max;
    double z_zernike;
  };


  struct RetValZ {
    float z_zernike;
    float z_max;
    float r_max;
    float z_at_0_0;
  };

  struct cornea {
    zernike_surface surface_anterior;
    zernike_surface surface_posterior;
    zernike_surface lenticule_surface_anterior;
    zernike_surface lenticule_surface_posterior;
    double pupil_radius ;
    double surface_distance;
  };

  struct coord_2D_cart {
    double x;
    double y;
  };

  struct coord_2D_spher {
    double r;
    double phi;
  };

  struct coord_3D_cart {
    double x;
    double y;
    double z;
  };

  struct zernike_coeff {
    int n;
    int m;
  };

  coord_2D_spher cart_to_spher_2D(coord_2D_cart);
  coord_2D_cart spher_to_cart_2D(coord_2D_spher);
  unsigned int factorial(unsigned int);
  int zernike_nm_to_j(zernike_coeff);
  zernike_coeff zernike_j_to_nm(int);
  RetValZ compute_surface(double, double, double, j_coeff_map_type, int, double);
  RetValZ compute_surface_shifted(double , double , double , j_coeff_map_type , int, double, double , double );
  double evaluate_zernike_surface(j_coeff_map_type, double, double, int); 
  double evaluate_zernike_surface_cart(j_coeff_map_type, coord_2D_cart, double, int );
  double zernike_term(int, int, double, double);
  zernike_surface_return create_zernike_surface(double, double, double, zernike_surface, int , double, double, double );

    }  // end zernike

#endif // ZERNIKE_H_
