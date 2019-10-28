#include <zernike.h>
// LOGGING
#include "easylogging++.h"


namespace zernike {
  
typedef std::map<int, double> j_coeff_map_type;

  // OTHER FUNCTIONS

  unsigned int factorial(unsigned int number) {
    if (number <= 1) {
      return 1; 
    }
    return factorial(number - 1) * number;
  }

  coord_2D_spher cart_to_spher_2D(coord_2D_cart coords){
    coord_2D_spher coords_spher;
    coords_spher.r   =  sqrt(coords.x * coords.x + coords.y * coords.y);
    if (coords_spher.r != 0.0){
      coords_spher.phi =  std::atan2(coords.y , coords.x);
    }
    else{
      coords_spher.phi =  0;
    }
    return coords_spher;
  }


  coord_2D_cart spher_to_cart_2D(coord_2D_spher coords){
    coord_2D_cart coords_cart;
    coords_cart.x = coords.r * std::cos(coords.phi);
    coords_cart.y = coords.r * std::sin(coords.phi);
    return coords_cart;
  }


  // ZERNIKE COEFFICIENT TRANSFORMATION

zernike_coeff zernike_j_to_nm(int j){
  zernike_coeff coeffs;
  coeffs.n = floor( sqrt(2*j+1) + 0.5 ) - 1;
  coeffs.m = 2*j - coeffs.n * (coeffs.n+2);
  return coeffs;
}

int zernike_nm_to_j(zernike_coeff nm){
  int n = nm.n;
  int m = nm.m;
  int j = (n*n + 2*n + m) / 2. ;
  return j;
}

  // ZERNIKE CONSTRUCTION


double zernike_radial_term(int n, int m, double rho){
  double radial_term = 0;

  if (n>=abs(m)) {
    int diff = n - abs(m);
    LOG(DEBUG) << "diff = " << diff ;
    if ( diff % 2 ==0 ) /* x is even */ {
      for (int k=0; k<= (diff / 2); k++){
        int krondelta = 0;
        if (m==0){
          krondelta = 1;
        }
        double normalisation = sqrt( 2*(n+1) / (1+krondelta) );
        double radial_term_iteration_1 = pow(-1,k) * factorial(n - k) * pow(rho, (n-2*k));
        double radial_term_iteration_2 = factorial(k) * factorial((n+abs(m))/2-k) * factorial((n-abs(m))/2-k);
        radial_term = radial_term + normalisation * radial_term_iteration_1 / radial_term_iteration_2;
        LOG(DEBUG) << "k = " << k << ", radial term = " << radial_term ;
      }
      
    }
  }
  LOG(DEBUG) << "zernike radial term: " << n << " , " << m << " -> " << radial_term ;
  return radial_term;
}


double zernike_azimuthal_term(int m, double phi){
  double azimuthal_term = 0;
  if (m>0) {
    azimuthal_term =   std::cos( std::abs(m) * phi);
  }
  else if (m<0) {
    azimuthal_term =   std::sin( std::abs(m) * phi);
  }
  else{
    azimuthal_term = 1;
  }
  LOG(DEBUG) << "zernike azimuth term: " << m << " -> " << azimuthal_term;
  return azimuthal_term;
}


double zernike_term(int n, int m, double rho, double phi){
  return zernike_radial_term(n, m, rho) * zernike_azimuthal_term(m, phi);
}


  // SURFACE COMPUTATION

RetValZ compute_surface(double x, double y, double z, j_coeff_map_type j_coeff_map, int j_max, double pupil_radius){  
  // precompute z_max and z_at 0 0  for each surface
  coord_2D_cart coord_2D;
  coord_2D.x = x;
  coord_2D.y = y;
  double z_zernike = evaluate_zernike_surface_cart(j_coeff_map, coord_2D, pupil_radius, j_max);
  
  coord_2D_cart coord_2D_max;
  coord_2D_max.x = 0.0;
  coord_2D_max.y = pupil_radius;
  double z_max      = evaluate_zernike_surface_cart(j_coeff_map, coord_2D_max, pupil_radius, j_max);
  
  double z_at_0_0 = evaluate_zernike_surface(j_coeff_map, 0.0, 10.0, j_max); 
  
  // x,y within pupil opening
  double max_r = sqrt(x * x + y * y);
  
  LOG(DEBUG) << "x = " << x << ", y = " << y << ", z = " << z << " => z_zernike = " << z_zernike << ", z_max = " << z_max << ", z_at_0_0 = " << z_at_0_0;
  
  // 1) move apex to z=0 
  double z_max_shifted, z_zernike_shifted;
  z_max_shifted 	    = z_max 	    - z_at_0_0;
  z_zernike_shifted 	= z_zernike 	- z_at_0_0;
  // 2) Orientation (opening towards positive z)
  RetValZ result;
  result.r_max 	    = max_r;
  result.z_at_0_0   = z_at_0_0; 
  //Orientation:
  //  (1) surface opened towards positive z  <=> j_4 > 0
  //  (2) surface opened towards negative z  <=> j_4 < 0 
  if (j_coeff_map[4] < 0){ // flip orientation
    result.z_zernike = - z_zernike_shifted;
    result.z_max     = - z_max_shifted; 
  }
  else{			   // leave orientation as is	
    result.z_zernike = z_zernike_shifted;
    result.z_max     = z_max_shifted; 
  }
  return result;
}

  RetValZ compute_surface_shifted(double x, double y, double z, j_coeff_map_type j_coeff_map, int j_max, double pupil_radius, double offset_from_z_0, double shift){
    RetValZ result = compute_surface(x, y, z, j_coeff_map, j_max, pupil_radius);
    result.z_zernike = result.z_zernike + shift + offset_from_z_0;
    result.z_max     = result.z_max     + shift + offset_from_z_0;
    return result;
  }


double evaluate_zernike_surface(j_coeff_map_type j_coeff_map, double rho, double phi, int j_coeff_max){
  double Z_from_zernike = 0.0;
  for (int j=0; j<= j_coeff_max; j++)
  {
    double coeff = 0.0;
    LOG(DEBUG) << "j = " << j << "; j_coeff_map.count(j) = " << j_coeff_map.count(j);
    if (j_coeff_map.count(j)) {
      coeff = j_coeff_map[j];
      LOG(DEBUG) << "j = " << j << ", coeff = " << coeff;
      }
    zernike_coeff nm_coeff = zernike_j_to_nm(j);
    LOG(DEBUG) << "j = " << j << " -> n = " << nm_coeff.n << " , m = " << nm_coeff.m ;
    double zernike = zernike_term(nm_coeff.n, nm_coeff.m, rho, phi);
    LOG(DEBUG) << "j = " << j << " : zernike term = " << zernike << " , coeff = " << coeff ;
    Z_from_zernike = Z_from_zernike + coeff * zernike;
    LOG(DEBUG) << "z_from_zernike = " << Z_from_zernike;
  }
  return Z_from_zernike;
}


double evaluate_zernike_surface_cart(j_coeff_map_type j_coeff_map, coord_2D_cart coord_cart, double pupil_radius, int j_coeff_max){
  /* Convert cartesian 2D coordinate to spherical coordinate */
  coord_2D_spher coord_spher = cart_to_spher_2D(coord_cart);
  /* TEST: output and back transformation */
  LOG(DEBUG) << "Cartesion: " << coord_cart.x << " , " << coord_cart.y;
  LOG(DEBUG) << "Spherical  " << coord_spher.r << " , " << coord_spher.phi;
  coord_2D_cart coord_cart_2 = spher_to_cart_2D(coord_spher);
  LOG(DEBUG) << "Cartesion (from spherical): " << coord_cart_2.x << " , " << coord_cart_2.y ;
  /* Compute rho = normalised polar radius */
  double rho = coord_spher.r/pupil_radius;
  /* Evaluate Z */
  double Z =  evaluate_zernike_surface(j_coeff_map, rho, coord_spher.phi, j_coeff_max);
  LOG(DEBUG) << "Z(" << "rho="<< rho <<", phi=" << coord_spher.phi << ") = " << Z ;
  return Z;
}


  zernike_surface_return create_zernike_surface(double x, double y, double z, zernike_surface surface, int max_number_zernike, double offset_from_z_0, double shift, double surface_thickness){
    /*
     * double offset_from_z_0:    arbitrary offset away from z=0; apex at z=0 seems to cause meshing issues
     * double shift:              shift -> difference between shift anterior and shift posterior is true thickness of lense
     * double surface_thickness:  this is the 'thickness',
     *                            i.e. max height of half a lense; should be larger than the shift difference betwee anterior and posterior
    */

    RetValZ result = compute_surface_shifted(x, y, z, surface.zernike_coefficients, max_number_zernike, surface.pupil_radius, offset_from_z_0, shift);
    double diff_zernike = z - ( result.z_zernike );
    double diff_max     = z - ( result.z_max + surface_thickness) ; 
    zernike_surface_return return_value;
    return_value.z_at_0_0 = result.z_at_0_0;
    return_value.z_max    = result.z_max;
    return_value.z_zernike= result.z_zernike;
    if ((diff_zernike >= 0)  && (diff_max < 0) && (result.r_max <= surface.pupil_radius)) {
      return_value.inout = -1;
    }
    else{
      return_value.inout = 1;
    }
    return return_value;
  }


} // end zernike
