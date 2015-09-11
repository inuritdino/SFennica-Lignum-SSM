#include <lpfgall.h>
#include <math.h>
#include <ctime>
#include <climits>
#ifndef FSTREAM_H
#include <fstream>
#define FSTREAM_H 1
#endif
//Random distr's functions
#include "randist.cpp"
//LIGNUM-specific constants
#include "user.h"//User specified input
//Parameters
#if PERTTUNEN1998 != 0//Include parameters from Perttunen 1998.
#include "Perttunen1998.par"
#elif SIEVANEN2008 != 0//Include parameters from Sievanen 2008.
#include "Sievanen2008.par"
#else
#include "lgmconst.h"//Default parameters
#endif //PERTUNEN1998 != 0 AND SIEVANEN2008 != 0
//LIGNUM types and classes defined
#include "lgmtypes.h"
//Main algorithms and routines
#include "methods.cpp"
//Shadow propagation model functions and constants
#include "shadow_propagation.cpp"
//Firmament radiation
//#include "firmament_radiation.cpp"
//MTG file functions
#include "mtg.cpp"
//Scatter/distribution generation
#include "scatter.cpp"
//GSL include
//#include "test.h"
//max number of cyl's
#define MAXNCYL 100000

// **************************************
// INITIALIZATION of LIGNUM specific data

// Root mass
float Wr = 0.0;// Aggregated root mass for the whole tree

//Number of iterations and timer
int nIter,Year;
std::clock_t startT,stopT;
std::vector<int> futileYears;

//Cyls information
std::vector<CylData> cyls;
int cyl_resize_count = 0;
int nCyl = 1;
int nShoot = 0;//number of shoots of this year
int nRemovedSpace = 0;
int nRemovedShed = 0;
int nRemovedMass = 0;
int nRemovedRupt = 0;
int nRemovedEnv = 0;
int nRemovedEnvTot = 0;
float maxRupt = 0.0;

//First cyl's starting position, length and radius, default orientation
V3f sta(0.0,0.0,0.0);
V3f Dir(0.0,1.0,0.0);

//Voxel space
std::vector<int> gSize (3);
VoxelBox vb(VSX,VSY,VSZ,sta);
Array3D vGrid(GSX,GSZ,GSY,0.0);//Grid of voxels(type is float)

//Light conditions
std::vector<float> fL;

//Respiration and photosynthetic productions
float P,M;

//Main balance equation solution entities
//Polynomial coefficients of the equation, 3rd power is max
float A[4] = {0.0};//A[0] - free term, A[1]*Lambda^1 etc.

//Main variable of the equation
float Lambda;

//FATAL ERRORS tracker
int fatal_error = 0;
//LOG file stream
std::ofstream logfile;

//Module declarations: Cyl and Bud
module Bud;//Cyl's bud complex, gets substituted with T- and L-buds.
module Cyl(int);//internode=cylinder

//Declare functions
float photosynthesis(float [], int);
int remove_buds(bool&, bool&, int);
int fl_ip(float [],std::vector<float> &, const int);
float tot_sap_mass(int);
float tot_fol_mass(int);
float pf_emp(int);
int age_increment(const int);
int num_L_buds(int);
bool is_empty_array(std::vector<bool>&);
float sap_area(int);
float sap_area_fol(int);
void check_num_cyls(int);
bool traversed_in_space(int);
bool branch_grow_lr(int);
int branch_shed_by_mass();
int branch_shed_by_rupture();
float rupture_modulus(int,std::vector<int>&);
V3f mass_center(std::vector<int> &);
float mass_of_cyl(int);
float mass_of_cyls(std::vector<int> &);
int get_all_offsets(int,std::vector<int> &);
float cyls_ground_proj_area(std::vector<int> &);
int num_deleted(int,std::vector<CylData>&);
int overall_info(int,int,std::vector<CylData>&,V3f,float,std::vector<Branch>&);
int pipe_model(std::vector<CylData>&,int,int,float[],\
	       std::vector<float>&,std::vector<float>&);
int pipe_model_debug(std::vector<CylData>&,int,std::vector<bool>&);
int pipe_model_debug_print(int,bool);
bool is_same_age(std::vector<bool>&);
int memb_in_layer(std::vector<bool>&);
int postpone_layer_members(int,std::vector<bool>&,std::vector<int>&);
int print_deleted_cyls(int);
