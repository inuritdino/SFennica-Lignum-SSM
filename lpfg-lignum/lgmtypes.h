#include <vector>

//Bit shift operator for V3f type
inline std::ostream& operator<< (std::ostream& out, V3f v){
  out << "(" << v.x << "," << v.y << "," << v.z << ")";
  return out;
}

//Declare the CylData class
class CylData {
public:
  CylData(): age(0),order(0),length(0.0),radius(0.0),start(0.0,0.0,0.0),\
	     end(0.0,0.0,0.0),axis(0.0,0.0,0.0),Wf0(0.0),Wf(0.0),\
	     hwR(0.0),parent(-1),extension(-1),children(){}
  void setAgeOrder(int A, int O){//Sets age and order info
    age = A;
    order = O;
    is_deleted = false;
    is_dead = false;
  }
  void setLength(float L){//Sets radius and all related entities
    length = L;
    lr = LR;
    radius = L * lr;//fixed ratio
    /* if(isnan(radius)){ */
    /*   std::cout << "ERROR: init with NaN radius " << std::endl; */
    /* } */
    af = AF;
    Wf0 = af * 2 * M_PI * length * radius;//Lignum formula
    Wf = Wf0;
    ksi = KSI;
    hwR = sqrt(ksi) * radius;
    //end = start + length*axis;//redefine the end point
  }
  void setAxis(V3f sta, V3f Ax){//Sets axis of the cyl
    start = sta;
    axis = Ax;
    end = start + length*axis;
  }
  void setParent(int p){//Sets the parent for the cyl
    parent = p;
  }
  void setChild(int new_child){//Adds a child of the cyl
    children.push_back(new_child);//
  }
  void setExtension(int ext){//Set the extension of the cyl
    extension = ext;
  }
  int age;// cyl age
  int order;// order of the branch cyl belongs to
  float length;// length of the cyl
  float radius;// radius of the cyl
  V3f start;// starting point
  V3f end;// ending point
  V3f axis;// axis, namely, <V3f>start - <V3f>end (a bit redundant)
  float Wf0;// initial foliage mass
  float Wf;// foliage mass of the cyl
  float hwR;// radius of the heartwood of the cyl
  int parent;// cyl's parent, i.e. a cyl it connects to directly
  int extension;//extension of the cyl
  std::vector<int> children;//children of the cyl
  bool is_deleted;//is this cyl deleted?
  bool is_dead;//is dead, i.e. included in the final structure
  float lr;//lr ratio for the cyl
  float af;//af const for the cyl
  float ksi;//ksi const for the cyl
};
std::ostream& operator<< (std::ostream& out, CylData cyl){
  out << cyl.age << " ";
  out << cyl.order << " ";
  out << cyl.length << " ";
  out << cyl.radius << " ";
  out << cyl.start << " ";
  out << cyl.end << " ";
  out << cyl.axis << " ";
  out << cyl.Wf0 << " ";
  out << cyl.Wf << " ";
  out << cyl.hwR << " ";
  out << cyl.parent << " ";
  out << cyl.extension << " ";
  out << "(";
  for(int i=0;i<cyl.children.size();i++){
    out << cyl.children[i];
    if(i < cyl.children.size()-1)
      out << ",";
  }
  out << ") ";
  out << cyl.is_deleted;
  out << std::endl;
  return out;
  /* if(cyl.is_deleted) */
  /*   out << "\tDELETED" << std::endl; */
  /* out << std::endl << "\tAge: " << cyl.age << std::endl; */
  /* out << "\tOrder: " << cyl.order << std::endl; */
  /* out << "\tLength: " << cyl.length << std::endl; */
  /* out << "\tRadius: " << cyl.radius << std::endl; */
  /* out << "\tStart: " << cyl.start << std::endl; */
  /* out << "\tEnd: " << cyl.end << std::endl; */
  /* out << "\tAxis: " << cyl.axis << std::endl; */
  /* out << "\tWf: " << cyl.Wf << std::endl; */
  /* out << "\tWf0: " << cyl.Wf0 << std::endl; */
  /* out << "\thwR: " << cyl.hwR << std::endl; */
  /* out << "\tparent: " << cyl.parent << std::endl; */
  /* //out << "\tis_child: " << cyl.is_child << std::endl; */
  /* out << "\textension: " << cyl.extension << std::endl; */
  /* out << "\tchildren: "; */
  /* for(int i=0; i<cyl.children.size(); i++){ */
  /*   out << cyl.children[i] << " "; */
  /* } */
  /* out << std::endl; */

  return out;
}

//Declare the Branch class
class Branch{
 public:
  std::vector<int> cyl_ind;
  int order;
  int parent_br;
};

std::ostream& operator<< (std::ostream& out, Branch br){
  out << "\torder: " << br.order << std::endl;
  out << "\tparent: " << br.parent_br << std::endl;
  out << "\tcyl's: ";
  for(int i = 0; i < br.cyl_ind.size(); i++){
    out << br.cyl_ind[i] << " ";
  }
  out << std::endl;
  return out;
}

//Declare the class for the Voxel space/grid
class VoxelBox {
public:
  VoxelBox(float xd = 0.03, float yd = 0.03, float zd = 0.03, \
	   V3f cent = V3f(0.0,0.0,0.0), int gsx = GSX, int gsz = GSZ)
  {
    vdx = xd;
    vdy = yd;
    vdz = zd;
    offset.x = cent.x - vdx/2.0 - ceil(gsx/2) * vdx;//NOTE: cent is by default (0,0,0) origin
    offset.y = cent.y - vdy/2.0;
    offset.z = cent.z - vdz/2.0 - ceil(gsz/2) * vdz;
  }
  //Voxel dimensions, in m
  float vdx;
  float vdy;
  float vdz;
  //VoxelBox offset, the corner point, the "leftmost" point of the grid
  V3f offset;
};

//Declare a class for 3D arrays of floats
class Array3D {
 public:
  size_t m_width, m_height;
  std::vector<float> m_data;
public:
  Array3D(size_t x, size_t y, size_t z, float init = 0.0):
    m_width(x), m_height(y), m_data(x*y*z, init)
  {}
  float& operator() (size_t x, size_t y, size_t z) {
    return m_data.at(x + y*m_width + z*m_width*m_height);
  }
  void Reset(float reinit = 0.0){
    m_data.assign(m_data.size(),reinit);
  }
};

//Class Line
class Line {
public:
  V3f point;//a point
  V3f axis;//a unit vector pointing the direction
  //Any point on the line is defined as: P = start + C*axis, where C is a coefficient
  Line(V3f p, V3f a): point(p), axis(a) {}
};


//Class Firmament

class Firmament{
public:
  Firmament(int no_incl=9, int no_azim=9,float rad_plane=1200.0);
  /* void resize(int no_incl, int no_azim, double diffuse_rad_plane); */
  /* void setDiffuseRadiation(const double rad); */
  /* void setDirectRadiation(const double rad) { directRadPlane = rad; } */
  /* void setSunPosition(const vector<double>& v); */
  /* vector<double>  getSunPosition() { return sunPosition; } */
  /* MJ directRadiation(vector<double>& direction); */
  /* MJ diffuseRadiationSum( const vector<double>& v); */
  /* int numberOfRegions() const { return numOfSectors; } */
  /* MJ diffuseRegionRadiationSum(int n, vector<double>& direction)const; */
  /* MJ diffuseHalfRegionRadiationSum(int n, vector<double>& direction)const; */
  /* MJ directHalfRegionRadiationSum(vector<double>& direction); */
  /* MJ diffusePlaneSensor(void) { return diffuseRadPlane; } */
  /* MJ diffuseBallSensor(void) { return diffuseRadBall; } */
  /* MJ diffuseForestRegionRadiationSum(int n, float z, float x,  float la, float ke, */
  /* 				     float H, float Hc, */
  /* 				     vector<double>& direction,double density); */
/*   void outDiff() { */
/*     cout << diffuseRad << endl; */
/*   } */
/*   void outAz() { */
/*     cout << zoneAzims << endl; */
/*   } */
/*   double getInclination(int n); */
/*   double getAzimuth(int n); */
/*   int getInclinationIndex(int n) { */
/*     if(n < 0 || n > (numOfSectors - 2) )  // numOfSectors - 1 == zenith */
/*       return -1; */
/*     else */
/*       return inclinationIndex[n]; */
/*   } */

/*   int getAzimuthIndex(int n) { */
/*     if(n < 0 || n > (numOfSectors - 2) ) // numOfSectors - 1 == zenith */
/*       return -1; */
/*     else */
/*       return azimuthIndex[n]; */
/*   } */
/*   double getSectorArea(int n) { */
/*     if(n < 0 || n > num_of_incl - 1) */
/*       return -1.0; */
/*     else */
/*       return areasByInclination[n]; */
/*   } */
 
/*   int getAzimDivision(int n) */
/* { */
/*     if(n < 0 || n > num_of_incl - 1) */
/*       return -1; */
/*     else */
/*       return azimDivisions[n]; */
/* } */
  

/*   int getNoOfAzimuths() { return num_of_azim; } */
/*   int getNoOfInclinations() { return num_of_incl; } */
/*   PositionVector getDirection(int n)const; */

/*   void outInclinations() { */
/*     int line = 1; */
/*     for(int i = 0; i < num_of_incl; i++) { */
/*       cout << inclinations[i] << " "; */
/*       line++; */
/*       if(line == 10) { */
/* 	cout << endl; */
/* 	line = 1; */
/*       } */
/*     } */
/*     if(line != 1) */
/*       cout << endl; */
/*   } */

/*   vector<pair<double,double> > getIncAz() { */
/*     vector<pair<double,double> > result(numOfSectors); */
/*     for(int i = 0; i < numOfSectors; i++) { */
/*       pair<double,double> p(getInclination(i),getAzimuth(i)); */
/*       result[i] = p; */
/*     } */
/*     return result; */
/*   } */

//protected:
  int num_of_incl;
  int num_of_azim;
  float diffuseRadScale;
  float directRadPlane, diffuseRadPlane, diffuseRadBall;
  float diffuseRadZenith;
  std::vector<float> sunPosition;

  std::vector<std::vector<float> > zoneAzims;
  std::vector<float> inclinations;
  std::vector<int> azimDivisions;
  std::vector<std::vector<float> > diffuseRad; //stores the radiation of sectors
  std::vector<int> inclinationIndex;   //tabulation of index of inclination in
                                   //diffuseRad as a function of running segment no.
  std::vector<int> azimuthIndex;       //as inclinationIndex but for azimuths;
  std::vector<float> dir_x;           //x-component of std::vector pointing to center
                                  //of ith sector (indexed by number of sector)
  std::vector<float> dir_y;          //as dir_x but for y
  std::vector<float> dir_z;          //as dir_x but for z
  std::vector<float> areasByInclination; //areas sectors in for each inclination
  float thetZ;
  int numOfSectors;
  float deltaIncl, halfDeltaIncl;
  float standDensity;
};
