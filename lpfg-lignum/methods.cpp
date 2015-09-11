#define TOL (1e-05)
//Functions' declarations
float sap_cross_area(std::vector<CylData> &, int);
int propagate_branch_backward(std::vector<CylData> &,int ,std::vector<int>&, \
			      std::vector<int> & , bool &);
int branch_shedding(std::vector<CylData> & , int );
int shift_branch_offsets(std::vector<CylData> &,int);
int branch_bending(std::vector<CylData> &,std::vector<int> &);
int branch_bending_wrapper(std::vector<CylData> & );
int move_cyl_prop(std::vector<CylData> &, int, int);
int fill_cyl_gaps(std::vector<CylData> &,int);
int delete_cyls(std::vector<CylData> &, std::vector<int> &, int);
float tree_height(std::vector<CylData> & , int , float );
float breast_height_diam(std::vector<CylData> & , int , V3f );
float min_dist_2lines(V3f&,V3f&,float&, float&, Line&, Line&);
float min_dist_line_seg(V3f& , V3f& ,float& ,float& ,Line& , V3f& , V3f& );
V3f lpfg_to_normal_orientation(V3f& );
int modify_coord_xy(V3f, V3f &, V3f &);
int modify_coord_xz(V3f , V3f & , V3f & , const V3f );
float projection_angle(V3f , V3f , V3f );
V3f rotate_v3f(V3f &,V3f& ,float);
double my_lambda_function (double , float []);
double r8_abs ( double );
double r8_epsilon ( );
double my_zero ( double , double, double , float [] );

//Some inline handy functions
inline float rad_to_deg(float rad){return ((180.0/M_PI)*rad);}
inline float deg_to_rad(float deg){return ((M_PI/180.0)*deg);}
//Set a random direction branching from the vector v, for rotations purposes
inline V3f rnd_dir(V3f v){
  V3f out(v.x+ran(2.0)-1,v.y+ran(2.0)-1,v.z+ran(2.0)-1);
  return (out.Normalize());
}

// ****************** SAPWOOD AREA **************************
float sap_cross_area(std::vector<CylData> & cyls, int c)
{
  return (M_PI*(cyls[c].radius*cyls[c].radius-cyls[c].hwR*cyls[c].hwR));
}

// ****************** BRANCH SHEDDING (BY FOLIAGE) ****************
int propagate_branch_backward(std::vector<CylData> & cyls,int tip_cyl,	\
			      std::vector<int>& child_root,		\
			      std::vector<int>& br_cyls, bool& no_foliage)
{//Propagate the branch from its tip TIP_CYL to the base=root.
  //Determine children branches root cyls in CHILD_ROOT.
  //Determine the cyl's of the branch itself in BR_CYLS.
  //Stop propagating, whenever(!) a cyl in a branch has some foliage left.
  //Return NO_FOLIAGE bool.
  if((cyls[tip_cyl].extension != -1) || (cyls[tip_cyl].children.size() != 0)){
    //not a tip-cyl return silently
    return 1;
  }
  if(cyls[tip_cyl].is_deleted)
    return 2;

  //Make zero size for the output vectors
  child_root.resize(0);
  br_cyls.resize(0);
  //Check the current root cyl
  int cc = tip_cyl;
  br_cyls.push_back(cc);
  if(cyls[cc].Wf < FOLTHR)
    no_foliage = true;
  else
    no_foliage = false;
  while( no_foliage && (cyls[cc].parent > -1) &&\
	 (cyls[cyls[cc].parent].extension == cc))
    {//No foliage and CC is an extension of its parent
      cc = cyls[cc].parent;//Shift to the parent
      br_cyls.push_back(cc);
      //Identify the children and form the output array
      for(int i=0; i<cyls[cc].children.size(); i++){
	child_root.push_back(cyls[cc].children[i]);
      }
      if(cyls[cc].Wf < FOLTHR)
	no_foliage = true;
      else
	no_foliage = false;
    }

  return 0;
}

int branch_shedding(std::vector<CylData> & cyls, int nCyl, int & nRemovedShed)
{//Shed the branches without foliage
  //Find the tip cyl's(highest order or recently added ones, most probably shoots)
  std::clock_t startT = clock();
  bool deleted = true;
  bool no_foliage;
  std::vector<int> child_root;
  std::vector<int> cyls_to_shed;
  std::vector<bool> protected_cyls (nCyl,false);
  unsigned int k;//number of tip 1st order cyl's to shed
  float shed_pois_mu;
  while(deleted)
    {//Stop iterations when there were no deleted cyls/branches (deleted == false)
      deleted = false;
      for(int i=0; i<nCyl; i++){
	//do not consider deleted cyl's, protected in current year or dead
	if(cyls[i].is_deleted || protected_cyls[i] || cyls[i].is_dead)
	  continue;
	if(!propagate_branch_backward(cyls,i,child_root,cyls_to_shed,no_foliage)){
	  //If there are no children and no foliage in the branch cyls_to_shed, delete it
	  if((child_root.size() == 0) && no_foliage){
	    if(cyls[i].age >= SHEDYEAR){
	      //Protected in prev years may come here
	      if(cyls[i].order == 1){
		//First-order branches special dealing
		//1. Generate random number of cyls to shed starting from the tip
		//Uniform distribution OR Poisson with parameterized MEAN.
		// k is the num of cyl's to shed, starting from tip
		if(SHEDDIST == 1){//Uniform distribution
		  k = (unsigned int)ceil( ran((float) cyls_to_shed.size() ) );
		}
		else if(SHEDDIST == 2){//Poisson with MEAN = SHEDMU
		  if((float)SHEDMU > 10.0){
		    std::cout << "Warning: SHEDMU cannot be larger than 10. Setting to 10." \
			      << std::endl;
		    shed_pois_mu = 10.0;
		  }
		  else
		    shed_pois_mu = (float)SHEDMU;
		  do{
		    k = ran_poisson(shed_pois_mu);
		    //std::cout << "k_poiss = " << k << std::endl;
		  }
		  while ( k > cyls_to_shed.size() );
		}
		else{//Otherwise, shed everything
		  k = (unsigned int)cyls_to_shed.size();
		}
		// std::cout << "Shed " << k << "/" << cyls_to_shed.size() << std::endl;
		// for(int j = 0; j < cyls_to_shed.size(); j++){
		//   std::cout << cyls_to_shed[j] << " ";
		// }
		// std::cout << std::endl;

		//2. Mark the non-shed cyls as dead and reshape the cyls_to_shed
		//k = 1;//debugging condition
		//Mark the non-shed cyls as dead, thus protecting them from shedding ever
		for(int j = k; j < cyls_to_shed.size(); j++)
		  cyls[cyls_to_shed[j]].is_dead = true;
		//Delete the non-shed ones from cyls_to_shed
		cyls_to_shed.erase(cyls_to_shed.begin()+k,cyls_to_shed.end());
		
		//cyls_to_shed.size() MUST be equal to k
		// std::cout << "#cyls to shed " << cyls_to_shed.size() << std::endl;
		// for(int j = 0; j < cyls_to_shed.size(); j++){
		//   std::cout << cyls_to_shed[j] << " ";
		// }
		// std::cout << std::endl;

		if(cyls_to_shed.size() != k){
		  std::cout << "FATAL ERROR: asked " << k << " to shed, but " \
			    << cyls_to_shed.size() << " will!" << std::endl;
		}
	      }
	      //Shed the cyls
	      delete_cyls(cyls,cyls_to_shed,nCyl);
	      nRemovedShed += cyls_to_shed.size();
	    }
	    else{
	      //Make protected if not supposed to be shed
	      for(int j = 0; j < cyls_to_shed.size(); j++){
		protected_cyls[cyls_to_shed[j]] = true;
	      }
	    }
	    deleted = true;
	    break;
	  }
	  else{
	    continue;
	  }
	}
      }
    }
  // std::cout << "Shedding time: " << \
  //   ((double)clock()-(double)startT)/(double)CLOCKS_PER_SEC << " sec." << std::endl;
  return 0;
}

//Some technical functions
bool do_belong_to_array(std::vector<int>&array, int some)
{
  for(int i = 0; i < array.size(); i++){
    if(some == array[i])
      return true;
  }
  return false;
}

// ***************** BRANCH BENDING ******************************

int shift_branch_offsets(std::vector<CylData> &cyls,int cyl)
{//Shifts the offset cyl's of the branch due to the change of orientation of CYL also
 //belonging to the branch
  std::vector<int> parent_cyls(1,cyl);
  std::vector<int> offset_cyls;
  std::vector<int> to_remove (0);
  int cc,j = 0,i;
  while( j < parent_cyls.size() ){
    cc = parent_cyls[j];
    //Get all offsets
    if(cyls[cc].extension > -1)
      offset_cyls.push_back(cyls[cc].extension);
    for(i=0; i<cyls[cc].children.size(); i++)
      offset_cyls.push_back(cyls[cc].children[i]);
    //Redefine the start and end points of the offset
    for(i=0; i<offset_cyls.size(); i++){
      cyls[offset_cyls[i]].setAxis(cyls[cc].end,cyls[offset_cyls[i]].axis);
      parent_cyls.push_back(offset_cyls[i]);
    }
    offset_cyls.resize(0);
    j++;
  }

  return 0;
}

int branch_bending(std::vector<CylData> & cyls,std::vector<int> &br_cyls)
{//Provides for a heuristic procedure for the branch bending.

  if(br_cyls.size() == 1)//no need to proceed
    return 2;
  //The first cyl in br_cyls is the base cyl of the branch
  //Define the rotation axis by cross prod between the parent and 1st cyl of the
  //branch.

  V3f rot_axis,cyl_axis,old_axis;
  //Define rotation axis as the rotation in the plane formed by the 1st cyl of the
  //branch and vertical axis=(0,1,0).
  rot_axis = V3f(0.0,1.0,0.0) % cyls[br_cyls[0]].axis;
  // rot_axis =								\
  //   cyls[cyls[br_cyls[0]].parent].axis % cyls[br_cyls[0]].axis;
  rot_axis.Normalize();
  if(rot_axis.Length() < TOL)//rot_axis is vertical => no bending
    return 1;
  
  float ang_incr;
  float first_bra_ang;
  
  // Start from the root of the branch (its first cyl)
  for(int i=0; i<br_cyls.size()-1; i++){
    first_bra_ang =							\
      acos(cyls[cyls[br_cyls[0]].parent].axis * cyls[br_cyls[i]].axis);
    first_bra_ang = rad_to_deg ( first_bra_ang );
    if(ZETASD < TOL)
      ang_incr = (float)ANGINCR;
    else
      ang_incr = ran_gauss_any((float)ANGINCR,(float)ZETASD);
    if((ang_incr + first_bra_ang) > MAXANG){
      ang_incr = MAXANG - first_bra_ang;
    }
    cyl_axis = rotate_v3f(cyls[br_cyls[i]].axis,rot_axis,deg_to_rad(ang_incr));
    cyls[br_cyls[i]].setAxis(cyls[br_cyls[i]].start,cyl_axis);
    shift_branch_offsets(cyls,br_cyls[i]);
  }

  return 0;
}

int branch_bending_wrapper(std::vector<CylData> & cyls)
{//Perform the branch bending in the end of each year by observing the overall
 //structure of the tree.
  // Go along the trunk and initialize the 1st order branches
  int cc = 0;
  std::vector<int> br;//branch cyl's
  std::vector<int> bcyls;//base cyl's - initiators of branches

  //Extract 1-order bases
  while(cyls[cc].extension > -1){
    for(int i=0;i<cyls[cc].children.size();i++)
      bcyls.push_back(cyls[cc].children[i]);
    cc = cyls[cc].extension;
  }
  for(int i=0;i<cyls[cc].children.size();i++)
    bcyls.push_back(cyls[cc].children[i]);

  //Extract branches and bend them
  while(bcyls.size()){
    cc = bcyls[0];//Rename the first cyl in BCYLS
    br.push_back(cc);

    while(cyls[cc].extension > -1){
      for(int i=0;i<cyls[cc].children.size();i++)
	bcyls.push_back(cyls[cc].children[i]);
      cc = cyls[cc].extension;
      br.push_back(cc);
    }
    for(int i=0;i<cyls[cc].children.size();i++)
      bcyls.push_back(cyls[cc].children[i]);
    
    branch_bending(cyls,br);
    br.resize(0);
    bcyls.erase(bcyls.begin());//Remove the processed base
  }

  return 0;
}


// ***************** DELETE CYLS **************************
int move_cyl_prop(std::vector<CylData> & cyls, int from, int to)
{//Move the properties of cyl FROM to cyl TO in the CYLS
  //We need to change the parent and children information of the cyl to move,
  //i.e. FROM.

  // 1. FROM has parent, change its extension or children information to TO
  if(cyls[from].parent > -1){
    if(cyls[cyls[from].parent].extension == from){//extension of the parent
      cyls[cyls[from].parent].extension = to;
    }
    else{//among the children of the parent
      for(int i = 0; i < cyls[cyls[from].parent].children.size(); i++){
	//Go through all children until find the FROM and replace it. Break then.
	if(cyls[cyls[from].parent].children[i] == from){
	  cyls[cyls[from].parent].children[i] = to;
	  break;
	}
      }
    }
  }
  // 2. FROM has extension, change its parent information to TO
  if(cyls[from].extension > -1){
    cyls[cyls[from].extension].parent = to;
  }
  // 3. FROM has children, change their parent information to TO
  for(int i = 0; i < cyls[from].children.size(); i++){
    cyls[cyls[from].children[i]].parent = to;
  }
  // 4. Move FROM to TO
  cyls[to] = cyls[from];
  // 5. Remove FROM
  cyls[from].is_deleted = true;

  return 0;
}

int fill_cyl_gaps(std::vector<CylData> & cyls,int nCyl)
{//Fill in the gaps of deleted cyls with not deleted ones by adjusting cyls.
  //Gaps are identified with is_deleted flag of the CylData.

  int j;//will indicate the closest to gap meaningful position
  int i = 0;
  while(i < nCyl){
    if(cyls[i].is_deleted){//gap
      j = i + 1;//next to the gap, possible meaningful position
    }
    while(cyls[j].is_deleted)//iterate forward through the gaps
      j++;
    //Now j contains the first after the i-th gap meaningful position
    //Move j'th cyl to i'th position
    move_cyl_prop(cyls,j,i);
    i++;
  }
  
  return 0;
}

int delete_cyls(std::vector<CylData> & cyls, std::vector<int> & cyls_to_remove, int nCyl)
{//Remove CYLS from the tree. Removing is done via a special flag in
 //CylData. Additionally, the parent's information of the removed cyl is cleared so
 //that the removed cyl is no longer among extension/children of the parent
  if(cyls_to_remove.size() == 0)
    return 1;
  int i=0;
  while( 1 ){
    //Delete the cyl
    // if(cyls[cyls_to_remove[i]].is_deleted){
    //   std::cerr << "Warning: cyl is already removed: " << cyls_to_remove[i] << std::endl;
    // }
    cyls[cyls_to_remove[i]].is_deleted = true;
    //Adjust topological information of the parent
    if( cyls[cyls[cyls_to_remove[i]].parent].extension == cyls_to_remove[i] )
      cyls[cyls[cyls_to_remove[i]].parent].extension = -1;//among extensions
    else{//among children
      for(int j=0; j<cyls[cyls[cyls_to_remove[i]].parent].children.size(); j++){
	if(cyls[cyls[cyls_to_remove[i]].parent].children[j] == cyls_to_remove[i]){
	  cyls[cyls[cyls_to_remove[i]].parent].children.\
	    erase(cyls[cyls[cyls_to_remove[i]].parent].children.begin()+j);
	  break;
	}
      }
    }
    i++;
    if(i == cyls_to_remove.size())
      break;
  }

  //fill_cyl_gaps(cyls,nCyl);
  
  return 0;
}

// ********************* TREE HEIGHT AND DIAMETER *********************
float tree_height(std::vector<CylData> & cyls, int nCyl, float lowest)
{//Calculate the tree height
  float H = 0.0;
  for(int i=0;i<nCyl;i++){
    if(cyls[i].is_deleted)
      continue;
    if(H < cyls[i].end.y)
      H = cyls[i].end.y;
  }

  return (H-lowest);
}

float breast_height_diam(std::vector<CylData> & cyls, int nCyl, V3f sta)
{//Calculate the breast height (1.3 m) diameter of the tree.
  float D = 0.0;
  for(int i=0; i<nCyl; i++){
    if(cyls[i].is_deleted || cyls[i].order != 0)
      continue;
    if( (fabs(cyls[i].start.y-sta.y) < 1.3) && fabs(cyls[i].end.y-sta.y) >= 1.3 ){
      D = 2 * cyls[i].radius;
      break;
    }
  }
  
  return D;
}


// *************** DISTANCES BETWEEN LINES *************************
float min_dist_2lines(V3f& x1,V3f& x2,float& c1, float& c2, Line& l1, Line& l2)
{//Calculate the minimum distance between two lines (returned).
  // X1 - is a point on Line L1 through which the minimal distance line passes
  // X2 - the same for Line L2
  // C1, C2 are the line coefficients
  // Any point on a line can be represented as: P = Line.start + C*Line.axis .
  // This code was adopted from the geomalgorithms.com. See at:
  // http://geomalgorithms.com/a07-_distance.html
  V3f w0 = l2.point - l1.point;
  float a = l1.axis * l1.axis;
  float b = l1.axis * l2.axis;
  float c = l2.axis * l2.axis;
  float d = l1.axis * w0;
  float e = l2.axis * w0;

  if( fabs(a * c - b * b) < 1e-05 ){//collinear=parallel lines
    c1 = 0.0;//fix point on the 1st line
    c2 = - d / b;
  }
  else{
    c1 = (-b * e + c * d) / (a * c - b * b);
    c2 = (-a * e + b * d) / (a * c - b * b);
  }
  x1 = l1.point + c1*l1.axis;
  x2 = l2.point + c2*l2.axis;
  
  V3f out = x2 - x1;
  return (out.Length());
}

float min_dist_line_seg(V3f& x1, V3f& x2,float& c1,float& c2, \
			Line& ln, V3f& q1, V3f& q2)
{//Calculate minimum distance between line and segment (returned);
  //Notations are similare to that of MIN_DIST_2LINES() function.
  //Except:
  //Q1 and Q2 denote the ends of the segment
  //When the shortest path is beyond the segment's span, the perpendicular to the
  //segment is calculated (not to the line) for recalculation of the path.
  //This path will be crossing either ends of the segment
  
  Line ls(q1,q2-q1);//the segment's line
  min_dist_2lines(x1,x2,c1,c2,ln,ls);
  // if C2 < 0 or > 1 than the shortest path beyond the segment's span
  if(c2 < 0.0){
    //The beginning of the segment
    c2 = 0.0;
    x2 = q1;
    c1 = - ((ln.point-q1)*ls.axis) / (ln.axis*ls.axis);
    x1 = ln.point + c1*ln.axis;
  }
  if(c2 > 1.0){
    //The end of the segment
    c2 = 1.0;
    x2 = q2;
    c1 = - ((ln.point-q2)*ls.axis) / (ln.axis*ls.axis);
    x1 = ln.point + c1*ln.axis;
  }
  V3f out = x2 - x1;
  return (out.Length());
}


// ******************* 3D ORIENTATION ****************************
V3f lpfg_to_normal_orientation(V3f& vector)
{
  V3f out = V3f(vector.x,-vector.z,vector.y);
  return out;
}

int modify_coord_xy(V3f vector, V3f & newX, V3f & newY)
{//Modifies the original X,Y,Z coord system such that
  // 1. newX = vector's projection onto XY-plane
  // 2. Y-axis = newY is rotated to correspond to the previous change
  // 3. Z is left unchanged.
  // There is a special procedure served when vector coincides with +/-Z direction.
  // The function account on the special LPFG orientation, namely,
  // X = Xlpfg, Y = -Zlpfg, Z = Ylpfg.
  // NOTE: series of modifications are applied to vector, but it is NOT returned.

  //For reference
  //V3f X = V3f(1.0,0.0,0.0);
  V3f Y = V3f(0.0,0.0,-1.0);
  V3f Z = V3f(0.0,1.0,0.0);
  //projection of vector onto XY
  float tmp = vector.y;
  vector.y = 0.0;
  vector.Normalize();
  if(vector.Length() < TOL){//Vertical vector = Z/-Z
    std::cout << "WARNING: VERTICAL VECTOR!" << std::endl;
    newX = vector;
    Z = newX % Y;
    newY = Y;
  }
  else{//Normal procedure
    //new X-axis
    newX = vector;
    //new Y-axis
    newY = Z % newX;//Y = cross product of Z and X
  }
  // ERRORS CHECK
  //Orthogonal check
  if((fabs(newX*newY) > TOL) || (fabs(newX*Z) > TOL)\
     || (fabs(Z*newY) > TOL)){
    std::cout << "Error(modify_coord_xy): Basis is not orthogonal: " \
	      << newX << ", " << newY << ", " << Z << std::endl;
    return 1;
  }
  //Length check
  if((newX.Length() > (1.0 + TOL)) || (newY.Length() > (1.0 + TOL))){
    std::cout << "Error(modify_coord_xy): Basis length is not one." << std::endl;
    return 2;
  }
  return 0;
}

int modify_coord_xz(V3f vector, V3f & newX, V3f & newZ, const V3f newY)
{//Modifies the coordinate system obtained by MODIFY_COORD_XY() function such that
  // 1. newX = vector
  // 2. Z-axis is rotated to correspond to the previous change
  // 3. Y(newY) is left unchanged (see modify_coord_xy)
  // The function account on the special LPFG orientation, namely,
  // X = Xlpfg, Y = -Zlpfg, Z = Ylpfg.
  // NOTE: series of modifications are applied to vector, but it is NOT returned.

  //For reference
  // V3f X = V3f(1.0,0.0,0.0);
  // V3f Y = V3f(0.0,0.0,-1.0);
  // V3f Z = V3f(0.0,1.0,0.0);
  //new X-axis
  newX = vector;
  //new Z-axis
  newZ = newX % newY;
  // ERRORS CHECK
  //Orthogonal check
  if((fabs(newX*newY) > TOL) || (fabs(newX*newZ) > TOL)\
     || (fabs(newZ*newY) > TOL)){
    std::cout << "Error(modify_coord_xz): Basis is not orthogonal: " \
	      << newX << ", " << newY << ", " << newZ << std::endl;
    return 1;
  }
  //Length check
  if((newX.Length() > (1.0 + TOL)) || (newZ.Length() > (1.0 + TOL))){
    std::cout << "Error(modify_coord_xz): Basis length is not one." << std::endl;
  }
  return 0;
}

float projection_angle(V3f vector, V3f ortho1, V3f ortho2)
{//Calculates the angle of the projection of the vector onto the plane determined by
 //ortho1 and ortho2 with the ortho1 direction.

  //Plane's normal
  V3f normal = ortho1 % ortho2;
  //Find the projection by 2 consecutive cross products
  vector = vector % normal;
  vector = normal % vector;
  vector.Normalize();
  //Find angle
  return (atan2f(vector * ortho2, vector * ortho1));
}

V3f rotate_v3f(V3f &vector,V3f& rot_axis,float angle)
{
  // Rotate vector around rot_axis for the specified angle(in rad). This is done
  // only when vector operations on <V3f> type are defined (they are defined by
  // default).

  //Cosine and sine of the angle
  float C = cos((double)angle);
  float S = sin((double)angle);
  //Define the rotation matrix by rows (3 in total)
  rot_axis.Normalize();//normalize the rotation matrix
  float X = rot_axis.x;
  float Y = rot_axis.y;
  float Z = rot_axis.z;
  V3f row1(X*X+(1-X*X)*C,			\
	   X*Y*(1-C)-Z*S,\
	   X*Z*(1-C)+Y*S);
  V3f row2(X*Y*(1-C)+Z*S,\
	   Y*Y+(1-Y*Y)*C,\
	   Y*Z*(1-C)-X*S);
  V3f row3(X*Z*(1-C)-Y*S,\
	   Y*Z*(1-C)+X*S,\
	   Z*Z+(1-Z*Z)*C);
  //vector rotated is a result of matrix multiplication of Rotation matrix and the
  //vector that needs to be rotated. Thus, we can define the resulting vector with
  //components obtained by dot products of the corresponding matrix and old vector
  //components.
  V3f out;
  out.x = row1 * vector;
  out.y = row2 * vector;
  out.z = row3 * vector;
  
  return out;
}

// ******************* LAMBDA-EQ SOLUTION ****************************
double my_lambda_function (double x, float A[])
{
  return (A[0]+A[1]*x+A[2]*x*x+A[3]*x*x*x);
}

// ==== THE CODE BELOW IS FOR BRENT's METHOD ====
// OBTAINED FROM: http://people.sc.fsu.edu/~jburkardt/cpp_src/brent/brent.html
// ALTHOUGH A BIT ADOPTED

//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  }
  else
  {
    value = - x;
  }
  return value;
}

//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}


//****************************************************************************80

double my_zero ( double a, double b, double t, float Acoef[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZERO seeks the root of a function F(X) (my_lambda_function) in an
//    interval [A,B].
//
//  Discussion:
//
//    The interval [A,B] must be a change of sign interval for F.
//    That is, F(A) and F(B) must be of opposite signs.  Then
//    assuming that F is continuous implies the existence of at least
//    one value C between A and B for which F(C) = 0.
//
//    The location of the zero is determined to within an accuracy
//    of 6 * MACHEPS * r8_abs ( C ) + 2 * T.
//
//    Thanks to Thomas Secretin for pointing out a transcription error in the
//    setting of the value of P, 11 February 2013.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 February 2013
//
//  Author:
//
//    Original FORTRAN77 version by Richard Brent.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Richard Brent,
//    Algorithms for Minimization Without Derivatives,
//    Dover, 2002,
//    ISBN: 0-486-41998-3,
//    LC: QA402.5.B74.
//
//  Parameters:
//
//    Input, double A, B, the endpoints of the change of sign interval.
//
//    Input, double T, a positive error tolerance.
//
//    Input, func_base& F, the name of a user-supplied c++ functor
//    whose zero is being sought.  The input and output
//    of F() are of type double.
//
//    Output, double ZERO, the estimated value of a zero of
//    the function F.
//
{
  double c;
  double d;
  double e;
  double fa;
  double fb;
  double fc;
  double m;
  double macheps;
  double p;
  double q;
  double r;
  double s;
  double sa;
  double sb;
  double tol;
//
//  Make local copies of A and B.
//
  sa = a;
  sb = b;
  fa = my_lambda_function ( sa, Acoef );
  fb = my_lambda_function ( sb, Acoef );

  c = sa;
  fc = fa;
  e = sb - sa;
  d = e;

  macheps = r8_epsilon ( );

  for ( ; ; )
  {
    if ( r8_abs ( fc ) < r8_abs ( fb ) )
    {
      sa = sb;
      sb = c;
      c = sa;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    tol = 2.0 * macheps * r8_abs ( sb ) + t;
    m = 0.5 * ( c - sb );

    if ( r8_abs ( m ) <= tol || fb == 0.0 )
    {
      break;
    }

    if ( r8_abs ( e ) < tol || r8_abs ( fa ) <= r8_abs ( fb ) )
    {
      e = m;
      d = e;
    }
    else
    {
      s = fb / fa;

      if ( sa == c )
      {
        p = 2.0 * m * s;
        q = 1.0 - s;
      }
      else
      {
        q = fa / fc;
        r = fb / fc;
        p = s * ( 2.0 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0 ) );
        q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 );
      }

      if ( 0.0 < p )
      {
        q = - q;
      }
      else
      {
        p = - p;
      }

      s = e;
      e = d;

      if ( 2.0 * p < 3.0 * m * q - r8_abs ( tol * q ) &&
        p < r8_abs ( 0.5 * s * q ) )
      {
        d = p / q;
      }
      else
      {
        e = m;
        d = e;
      }
    }
    sa = sb;
    fa = fb;

    if ( tol < r8_abs ( d ) )
    {
      sb = sb + d;
    }
    else if ( 0.0 < m )
    {
      sb = sb + tol;
    }
    else
    {
      sb = sb - tol;
    }

    fb = my_lambda_function ( sb, Acoef );

    if ( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) )
    {
      c = sa;
      fc = fa;
      e = sb - sa;
      d = e;
    }
  }
  return sb;
}

