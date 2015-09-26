#ifndef FSTREAM_H
#include <fstream>
#define FSTREAM_H 1
#endif


int find_branch_by_order(std::vector<Branch> & br, int order, \
			 std::vector<int> & br_ind)
{
  br_ind.resize(0);
  for(int i = 0; i < br.size(); i++){
    if(br[i].order == order)
      br_ind.push_back(i);
  }
  // for(int i = 0; i < br_ind.size(); i++)
  //   std::cout << br_ind[i] << " ";
  // std::cout << std::endl;
  
  return 0;
}

int parent_branch(std::vector<CylData> & cyls, std::vector<Branch> & br, int cc)
{//Find the parent branch of a branch cotaining cyl CC as the base cyl.
  if(br.size() == 0)
    return -1;
  
  //std::cout << "Got here, branche " << br_ind << std::endl;
  //Parent cyl of the 1st cyl of the branch

  int par_cyl = cyls[cc].parent;
  if(par_cyl == -1)
    return -1;

  for(int i = 0; i < br.size(); i++){
    for(int j = 0; j < br[i].cyl_ind.size(); j++){
      if(par_cyl == br[i].cyl_ind[j])
	return i;
    }
  }
  
  return -2;//should not be possible
}

int scatter_output(std::vector<CylData> & cyls, int nCyl, std::vector<Branch> & br)
{
  std::ofstream scat;
  scat.open("scatter.dat",std::ofstream::out);
  std::vector<int> curr_brs;
  int order = 0;
  float len,rad,ang,gamma,zeta,tot_len;
  int n,m,k;
  V3f rax,newZ,newY,newX,newV;
  V3f newVprojXY,newVprojXZ;

  //TAPER
  find_branch_by_order(br,order,curr_brs);
  //std::cout << curr_brs.size() << std::endl;

  while(curr_brs.size() > 0){
    scat << "# order " << order << std::endl;
    // TAPER: Radius along a branch = F(length along the branch)
    scat << "# taper" << std::endl;
    for(int i = 0; i < curr_brs.size(); i++){
      len = 0.0;
      for(int j = 0; j < br[curr_brs[i]].cyl_ind.size(); j++){
	len += cyls[br[curr_brs[i]].cyl_ind[j]].length;
	rad = cyls[br[curr_brs[i]].cyl_ind[j]].radius;
	if(isnan(len) || isnan(rad)){
	  std::cerr << "Taper: instance of NaN" << std::endl;
	}
	else{
	  scat << len << " " << rad << std::endl;
	}
      }
    }
    scat << std::endl << std::endl;
    // BRA: branching angle of a branch
    scat << "# bra" << std::endl;
    if(order > 0){
      for(int i = 0; i < curr_brs.size(); i++){
	n = br[curr_brs[i]].cyl_ind[0];//1st cyl of the br
	ang = acos(cyls[n].axis*cyls[cyls[n].parent].axis);
	if(isnan(ang)){
	  std::cerr << "Bra: instance of NaN: " << std::endl;
	  std::cerr << cyls[n].axis << "; " << cyls[cyls[n].parent].axis <<\
	    ": " << cyls[n].axis*cyls[cyls[n].parent].axis << std::endl;
	}
	else{
	  scat << ang << std::endl;
	}
      }
    }
    scat << std::endl << std::endl;
    // LCHI_LAPAR: Length of the child branch = F(length along the parent branch)
    scat << "# lchi_lapar" << std::endl;
    if(order > 0){
      for(int i = 0; i < curr_brs.size(); i++){
	n = br[curr_brs[i]].cyl_ind[0];//1st cyl of the br
	m = parent_branch(cyls,br,n);//parent branch of the branch i
	k = cyls[n].parent;//parent cyl of n
	len = 0.0;
	for(int j = 0; j < br[m].cyl_ind.size(); j++){
	  len += cyls[br[m].cyl_ind[j]].length;
	  if(br[m].cyl_ind[j] == k)
	    break;
	}
	tot_len = 0.0;
	for(int j = 0; j < br[curr_brs[i]].cyl_ind.size(); j++){
	  tot_len += cyls[br[curr_brs[i]].cyl_ind[j]].length;
	}
	if(isnan(len) || isnan(tot_len)){
	  std::cerr << "lchi_lapar: instance of NaN " << std::endl;
	}
	else{
	  scat << len << " " << tot_len << std::endl;
	}
      }
    }
    scat << std::endl << std::endl;
    // LCHI_BRA_LAPAR: LCHI+BRA = F(LAPAR)
    scat << "# lchi_bra_lapar" << std::endl;
    if(order > 0){
      for(int i = 0; i < curr_brs.size(); i++){
	n = br[curr_brs[i]].cyl_ind[0];//1st cyl of the br
	m = parent_branch(cyls,br,n);//parent branch of the branch i
	k = cyls[n].parent;//parent cyl of n
	len = 0.0;
	for(int j = 0; j < br[m].cyl_ind.size(); j++){
	  len += cyls[br[m].cyl_ind[j]].length;
	  if(br[m].cyl_ind[j] == k)
	    break;
	}
	tot_len = 0.0;
	for(int j = 0; j < br[curr_brs[i]].cyl_ind.size(); j++){
	  tot_len += cyls[br[curr_brs[i]].cyl_ind[j]].length;
	}
	ang = acos(cyls[n].axis*cyls[cyls[n].parent].axis);
	if(isnan(len) || isnan(tot_len) || isnan(ang)){
	  std::cerr << "lchi_bra_lapar: instance of NaN " << std::endl;
	}
	else{
	  scat << len << " " << tot_len << " " << ang << std::endl;
	}
      }
    }
    scat << std::endl << std::endl;
    // CURV: 2-angle orientation in space = F(length along the branch)
    scat << "# curv" << std::endl;
    for(int i = 0; i < curr_brs.size(); i++){
      len = 0.0;
      tot_len = 0.0;
      for(int j = 0; j < br[curr_brs[i]].cyl_ind.size(); j++)
	tot_len += cyls[br[curr_brs[i]].cyl_ind[j]].length;
      for(int j = 1; j < br[curr_brs[i]].cyl_ind.size(); j++){
	//NOTE: X,Y,Z-axes are X,-Z,Y-axes in LPFG
	//NOTE: newX,newY,newZ-axes are new basis for X,Y,Z (in normal orientation).

	//****************************************************************
	// I. GLOBAL ADJUSTMENT OF THE COORD SYSTEM
	// newX = cyls[br[curr_brs[i]].cyl_ind[j-1]].axis;
	// newY = V3f(0.0,1.0,0.0);
	// newZ = V3f(0.0,0.0,1.0);
	// newV = cyls[br[curr_brs[i]].cyl_ind[j]].axis;
	// rax =  V3f(1.0,0.0,0.0) % newX;//rot. axis
	// rax.Normalize();
	// //************************************************************************
	// //*************************** IMPORTANT NOTE *****************************
	// //************************************************************************
	// // newX, newY and newZ form the new basis
	// // We need to define the projections on newXnewY and newXnewZ planes.
	// // Thus, we need to express every vector in the new basis, i.e. we need to
	// // solve the next equation to find Vnew (see below):
	// // newV = M * oldV, where oldV is the vector with coordinates in the original
	// // basis, whereas newV is that in the new basis. M is the inverse of the
	// // transformation matrix from old to new basis (in our case just rotation
	// // matrix, by which the transformation was made). We can note that M = R^(-1)
	// // = R', where R' is the rotation matrix obtained by inverting the angle of
	// // rotation: from Ang to -Ang.
	// // Thus, newV = rotate_v3f(oldV,rax,-ang);
	// // See for example here:
	// // https://www.math.hmc.edu/calculus/tutorials/changebasis/
	// //************************************************************************
	// if((rax.Length()) > 1e-05){
	//   ang = acos(newX * V3f(1.0,0.0,0.0));//Find the angle of rotation
	//   newY = rotate_v3f(newY,rax,ang);
	//   newZ = rotate_v3f(newZ,rax,ang);
	//   newV = rotate_v3f(newV,rax,-ang);//newV now expressed in the new basis
	// }
	// //Define the projections in the original basis
	// newVprojXY = V3f(newV.x,newV.y,0.0);//projection
	// newVprojXY = newVprojXY.x*newX + newVprojXY.y*newY;//in original basis
	// newVprojXY.Normalize();//normalize
	// newVprojXZ = V3f(newV.x,0.0,newV.z);
	// newVprojXZ = newVprojXZ.x*newX + newVprojXZ.z*newZ;
	// newVprojXZ.Normalize();
	// //Convert from the LPFG orientation to the normal, where Z- is UP direction
	// newX = lpfg_to_normal_orientation(newX);
	// newY = lpfg_to_normal_orientation(newY);
	// newZ = lpfg_to_normal_orientation(newZ);
	// newVprojXY = lpfg_to_normal_orientation(newVprojXY);
	// newVprojXZ = lpfg_to_normal_orientation(newVprojXZ);
	// //Gamma angle: XY-component (X-Z in LPFG)
	// gamma = atan2f(newVprojXY * newY, newVprojXY * newX);
	// //Zeta angle: XZ-component (XY in LPFG)
	// zeta = atan2f(newVprojXZ * newZ, newVprojXZ * newX);
	// END 
	//**************************************************************
	//
	// II. TWO-ROTATIONS METHOD
	// 1. Rotate around (Z) axis for X-axis to coincide with the j-1'th cyl(parent)
	// projection onto XY-plane. Rotate the Y-axis correspondingly and
	// do not change Z-axis. Calculate gamma in the newX-newY plane for j'th cyl(child).
	// 2. Rotate around (Y) axis for newX to coincide with the j-1'th cyl. Rotate
	// the Z-axis correspondingly and do not change Y-axis. Calculate zeta in the
	// newX-newZ plane for j'th cyl(child).
	// NOTE: X = X(lpfg), Y = -Z(lpfg), and Z = Y(lpfg)

	//1. Rotation around Z and XY projection angle gamma
	modify_coord_xy(cyls[br[curr_brs[i]].cyl_ind[j-1]].axis,newX,newY);
	gamma = projection_angle(cyls[br[curr_brs[i]].cyl_ind[j]].axis,newX,newY);
	
	//2. Rotation around Y and XZ projection angle zeta
	modify_coord_xz(cyls[br[curr_brs[i]].cyl_ind[j-1]].axis,	\
			newX,newZ,(const V3f)newY);
	zeta = projection_angle(cyls[br[curr_brs[i]].cyl_ind[j]].axis,newX,newZ);

	//**************************************************************

	//Relative length along the branch
	len += cyls[br[curr_brs[i]].cyl_ind[j-1]].length / tot_len;
	//Output
	if(isnan(gamma) || isnan(zeta) || isnan(len)){
	  std::cerr << "curv: instance of NaN" << std::endl;
	}
	else{
	  scat << gamma << " " << zeta << " " << len << std::endl;
	}
      }
    }
    scat << std::endl << std::endl;
    
    //+++++++++++++++++++++++++
    order++;
    find_branch_by_order(br,order,curr_brs);
  }

  scat.close();
  return 0;
}

int extract_branches(std::vector<CylData> & cyls, int nCyl, std::vector<Branch> & br)
{
  if(nCyl == 0)
    return 1;

  br.resize(0);
  
  int cc = 0;
  std::vector<int> base_cyls(1,cc);

  Branch br_tmp;
  int n_br = 0;
  
  while(base_cyls.size() > 0){
    cc = base_cyls[0];//Rename the first cyl in BASE_CYLS
    br_tmp.cyl_ind.resize(0);
    br_tmp.cyl_ind.push_back(cc);

    while(cyls[cc].extension > -1){
      for(int i=0;i<cyls[cc].children.size();i++)
	base_cyls.push_back(cyls[cc].children[i]);
      cc = cyls[cc].extension;
      br_tmp.cyl_ind.push_back(cc);
    }
    for(int i=0;i<cyls[cc].children.size();i++)
      base_cyls.push_back(cyls[cc].children[i]);

    base_cyls.erase(base_cyls.begin());

    //Copy to the branch array
    br_tmp.order = cyls[br_tmp.cyl_ind[0]].order;
    br_tmp.parent_br = parent_branch(cyls,br,br_tmp.cyl_ind[0]);
    br.push_back(br_tmp);
    n_br++;
  }

  // for(int i = 0; i < br.size(); i++){
  //   std::cout << "Branch " << i << ":" << std::endl;
  //   std::cout << br[i];
  // }

  return 0;
}
