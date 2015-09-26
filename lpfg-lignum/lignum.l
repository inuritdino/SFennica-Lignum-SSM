#include "lgm.h"

//====== AXIOM the initial seed for the growing tree 
derivation length: YEARS;
axiom: Cyl(0) Bud();

//====== CONTROL STATEMENTS
Start:
{
  //Get some space for cyl's
  cyls.resize((++cyl_resize_count)*MAXNCYL);
  //Init timer
  startT = clock();
  //Init the axiom cyl
  cyls[0].setAgeOrder(0,0);
  cyls[0].setLength(INITLEN);
  Wr += AR * cyls[0].Wf;
  cyls[0].setAxis(sta,Dir.Normalize());
  cyls[0].setParent(-1);//no parent for the axiom
  cyls[0].setExtension(-1);//init non-valid extension for the axiom
  cyls[0].is_deleted = false;
  //Init iterators
  nIter = 0;
  Year = 0;
  //Init rand generator
  //std::cout << "RNDSEED = " << RNDSEED << std::endl;
  if(RNDSEED)//Fixed random seed
    sran(RNDSEED);
  else{//if RNDSEED == 0, that means automatic(based on time) seeding
  time_t timer;
  long int x = time(&timer);
  if(x < LONG_MAX){
    sran(x);
    //std::cout << "Init sran with " << x << std::endl;
  }
  else{
    sran(LONG_MAX - (x % LONG_MAX));
    //std::cout << "Init sran by " << (LONG_MAX - (x % LONG_MAX)) << std::endl;
  }
  }

}
End:
{
  std::cout << "%%% Simulation finished. Year " << Year << " %%%" << std::endl;
  std::cout << "Total number of cyl's: " << nCyl <<  \
    "/" << num_deleted(nCyl,cyls) << " deleted (";
  std::cout << "by space: " << nRemovedSpace << ", shedding: " << \
    nRemovedShed << ", mass: " << nRemovedMass <<		    \
    ", rupture: " << nRemovedRupt << ")" << std::endl;
  std::cout << "Not grown due to environment: " << nRemovedEnvTot << std::endl;
  std::cout << "Max rupture: " << maxRupt << std::endl;
  // for(int i=0; i < nCyl; i++)
  //   std::cout << cyls[i];
  //std::cout << "Wr = " << Wr << std::endl;
  std::cout << "Diameter at base: " << 2*cyls[0].radius << std::endl;
  std::cout << "Breast height(1.3m) diameter: " << breast_height_diam(cyls,nCyl,sta) \
	    << std::endl;
  std::cout << "Height: " << tree_height(cyls,nCyl,sta.y) << std::endl;
  if(futileYears.size())
    std::cout << "Futile years: " << std::endl;
  for(int i=0; i<futileYears.size(); i++)
    std::cout << futileYears[i] << " ";
  if(futileYears.size())
    std::cout << std::endl;
  //Stop timer and report
  stopT = clock();
  std::cout << "Calculation time: " << \
    ((double)stopT - (double)startT)/(double)CLOCKS_PER_SEC << " sec." << std::endl;
  if(cyl_resize_count > 1)
    std::cout << "CylData resized " << cyl_resize_count-1 << " times." << std::endl;

  //Write the output structure in MTG file
  mtg_write(cyls,nCyl);
  //Random number generator check
  // float x[1000];
  // std::ofstream out;
  // out.open("rangauss.dat",std::ofstream::out);
  // for(int i=0; i<1000;i++){
  //   x[i] = ran_gauss(0.0,4.0);
  //   out << x[i] << std::endl;
  // }
  // out.close();
  std::vector<Branch> br;
  extract_branches(cyls,nCyl,br);
  scatter_output(cyls,nCyl,br);
  std::cout << "Extracted " << br.size() << " branches" << std::endl;
  //LOG file
  logfile.open("lignum.log",std::ofstream::out);
  logfile << fatal_error << std::endl;
  logfile.close();
  //Overall info to the file
  overall_info(Year,nCyl,cyls,sta,maxRupt,br);

}
StartEach:
{
  if(fatal_error)
    return ;
  nIter++;
  // Age increment: adds a year to the age of cyl's grown in previous years
  age_increment(nCyl);
  //Increment YEARs
  Year++;
  if(VERB)
    std::cout << "Year: " << Year << std::endl;
  // PHYSIOLOGY BALANCE CALCULATION
  //Respiration
  M = MF*tot_fol_mass(nCyl) + MR*Wr + MS*tot_sap_mass(nCyl);
  // Radiation conditions
  float ip[nCyl];
  float I[nCyl];//insolation energy (absolute) for each cyl
  fL.resize(nCyl);
  shadow_preset(vGrid,vb,sta,cyls,nCyl);//Pre set the shadow pattern
  propagate_shadow(I,ip,cyls,nCyl,vb,Year,vGrid,fatal_error);
  fl_ip(ip,fL,nCyl);
  // Photosynthesis
  P = photosynthesis(I,nCyl);
  if(VERB){
    std::cout << "P = " << P << ", M = " << M;
    std::cout << ", P - M = " << P-M << std::endl;
  }
}
EndEach:
{
  if(Year == YEARS || fatal_error){//THe last year simulation, do not produce anything
    return;
  }
  if(futileYears.size() != 0 && futileYears.back() == Year){
    //Perform senescence
    for(int i=0; i<nCyl; i++){
      if(cyls[i].is_deleted)
	continue;
      //Foliage shedding(this was accounted when calculated Lambda)
      cyls[i].Wf = pf_emp(i) * cyls[i].Wf0;
      //Sapwood senescence: sapwood turns into heartwood
      if(cyls[i].age > SSSTART){
	cyls[i].hwR = sqrt(cyls[i].hwR*cyls[i].hwR +			\
			   SS*(cyls[i].radius*cyls[i].radius-cyls[i].hwR*cyls[i].hwR));
      }
    }
    //Root senescence
    Wr -= SR * Wr;
    return;
  }
  if(VERB && nRemovedEnv)
    std::cout << "Removed " << nRemovedEnv << " buds due to environment" << std::endl;
  nRemovedEnv = 0;
  //*************
  //MAIN EQUATION:
  // A[0] + A[1]*Lambda + A[2]*Lambda^2 + A[3]*Lambda^3 = 0
  // The new cyl's(shoots) mass increments
  for(int i=nCyl-nShoot; i<nCyl; i++){
    A[2] += AR*cyls[i].Wf;//Root mass increment(Lambda^2)
    A[2] += cyls[i].Wf;//Foliage mass increment(Labmda^2)
    //Wood mass increment
    A[3] += M_PI*cyls[i].radius*cyls[i].radius*cyls[i].length*RHO;//Lambda^3 term
  }
  // (P-M) term goes to A[0] with minus
  A[0] += M-P;

  //**************
  //PIPE MODEL
  //*************
  std::vector<float> sa_final (nCyl-nShoot);
  std::vector<float> sa_current (nCyl-nShoot);
  pipe_model(cyls,nCyl,nShoot,A,sa_final,sa_current);

  //********************
  //SOLVE THE EQUATION
  //********************
  int lambda_iter = 0;
  float x0;
  // std::cout << "Equation: " << A[0] << " + " << A[1] << "*L" << " + " << A[2] <<\
  //   "*L^2" << " + " << A[3] << "*L^3" << " = 0" << std::endl;
  do{
    x0 = ran(5.0);
    //Lambda = my_zero(0.0,1.0-lambda_iter*0.05,0.0000001,A);
    Lambda = my_zero(0.0,x0,0.0001,A);
    lambda_iter++;
    if(lambda_iter == 20){
      break;
    }
  }
  while((fabs(A[0]+A[1]*Lambda+A[2]*Lambda*Lambda+A[3]*Lambda*Lambda*Lambda) \
	 > 0.0001) || Lambda < 0.0);
  lambda_iter = 0;
  while((fabs(A[0]+A[1]*Lambda+A[2]*Lambda*Lambda+A[3]*Lambda*Lambda*Lambda) \
	 > 0.0001) || Lambda < 0.0){
    //Lambda = my_zero(0.0,1.0+lambda_iter*0.05,0.0000001,A);
    x0 = ran(10.0);
    Lambda = my_zero(0.0,x0,0.0001,A);
    lambda_iter++;
    if(lambda_iter == 20){
      std::cout << "WARNING!" << std::endl;
      std::cout << "\tThe equation: ";
      std::cout << A[0] << " + " << A[1] << "*L" << " + " << A[2] << "*L^2";
      std::cout << " + " << A[3] << "*L^3" << " = " << 0 << std::endl;
      std::cout << "\tLambda = " << Lambda << ": ";
      std::cout << A[0]+A[1]*Lambda+A[2]*Lambda*Lambda+A[3]*Lambda*Lambda*Lambda	\
		<< " = 0.0 (after " << 2*lambda_iter << " iterations)." << std::endl;
      fatal_error = 1;
      break;
    }
  }

  if(VERB)
    std::cout << "Lambda = " << Lambda << std::endl;

  //Branch bending, should be done before removing cyl's
  branch_bending_wrapper(cyls);

  //Cyls to remove later on
  std::vector<int> cyls_to_remove(0);

  //Remove the old cyl's that traverse the ground and their offsets to avoid gaps
  std::vector<int> offsets(0);
  bool duplicate = false;
  for(int i = 0; i < nCyl-nShoot; i++){
    if((cyls[i].end.y < 0.0) && (!cyls[i].is_deleted)){
      get_all_offsets(i,offsets);
      for(int j = 0; j < offsets.size(); j++){
	if(!do_belong_to_array(cyls_to_remove,offsets[j])){
	  cyls_to_remove.push_back(offsets[j]);
	}
      }
      offsets.resize(0);
    }
  }
  
  //CORRECT the shoot length and radii according to the found Lambda
  int k;
  for(int i=nCyl-nShoot; i<nCyl; i++){
    //The final length for the cyl
    cyls[i].setLength(cyls[i].length * Lambda);
    cyls[i].setAxis(cyls[i].start,cyls[i].axis);//redefine the end point
    // Check the collision with other cyl's in space
    if(traversed_in_space(i)){//if true, do not grow
      if(!do_belong_to_array(cyls_to_remove,i))
	cyls_to_remove.push_back(i);
      continue;
    }
    if(cyls[i].length < MINLEN){
      if(!do_belong_to_array(cyls_to_remove,i))
	cyls_to_remove.push_back(i);
      continue;
    }
    //Increment root mass
    Wr += AR * cyls[i].Wf;
  }

  //******************
  //SENESCENCE
  //******************
  std::vector<bool> pipe_failed (nCyl-nShoot,false);
  //ADJUST the radius of older cyl's according to the found Lambda
  float tmp_rad;
  for(int i=0; i<nCyl-nShoot; i++){
    if(cyls[i].is_deleted)
      continue;
    //Sapwood senescence: sapwood turns into heartwood
    if(cyls[i].age > SSSTART){
      cyls[i].hwR = sqrt(cyls[i].hwR*cyls[i].hwR +			\
			 SS*(cyls[i].radius*cyls[i].radius-cyls[i].hwR*cyls[i].hwR));
    }

    //New radius
    tmp_rad = cyls[i].radius;//keep the current radius value
    if( (Lambda*Lambda*sa_final[i]/M_PI - sa_current[i]/M_PI + \
	 sap_area(i)/M_PI + sap_area_fol(i)/M_PI + cyls[i].hwR*cyls[i].hwR) >= 0.0){
      cyls[i].radius = \
	sqrt(\
	     Lambda*Lambda*sa_final[i]/M_PI - sa_current[i]/M_PI \
	     + sap_area(i)/M_PI	\
	     + sap_area_fol(i)/M_PI		\
	     + cyls[i].hwR*cyls[i].hwR\
	     );
      if(cyls[i].radius < tmp_rad){
	pipe_failed[i] = true;
	cyls[i].radius = tmp_rad;
      }
    }
    if(std::isnan(cyls[i].radius)){
      std::cout << "ERROR: secondary thickening with NaN radius!" << std::endl;
    }
    //Foliage shedding
    cyls[i].Wf = pf_emp(i) * cyls[i].Wf0;
  }
  //Root senescence
  Wr -= SR * Wr;

  //Debugging the pipe model conditions
  //pipe_model_debug(cyls,nCyl,pipe_failed);

  //Delete cyl's not approved by space and length limitations
  delete_cyls(cyls,cyls_to_remove,nCyl);
  nRemovedSpace += cyls_to_remove.size();
  
  //Branch shedding(no foliage in the branch ==> shed it)
  branch_shedding(cyls,nCyl,nRemovedShed);
  
  //Shedding by mass
  //branch_shed_by_mass();
  //branch_shed_by_rupture();


  //Nullify the count for the new shoots and equation coefficients
  nShoot = 0;
  A[0] = 0.0;
  A[1] = 0.0;
  A[2] = 0.0;
  A[3] = 0.0;

  //Print deleted cyls
  //print_deleted_cyls(nCyl);
}

//===== PRODUCTIONS

Cyl(parCyl) < Bud()://Shoots production
{
  if(Year == YEARS || fatal_error){//THe last year simulation
    return;
  }
  if((P-M) <= 0.0){
    std::cout << "FATAL ERROR: Balance is negative." << std::endl;
    fatal_error = 1;
    // if(futileYears.size() == 0 || futileYears.back() != Year){
    //   std::cout << "Balance is negative. Futile year " << Year << std::endl;
    //   futileYears.push_back(Year);
    // }
    return;
  }
  if(cyls[parCyl].is_deleted)//parent cyl is deleted do not produce
    return;
  if( sqrt((cyls[parCyl].end.x - sta.x)*(cyls[parCyl].end.x - sta.x) +	\
	   (cyls[parCyl].end.z - sta.z)*(cyls[parCyl].end.z - sta.z)) > \
      ((float)ENVSHRAD) ){
    //std::cout << "Cyl " << parCyl << " bud does not grow " << std::endl;
    nRemovedEnv++;
    nRemovedEnvTot++;
    produce Cyl(parCyl);
  }
  // if(!branch_grow_lr(parCyl))//Len/BaseRad not approved.
  //   return;

  //std::cout << "fL parent = " << fL[parCyl] << std::endl;
  
  //Check whether to produce the buds
  bool remove_l, remove_t;
  remove_buds(remove_l,remove_t,parCyl);

  //Determine number of L-buds
  int lbuds = 0;
  if(!remove_l)
    lbuds = num_L_buds(parCyl);

  //Produce the L-buds first
  float rollAng;
  if(lbuds > 0){
    rollAng = 360.0/(float)lbuds;
    // if(cyls[parCyl].order == 1)
    //   rollAng = 180.0/(float)lbuds;
  }

  //Axis stuff
  V3f cyl_axis,rdir,rot_axis;
  static V3f init_lat_dir_w1 = V3f(1.0,0.0,0.0);
  
  for(int i=0; i<lbuds; i++){
    cyls[nCyl].setAgeOrder(0,cyls[parCyl].order+1);//set age and order
    cyls[nCyl].is_deleted = false;
    cyls[nCyl].setParent(parCyl);//this is a child of parCyl, since L-bud
    cyls[nCyl].setExtension(-1);//init the non-valid extension
    cyls[parCyl].setChild(nCyl);// add child to the parent
    //IMPORTANT: here length, radius and related entities are set prefactored with
    //Lambda to be found later, e.g., if radius = R, the actual radius will be reset
    //later to R*Lambda. The same for the T-buds expansion below.
    cyls[nCyl].setLength(fL[parCyl]*(1 - (cyls[parCyl].order+1)*Q));//ORDER + 1?
    //Rotate to get axis
    if(i == 0){
      if(cyls[parCyl].order == 0){
	//std::cout << "Init. lateral dir. before: " << init_lat_dir_w1 << std::endl;
	rdir = init_lat_dir_w1;
	init_lat_dir_w1 = \
	  rotate_v3f(init_lat_dir_w1,cyls[parCyl].axis,deg_to_rad(INITROTANG));
	//std::cout << "Init. lateral dir. after: " << init_lat_dir_w1 << std::endl;
      }
      else{
	//rdir = rnd_dir(cyls[parCyl].axis);
	rot_axis = V3f(0.0,1.0,0.0);//temporarily use rot_axis
	// std::cout << "Scalar: " << rot_axis << " * " << cyls[parCyl].axis \
	// 	  << " = " << rot_axis * cyls[parCyl].axis << std::endl;
	if(fabs((rot_axis * cyls[parCyl].axis) - 1.0) < 1e-05)
	  rot_axis = V3f(1.0,0.0,0.0);
	rdir = rotate_v3f(cyls[parCyl].axis,rot_axis,1.0);
	// std::cout << "Rotate y around y: " << \
	//   rotate_v3f(rot_axis,rot_axis,1.0) << std::endl;
      }
      rot_axis = cyls[parCyl].axis % rdir;
      cyl_axis = rotate_v3f(cyls[parCyl].axis,rot_axis,deg_to_rad(BRAANG));
    }
    else{
      rot_axis = cyls[parCyl].axis;
      cyl_axis = rotate_v3f(cyl_axis,cyls[parCyl].axis,deg_to_rad(rollAng));
    }
    if(fabs(cyl_axis.Length()-1) > 1e-04)
      std::cout << "ERROR: New L-shoot length = " << cyl_axis.Length() << std::endl;
    
    cyls[nCyl].setAxis(cyls[parCyl].end,cyl_axis);
    nCyl++;
    nShoot++;
    check_num_cyls(nCyl);
    nproduce SB Cyl(nCyl-1) Bud() EB;
  }
  if(remove_t){
    //produce ;//produce nothing, this also returns from the function
    produce Bud;//OR produce just Bud, delayed production
  }

  //Then, produce the T-bud
  cyls[nCyl].setAgeOrder(0,cyls[parCyl].order);//set age and order
  cyls[nCyl].is_deleted = false;
  cyls[nCyl].setParent(parCyl);//this is extension of the parCyl
  cyls[nCyl].setExtension(-1);//init the non-valid extension
  cyls[parCyl].setExtension(nCyl);//set the extension of the parCyl
  cyls[nCyl].setLength(fL[parCyl]*(1 - cyls[parCyl].order*Q));
  if(fabs(cyls[parCyl].axis.Length()-1) > 1e-04){
    std::cout << "ERROR: New T-shoot length = " << cyls[parCyl].axis.Length() << \
      std::endl;
  }
  V3f newX,newY,newZ;
  float ang;
  cyl_axis = cyls[parCyl].axis;
  //This random rotation introduces some variability to ZETA as well
  if(cyls[parCyl].order > 0){
    modify_coord_xy(cyl_axis,newX,newY);
    rot_axis = newX % newY;
    ang = ran_gauss_any(0.0,GAMMASD);
    cyl_axis = rotate_v3f(cyl_axis, rot_axis, deg_to_rad(ang));
    cyl_axis.Normalize();
  }
  cyls[nCyl].setAxis(cyls[parCyl].end,cyl_axis);
  nCyl++;
  nShoot++;
  check_num_cyls(nCyl);
  //std::cout << "Produced " << nShoot << " shoots" << std::endl;
  produce Cyl(nCyl-1) Bud;
}

//===== INTERPRETATIONS

interpretation:
// Bud():
// {
//   produce SetColor(1) Sphere0;
// }

Cyl(cyl):
{
  if((cyl == 0) && LABELFLAG){
    float height = 0.0;
    int highest_cyl = 0;
    for(int i=0;i<nCyl;i++){
      if(cyls[i].is_deleted)
    	continue;
      if(height < cyls[i].end.y){
    	height = cyls[i].end.y;
    	highest_cyl = i;
      }
    }
    height = height - sta.y;
    std::string h_label("Height = ");
    h_label.append(std::to_string(height));
    std::string d_label("Base diameter = ");
    d_label.append(std::to_string(2*cyls[0].radius));
    nproduce SetColor(3) MoveTo(cyls[highest_cyl].end.x,1.01*cyls[highest_cyl].end.y, \
    		    cyls[highest_cyl].end.z) \
      Label(h_label.c_str());
    nproduce MoveTo(cyls[0].start.x,-0.05,cyls[0].start.z) \
      Label(d_label.c_str());
    nproduce MoveTo(1.0,0.0,1.0) SetWidth(0.03)	\
      StartGC() G(1.0) EndGC() Label("1.0m");
    //std::cout << "Base diameter: " << d_label << std::endl;
    // //Plot the neghboring pivots as shadowers
    // nproduce CurrentTexture(0) SetColor(2) MoveTo(0.1,0.0,0.5) SetHead(0.0,1.0,0.0,1.0,0.0,0.0) \
    //   SetWidth(cyls[0].radius) StartGC() G(height) EndGC();
  }
  if(!(cyls[cyl].is_deleted)){
    V3f rdir = V3f(0.0,1.0,0.0);
    if(fabs((rdir * cyls[cyl].axis) - 1.0) < TOL)
      rdir = V3f(1.0,0.0,0.0);
    V3f rot_axis = cyls[cyl].axis % rdir;
    rot_axis.Normalize();
    V3f up = rotate_v3f(cyls[cyl].axis,rot_axis,M_PI/2.0);//any perpendicular vector
    nproduce CurrentTexture(0) CurrentContour(0) SetColor(2) ContourSides(10) \
      MoveTo(cyls[cyl].start.x,cyls[cyl].start.y,cyls[cyl].start.z) \
      SetHead(cyls[cyl].axis.x,cyls[cyl].axis.y,cyls[cyl].axis.z,up.x,up.y,up.z) \
      SetWidth(2*cyls[cyl].radius) StartGC() G(cyls[cyl].length) EndGC();
    //************ NEEDLE DRAWING ****************
    if(FOLFLAG){
      V3f n_st,n_st0,n_he;
      //int divide = floor(NUMNEEDLELONG*cyls[cyl].Wf/cyls[cyl].Wf0);
      int divide = floor(cyls[cyl].length/(((float)NEEDLEDENMAX-(float)NEEDLEDENMIN)*(cyls[cyl].Wf/cyls[cyl].Wf0)));
      for(int j = 0; j < divide; j++){
	//Needle start position
	n_st0 = cyls[cyl].start + ((j+1.0)/(1.0+(double)divide))*cyls[cyl].length*cyls[cyl].axis;
	for (int i = 0; i < NUMNEEDLECIRC; i++){
	  //Adjust start position
	  n_st = n_st0 + cyls[cyl].radius*(rotate_v3f(cyls[cyl].axis,rot_axis,M_PI/2.0));
	  //Needle heading direction
	  n_he = rotate_v3f(cyls[cyl].axis,rot_axis,NEEDLEANG);
	  n_he.Normalize();

	  //Define up vector (reuse up)
	  rdir = V3f(0.0,1.0,0.0);
	  if(fabs((rdir * n_he) - 1.0) < TOL)
	    rdir = V3f(1.0,0.0,0.0);
	  up = n_he % rdir;
	  up.Normalize();
	  up = rotate_v3f(n_he,up,M_PI/2.0);
	  nproduce SetColor(1) MoveTo(n_st.x,n_st.y,n_st.z)			\
	    SetHead(n_he.x,n_he.y,n_he.z,up.x,up.y,up.z)			\
	    SetWidth(NEEDLEDIAM) G(NEEDLELEN);
	  rot_axis = rotate_v3f(rot_axis,cyls[cyl].axis,2.0*M_PI/(NUMNEEDLECIRC));
	}
      }
    }
    //**************************************************
    // DEBUGGING THE SHADOW PROPAGATION. VISUALIZATION.
    // V3f cyl_point;
    // int xv,yv,zv;
    // float x,y,z;
    // float axis_fraction = 1.0/(NPOINCYL+1.0);
    // std::cout << "Axis fraction = " << axis_fraction << std::endl;
    // for(int j=0;j<NPOINCYL;j++){
    //   cyl_point = cyls[cyl].start + (j+1)*axis_fraction*cyls[cyl].length*cyls[cyl].axis;
    //   //Define the voxels over X,Y and Z spans. 0-based numbering.
    //   //std::cout << "Cyl center = " << cyl_point << std::endl;
    //   //std::cout << "Offset = " << vb.offset << std::endl;
    //   xv = ceil( (cyl_point.x - vb.offset.x) / vb.vdx ) - 1;
    //   yv = ceil( (cyl_point.y - vb.offset.y) / vb.vdy ) - 1;
    //   zv = ceil( (cyl_point.z - vb.offset.z) / vb.vdz ) - 1;
    //   //std::cout << xv << "," << yv << "," << zv << std::endl;
    //   for(int iy = (yv); iy >= (yv-QMAX); iy--){
    // 	for(int ix = (xv-(yv-iy)); ix <= (xv+(yv-iy)); ix++){
    // 	  for(int iz = (zv-(yv-iy)); iz <= (zv+(yv-iy)); iz++){
    // 	    //std::cout << "Voxels: " << ix << "," << iy << "," << iz << std::endl;
    // 	    if( (ix < 0) || (ix >= GSX) || (iy < 0) || \
    // 		(iy >= GSY) || (iz < 0) || (iz >= GSZ) )
    // 	      continue;
    // 	    x = vb.offset.x + ix*VSX + VSX/2;
    // 	    y = vb.offset.y + iy*VSY + VSY/2;
    // 	    z = vb.offset.z + iz*VSZ + VSZ/2;
    // 	    //std::cout << "(x,y,z)=" << x << "," << y << "," << z << std::endl;
    // 	    nproduce MoveTo(x,y,z) Sphere(VSX/2.0);
    // 	  }
    // 	}
    //   }
    // }
    produce ;
  }
  else
    produce ;
}

// ***************************
//===== FUNCTION SECTION =====

float photosynthesis(float I[], int nCyl)
{
  float P = 0.0;
  for(int i=0; i<nCyl; i++){
    if(cyls[i].is_deleted)
      continue;
    if(cyls[i].Wf > TOL)
      P += PR * I[i];
  }
  return P;
}

int remove_buds(bool& remove_L, bool& remove_T, int cyl)
{//Determines whether to remove L- and/or T-buds for the parent cyl
  //NOTE: (1-q*order) < 0 implies also (1-q*(order+1)) < 0.
  if( (fL[cyl] < 1e-05) ){ // <= condition
    //std::cout << "fL condition success: " << fL[cyl] << ", cyl " << cyl << std::endl;
    remove_L = true;
    remove_T = true;
  }
  else if( (1 - Q*cyls[cyl].order) < 1e-05 ){
    //std::cout << "Order condition success: " << (1 - Q*cyls[cyl].order) << ", cyl " << cyl << std::endl;
    remove_L = true;
    remove_T = true;
  }
  else if( (1 - Q*(cyls[cyl].order+1)) < 1e-05){//and fL > 0
    //std::cout << "Order condition success for L-buds: " << (1 - Q*(cyls[cyl].order+1)) << ", cyl " << cyl << std::endl;
    remove_L = true;
    remove_T = false;
  }
  else{
    //std::cout << "All buds survive" << ", cyl " << cyl << std::endl;
    remove_L = false;
    remove_T = false;
  }
  return 0;
}

int fl_ip(float ip[],std::vector<float> &fL,const int size)
{
  for(int i=0; i<size; i++){
    fL[i] = (ip[i] - FLA) / (1 - FLA);
    if(fL[i] < 0.0)
      fL[i] = 0.0;
  }
  return 0;
}

float tot_sap_mass(int nCyl)
{
  float x = 0.0;
  for(int i=0; i<nCyl; i++){
    if(cyls[i].is_deleted)
      continue;
    x += RHO * cyls[i].length * M_PI * \
      (cyls[i].radius*cyls[i].radius - cyls[i].hwR*cyls[i].hwR);
    if(std::isnan(x)){
      std::cout << "ERROR: NaN for cyl " << i << ", " << cyls[i].length << ",";
      std::cout << cyls[i].radius << "," << cyls[i].hwR << "," << cyls[i].lr << std::endl;
      break;
    }
  }
  return x;
}

float tot_fol_mass(int nCyl)
{
  float x = 0.0;
  for(int i=0; i<nCyl; i++){
    if(cyls[i].is_deleted)
      continue;
    x += cyls[i].Wf;
  }
  return x;
}

float pf_instant(int cyl)
{//How much foliage is left compared to the current amount Wf
  if(cyls[cyl].age > 5)
    return 0.0;
  else if(cyls[cyl].age >= 4)
    return (0.3/0.5);
  else if(cyls[cyl].age >= 3)
    return (0.5/0.7);
  else if(cyls[cyl].age >= 2)
    return (0.7/0.9);
  else if(cyls[cyl].age >= 1)
    return 0.9;
  else
    return 1.0;
}

float pf_emp(int cyl)
{//How much foliage is left for a cyl compared to initial amount Wf0
  if(cyls[cyl].age >= 5)
    return 0.0;
  else if(cyls[cyl].age >= 4)
    return 0.3;
  else if(cyls[cyl].age >= 3)
    return 0.5;
  else if(cyls[cyl].age >= 2)
    return 0.7;
  else if(cyls[cyl].age >= 1)
    return 0.9;
  else
    return 1.0;
}

int age_increment(const int nCyl)
{
  for(int i=0; i<nCyl; i++)
    cyls[i].age++;
  return 0;
}

int num_L_buds(int c)
{
  //For debugging
  //nLbud = 1;
  
  if(fL[c] < (float)FLBUDS)
    return 0;
  int nLbud;
  
#if PERTTUNEN1998 != 0
  if(cyls[c].Wf > 0.06)
    nLbud = 7;
  else if(cyls[c].Wf > 0.03)
    nLbud = 6;
  else if(cyls[c].Wf > 0.01)
    nLbud = 5;
  else if(cyls[c].Wf > 0.005)
    nLbud = 4;
  else if(cyls[c].Wf > 0.002)
    nLbud = 3;
  else if(cyls[c].Wf > 0.0011)
    nLbud = 2;
  else if(cyls[c].Wf > 0.0002)
    nLbud = 1;
  else
    nLbud = 0;
#else //This includes SIEVANEN2008 conditions
  if(cyls[c].Wf > 0.0054)
    nLbud = 7;
  else if(cyls[c].Wf > 0.0042)
    nLbud = 6;
  else if(cyls[c].Wf > 0.003)
    nLbud = 5;
  else if(cyls[c].Wf > 0.0023)
    nLbud = 4;
  else if(cyls[c].Wf > 0.0017)
    nLbud = 3;
  else if(cyls[c].Wf > 0.0013)
    nLbud = 2;
  else if(cyls[c].Wf > 0.001)
    nLbud = 1;
  else
    nLbud = 0;
#endif
  
  return nLbud;
}

bool is_empty_array(std::vector<bool>& array)
{//Checks whether all elements in an array are false
  for(int i=0; i<array.size(); i++){
    if(array[i])
      return false;
  }
  return true;
}

float sap_area(int cyl)
{
  return (M_PI*(cyls[cyl].radius*cyls[cyl].radius-cyls[cyl].hwR*cyls[cyl].hwR));
}

float sap_area_fol(int i)
{
  return ((cyls[i].Wf-pf_emp(i)*cyls[i].Wf0)*(1-cyls[i].ksi)*cyls[i].lr/(2*cyls[i].af));
}

void check_num_cyls(int num_cyl)
{
  if(num_cyl >= (cyl_resize_count*MAXNCYL)){
    std::cout << "nCyl = " << num_cyl << ". Resizing the cyl structure." << std::endl;
    cyls.resize((++cyl_resize_count)*MAXNCYL);
  }
  
}

bool traversed_in_space(int cyl)
{//Is the CYL traversed in space with other cyl's if grown or goes below the ground
  //do not allow growing below the ground

  //Below the ground
  if(cyls[cyl].end.y < 0.0)
    return true;
  //Traversing in space with other cyl's
  float D,c1,c2;
  V3f x1,x2;
  Line ln(cyls[cyl].start,cyls[cyl].end-cyls[cyl].start);
  V3f cyl_center = 0.5*(cyls[cyl].start + cyls[cyl].end);
  bool child;
  for(int i=0; i<nCyl; i++){
    //Exclude deleted, cyl's parent, child/extension, itself, and cyls with common parent
    if(cyls[i].is_deleted || (i == cyl) || (i == cyls[cyl].parent) || (i == cyls[cyl].extension))
      continue;
    child = false;
    for(int j=0; j<cyls[cyl].children.size(); j++){
      if(i == cyls[cyl].children[j]){
	child = true;
	break;
      }
    }
    if(child)
      continue;
    //if( i >= nCyl-n_shoots ){//Shoots of the same year
      // std::cout << i << " vs. " << cyl << ": " << \
      // 	cyls[i].parent << " vs. " << cyls[cyl].parent << std::endl;
      if(cyls[i].parent == cyls[cyl].parent)
	continue;
      //}
    //Check the space only in vicinity of the current cyl. The vicinity is a sphere
    //around target cyl with center=center of cyl and radius=cyl's length.
    if( ((cyls[i].end-cyl_center).Length() > cyls[cyl].length) && \
	((cyls[i].start-cyl_center).Length() > cyls[cyl].length) )
      continue;
    //std::cout << "Check space: " << cyl << " vs. " << i << std::endl;
    D = min_dist_line_seg(x1,x2,c1,c2,ln,cyls[i].start,cyls[i].end);
    if( D < (cyls[cyl].radius+cyls[i].radius) && \
	(c1 >= 0.0) && (c1 <= 1.0) )
      return true;
  }

  return false;
}

bool branch_grow_lr(int tip)
{// Decides whether to grow a branch based on the Len/BaseRadius ratio. The TIP is
 // the tip cyl of the branch.
  float len = 0.0;
  std::vector<int> branch(1,tip);
  len += cyls[tip].length;
  int cc = tip;
  while((cyls[cc].parent > -1) && \
	(cyls[cyls[cc].parent].extension == cc)){
    branch.push_back(cyls[cc].parent);
    cc = cyls[cc].parent;
    len += cyls[cc].length;
  }

  float ratio = len/cyls[branch.back()].radius;
  if( ratio > BRCRITLR ){
    // std::cout << "Branch: " << std::endl;
    // for(int i=0; i<branch.size(); i++)
    //   std::cout << branch[i] << " ";
    // std::cout << "<" << ratio << ">" << " not grow" << std::endl;
    return false;
  }
  else{
    return true;
  }
  
}

int branch_shed_by_mass()
{//Shed the branches when they get critical mass.
  //Take the base cyl's of the 1st order branches
  std::vector<int> base_cyls(0);
  int i=0,j;
  //Extract 1st order base_cyls
  while(cyls[i].extension > -1){
    for(j=0; j < cyls[i].children.size(); j++)
      base_cyls.push_back(cyls[i].children[j]);
    i = cyls[i].extension;
  }
  for(j=0; j < cyls[i].children.size(); j++)
    base_cyls.push_back(cyls[i].children[j]);
  //std::cout << "Found " << base_cyls.size() << " 1st order branches" << std::endl;
  std::vector<int> branch(0);
  float len = 0.0;
  std::vector<int> offsets;
  for(j=0; j < base_cyls.size(); j++){
    //std::cout << "Branch #" << j+1 << std::endl;
    branch.push_back(base_cyls[j]);
    //len += (cyls[branch.back()].end - cyls[branch.back()].start).Length();
    len += cyls[branch.back()].length;
    while(cyls[branch.back()].extension > -1){
      branch.push_back(cyls[branch.back()].extension);
      //len += (cyls[branch.back()].end - cyls[branch.back()].start).Length();
      len += cyls[branch.back()].length;
    }
    //std::cout << "Base cyl: " << base_cyls[0] << std::endl;
    //std::cout << "Len/BaseRad ratio: " << len / cyls[branch[0]].radius << std::endl;
    // std::cout << "Pressure = " << 9.8*mass_of_cyls(branch) << "/" <<\
    //   cyls_ground_proj_area(branch) << " = " <<			\
    //   9.8*mass_of_cyls(branch)/cyls_ground_proj_area(branch) \
    // 	      << std::endl;
    //std::cout << "Mass center: " << mass_center(branch) << std::endl;
    // SIMPLE criterion for shedding by mass: Critical ratio Len/BaseCylRad.
    if((len/cyls[branch[0]].radius) > BRCRITLR){
      offsets.resize(0);
      get_all_offsets(branch[0],offsets);
      // std::cout << "Deleting branch rooted in " << branch[0]\
      // 		<< ". Ratio: " << len/cyls[branch[0]].radius << std::endl;
      delete_cyls(cyls,offsets,nCyl);
      nRemovedMass += offsets.size();
    }
    branch.resize(0);
    len = 0.0;
  }
  return 0;
}

int branch_shed_by_rupture()
{
  std::ofstream ofile;
  ofile.open("rupture.dat");
  std::vector<int> offsets;
  float rupt;
  for(int i = 0; i < nCyl; i++){
    if(cyls[i].is_deleted)
      continue;
    if(cyls[i].order != 0){
      rupt = rupture_modulus(i,offsets);
      if(rupt > maxRupt)
	maxRupt = rupt;
      ofile << rupt << std::endl;
      if(rupt > CRITRUPT){
	// std::cout << \
	//   "Critical rupture "<< CRITRUPT << " exceeded for cyl " <<	\
	//   i << ": R = " << rupt << std::endl;
	nRemovedRupt += offsets.size();
	delete_cyls(cyls,offsets,nCyl);
      }
      offsets.resize(0);
    }
  }
  ofile.close();
  return 0;
}

float rupture_modulus(int init_cyl,std::vector<int> & offsets)
{
  get_all_offsets(init_cyl,offsets);
  offsets.erase(offsets.begin());//remove the initializer cyl
  std::vector<int> arm(1,init_cyl);//arm-set of cyl's up to branch point
  int cc = init_cyl;
  while(cyls[cyls[cc].parent].extension == cc){
    arm.push_back(cyls[cc].parent);
    cc = cyls[cc].parent;
  }
  float arm_length = 0.0;

  //std::cout << "Init cyl " << init_cyl << ": ";
  for(int i=0; i<arm.size(); i++){
    arm_length += cyls[arm[i]].length;
    //std::cout << arm[i] << "(" << cyls[arm[i]].length << ")" << " ";
  }
  float A = cyls_ground_proj_area(arm);
  
  float Mtot = mass_of_cyls(offsets);
  float P;
  cc = init_cyl;
  arm_length *= 2;
  // P = 1.5 * Mtot * arm_length / \
  //   (2*cyls[cc].radius*4*cyls[cc].radius*cyls[cc].radius);
  P = 9.8 * Mtot; P = P / (2*A);
  //std::cout << "R = " << P << std::endl;
  return P;
}

V3f mass_center(std::vector<int> & my_cyls)
{
  V3f mass_center(0.0,0.0,0.0);
  V3f cyl_center;
  for(int i=0; i<my_cyls.size(); i++){
    cyl_center = cyls[my_cyls[i]].start + 0.5*cyls[my_cyls[i]].axis;
    if(i == 0)
      mass_center = mass_of_cyl(my_cyls[i]) * cyl_center;
    else
      mass_center += mass_of_cyl(my_cyls[i]) * cyl_center;
  }
  mass_center = mass_center / mass_of_cyls(my_cyls);
  return mass_center;
}

float mass_of_cyl(int c)
{//Mass of the wood and foliage of the cyl.
  float m = 0.0;
  m = M_PI*cyls[c].radius*cyls[c].radius*cyls[c].length*RHO + cyls[c].Wf;
  return m;
}

float mass_of_cyls(std::vector<int> & my_cyls)
{
  float M = 0.0;
  for(int i=0; i<my_cyls.size(); i++){
    M += mass_of_cyl(my_cyls[i]);
  }
  return M;
}

int get_all_offsets(int base_cyl,std::vector<int> &offsets)
{//Get all offset of the base cyl (initiator of a branch), including all children's
 //branches. Deleted cyls do not have topology information (i.e. ext's and child'n,
 //so they will not be in the offsets.
  // if((cyls[base_cyl].parent > -1) &&				\
  //    (cyls[cyls[base_cyl].parent].extension == base_cyl)){
  //   std::cerr << "ERROR(all_offsets): not a base cyl" << std::endl;
  //   return 1;
  // }
  if(cyls[base_cyl].is_deleted){
    std::cerr << "ERROR(all_offsets): deleted cyl" << std::endl;
    return 2;
  }
  int cc;
  std::vector<int> bcyls(1,base_cyl);
  while(bcyls.size()){
    cc = bcyls[0];//Rename the first cyl in BCYLS
    offsets.push_back(cc);//Add it to the offsets
    //Get all ext's of the CC
    while(cyls[cc].extension > -1){
      offsets.push_back(cyls[cc].extension);
      //Add base cyls from the children array of CC
      for(int i=0;i<cyls[cc].children.size();i++)
	bcyls.push_back(cyls[cc].children[i]);
      //Change CC
      cc = cyls[cc].extension;
    }
    //For the last extension
    for(int i=0;i<cyls[cc].children.size();i++)
      bcyls.push_back(cyls[cc].children[i]);

    bcyls.erase(bcyls.begin());
  }
  return 0;
}


float cyls_ground_proj_area(std::vector<int> & my_cyls)
{//Ground projection area. Ground is XZ-plane.
  float A = 0.0;
  V3f proj;//projection vector
  float cos_value;//cosine of the angle between projection and axis
  for(int i=0; i<my_cyls.size(); i++){
    proj = V3f(cyls[my_cyls[i]].axis.x,0.0,cyls[my_cyls[i]].axis.z);
    if( (fabs(proj.x) < 1e-05) &&	\
	(fabs(proj.y) < 1e-05) &&	\
	(fabs(proj.z) < 1e-05) )
      cos_value = 0.0;//projection is zero
    else{
      proj.Normalize();
      cos_value = proj * cyls[my_cyls[i]].axis;
    }
    if(cos_value < 0.0){
      std::cout << "COS = " << cos_value << std::endl;
    }
    A += 2*cyls[my_cyls[i]].radius*cyls[my_cyls[i]].length*cos_value + \
      M_PI*cyls[my_cyls[i]].radius*cyls[my_cyls[i]].radius*sqrt(1-cos_value*cos_value);
  }

  return A;
}



int num_deleted(int nCyl, std::vector<CylData>& cyls)
{
  int out = 0;
  for(int i=0; i<nCyl; i++)
    if(cyls[i].is_deleted)
      out++;

  return out;
}

int overall_info(int Year, int nCyl, std::vector<CylData> & cyls, V3f sta,\
		 float maxRupt,std::vector<Branch> & br)
{
  std::ofstream info;
  info.open("info.dat",std::ofstream::out);
  //Headers
  info << "Year\t" << "nCyl\t" << "nDel\t" << "maxRupt\t" << "DiamBase\t" \
       << "BreastDiam\t" << "Height\t" << "nBr\t" << std::endl;
  info << Year << "\t" << nCyl << "\t" << num_deleted(nCyl,cyls) << "\t" << \
    maxRupt << "\t" << (2*cyls[0].radius) << "\t" << breast_height_diam(cyls,nCyl,sta) \
       << "\t" << tree_height(cyls,nCyl,sta.y) << "\t" << br.size() << "\t" << std::endl;

  info.close();
  return 0;
}

int pipe_model_debug(std::vector<CylData> &cyls,int nCyl,std::vector<bool>& pipe_failed)
{
  //Debugging the pipe model conditions
  float sap_offsets = 0.0;
  float sap_base = 0.0;
  int nodes = 0;
  bool is_failed = false;
  for(int i=0; i<nCyl; i++){
    if(cyls[i].is_deleted)
      continue;
    if(cyls[i].extension > -1 || cyls[i].children.size()){
      sap_base = sap_area(i);
      nodes++;
      // std::cout << nodes << ") " << i << "==>";
      if(pipe_failed[i])
	is_failed = true;
      if(cyls[i].extension > -1){
	sap_offsets += sap_area(cyls[i].extension);
	// std::cout << cyls[i].extension;
      }
      // std::cout << "; ";
      for(int j = 0; j < cyls[i].children.size(); j++){
	sap_offsets += \
	  sap_area(cyls[i].children[j]) + (1-PIPECOEFF)*sap_area(cyls[i].children[j]);
	// std::cout << cyls[i].children[j] << ",";
      }
      if(is_failed){
	pipe_model_debug_print(i,pipe_failed[i]);
	// std::cout << "Base " << i << " pipe failed: " << \
	//   (sap_base-sap_offsets) << std::endl;
	// std::cout << "base SA = " << sap_base << " > " \
	//  	  << sap_offsets << " = offsets SA";
      }
      else{
	if(fabs(sap_base-sap_offsets) > 1e-07){
	  std::cout << std::endl << "WARNING: sapwoods do not match: " << std::endl;
	  std::cout << "\t base " << i << ", diff: " <<\
	    (sap_base-sap_offsets) << std::endl;
	  pipe_model_debug_print(i,pipe_failed[i]);
	}
	// std::cout << "base SA = " << sap_base << " =?= " \
	// 	  << sap_offsets << " = offsets SA";
      }
      // std::cout << std::endl;
    }
    is_failed = false;
    sap_offsets = 0.0;
  }
  return 0;
}

int pipe_model_debug_print(int i,bool is_failed)
{
  float sap_offsets = 0.0;
  float sap_base = 0.0;
  if(cyls[i].extension > -1 || cyls[i].children.size()){
    sap_base = sap_area(i);
    std::cout << "Base cyl " << i << "==>";
    if(cyls[i].extension > -1){
      sap_offsets += sap_area(cyls[i].extension);
      std::cout << cyls[i].extension;
    }
    std::cout << "; ";
    for(int j = 0; j < cyls[i].children.size(); j++){
      sap_offsets += \
	sap_area(cyls[i].children[j]) + (1-PIPECOEFF)*sap_area(cyls[i].children[j]);
      std::cout << cyls[i].children[j] << ",";
    }
    if(is_failed){
      std::cout << "Base " << i << " pipe failed: " << \
	(sap_base-sap_offsets) << std::endl;
      std::cout << "base SA = " << sap_base << " > " \
		<< sap_offsets << " = offsets SA";
    }
    else{
      std::cout << "base SA = " << sap_base << " =?= " \
		<< sap_offsets << " = offsets SA";
    }
    std::cout << std::endl;
  }
  return 0;
}

int pipe_model(std::vector<CylData>&cyls,int nCyl,int nShoot,float A[],\
	       std::vector<float>& sa_final,std::vector<float>& sa_current)
{
  //TERMINOLOGY:
  //SHOOTS are cyl's produced this year with not known length/radius (dependent on Lambda)
  //LAYER is a set of cyl's having the same number of cyl's on the path from the base
  //cyl (0) to them.
  //TIPS layer is a layer of cyls immediately beneath the shoots produced this year.
  //
  //The sapwood area (SA) increment induced by the shoots onto the cyl's beneath is
  //equal to:
  //******************************************************
  //*****      SA_delta = SA_final - SA_current     ******
  //******************************************************
  //SA_final is dependent on Lambda, since for the tips SA_final = \sum (SA_shoots),
  //where as SA_shoots is function of Lambda^2 (A[2]).
  //SA_current is the current SA for a cyl (including account on the senescence and
  //release due to the foliage shedding),i.e. a constant (A[0]). More generally, it
  //is "everything that needs not to be grown"
  //SA_delta is equal to ZERO for the dead-end cyl's, i.e. those not having any
  //shoots produced this year. And hence all cyl's beneath the dead-end cyl's also do
  //not receive any radial growth from them.
  //IMPORTANT: SA_delta is propagating down the layers, and NOT SA_final or
  //SA_current alone. We separate the two, only because of (non)dependency on Lambda.
  //ALGORITHM outline:
  //1. Define the new SA introduced by the shoots
  //2. Define the tip-layer cyl's that receive the newly allocated SA
  //2a. Define the tip-layer of dead-ends.
  //3. Process the tips-layer that receive the growth from the newly allocated SA
  //4. Process the layers down from the tip-layer to the base to determine how much
  //   of SA needs to be grown
  // All growth in terms of MASS.
  
  //SA_final: initialized to zero for all cyl's. NOTE: SA_final = SA/Lambda^2
  sa_final.assign(nCyl-nShoot,0.0);
  //SA_current
  sa_current.assign(nCyl-nShoot,0.0);
  //SA_fol: SA associated with foliage release of the offsets of a cyl
  std::vector<float> sa_fol(nCyl-nShoot,0.0);
  //SA_ss: SA associated with sapwood senescence of the offsets of a cyl. NOTE: this
  //is matched automatically given the same global SS rate. However, if the rate is
  //different for different cyl's this should be matched as well
  std::vector<float> sa_ss(nCyl-nShoot,0.0);
  //Layer members: logical array of cyl belonging to the shoot-pathway or
  //dead-end-pathway  down
  std::vector<bool> shoot_layer(nCyl-nShoot,false);
  std::vector<bool> dead_layer(nCyl-nShoot,false);
  //std::vector<int> shoot_postponed(0);//there must be no shoot_postponed
  std::vector<int> dead_postponed(0);
  std::vector<int> shoot_postponed(0);
  //Processed cyl's
  std::vector<bool> processed_cyl(nCyl-nShoot,false);

  //***********************************
  // 1. New SA introduced by the shoots
  for(int i = nCyl-nShoot; i < nCyl; i++){
    if(cyls[cyls[i].parent].extension == i){//Extension 1-to-1 correspondence
      sa_final[cyls[i].parent] += 1.0 * sap_area(i);
    }
    else{
      sa_final[cyls[i].parent] += PIPECOEFF * sap_area(i);
    }
    //**************************************************
    //2. Meanwhile form the tip layer of the shoot-ends.
    shoot_layer[cyls[i].parent] = true;
  }
  //************************************
  //2a. The tip-layer of the dead-ends
  for(int i = 0; i < nCyl-nShoot; i++){
    if(cyls[i].is_deleted)
      continue;
    if(cyls[i].extension == -1 && cyls[i].children.size() == 0){
      dead_layer[i] = true;
    }
    //We do not need to grow anything for the tip layer
  }

  //*********************************************************
  // 3. The tips-layer processing (needs separate processing)
  //std::cout << "=== Processing the tip layer: " << std::endl;
  std::vector<int> shoot_tmp (0);
  std::vector<int> dead_tmp (0);
  for(int i = 0; i < nCyl-nShoot; i++){
    if(shoot_layer[i]){
      if(!processed_cyl[i])
	processed_cyl[i] = true;
      else
	std::cout << "Warning: already processed cyl" << std::endl;
      //How much needs to be grown in terms of MASS?
      //+final SA = f(Lambda^2)
      A[2] += sa_final[i] * cyls[i].length * RHO;
      //-current SA = Const (+ real current area - senescenced + foliage released)
      if(cyls[i].age > SSSTART)
	sa_ss[i] = SS * sap_area(i);
      sa_fol[i] = sap_area_fol(i);
      sa_current[i] = sap_area(i) - sa_ss[i] + sa_fol[i];
      A[0] -= sa_current[i] * cyls[i].length * RHO;
      //Remove the processed cyl from the layer stack
      shoot_layer[i] = false;
      //Form the next layer
      if(cyls[i].parent > -1){
	if(cyls[cyls[i].parent].extension == i){
	  sa_final[cyls[i].parent] += 1.0 * sa_final[i];
	  sa_current[cyls[i].parent] += 1.0 * sa_current[i];
	  sa_ss[cyls[i].parent] += 1.0 * sa_ss[i];
	  sa_fol[cyls[i].parent] += 1.0 * sa_fol[i];
	}
	else{
	  sa_final[cyls[i].parent] += PIPECOEFF * sa_final[i];
	  sa_current[cyls[i].parent] += PIPECOEFF * sa_current[i];
	  sa_ss[cyls[i].parent] += PIPECOEFF * sa_ss[i];
	  sa_fol[cyls[i].parent] += PIPECOEFF * sa_fol[i];
	}
	shoot_tmp.push_back(cyls[i].parent);
      }
    }
    if(dead_layer[i]){
      if(!processed_cyl[i])
	processed_cyl[i] = true;
      else
	std::cout << "Warning: already processed cyl" << std::endl;
      if(cyls[i].age > SSSTART)
    	sa_ss[i] = SS * sap_area(i);
      sa_fol[i] = sap_area_fol(i);
      dead_layer[i] = false;
      //Form the next layer and make sure the parent is not from the shoot path
      if(cyls[i].parent > -1){
    	if(cyls[cyls[i].parent].extension == i){
    	  sa_ss[cyls[i].parent] += 1.0 * sa_ss[i];
    	  sa_fol[cyls[i].parent] += 1.0 * sa_fol[i];
    	}
    	else{
    	  sa_ss[cyls[i].parent] += PIPECOEFF * sa_ss[i];
    	  sa_fol[cyls[i].parent] += PIPECOEFF * sa_fol[i];
    	}
    	dead_tmp.push_back(cyls[i].parent);
      }
    }
  }
  //The next layer
  for(int i = 0; i < shoot_tmp.size(); i++){
    shoot_layer[shoot_tmp[i]] = true;
  }
  for(int i = 0; i < dead_tmp.size(); i++){
    dead_layer[dead_tmp[i]] = true;
  }
  shoot_tmp.resize(0);
  dead_tmp.resize(0);

  //************************************
  // 4. The layers beneath the tip-layer
  int layer = 0;
  while(!is_empty_array(shoot_layer) || !is_empty_array(dead_layer)){
    layer++;//layer is equal to age+1 of the cyl's of the layer.
    //memb_in_layer(shoot_layer);
    //Postpone some layers, if do not match by age/layer number.
    postpone_layer_members(layer,dead_layer,dead_postponed);
    postpone_layer_members(layer,shoot_layer,shoot_postponed);
    if(!is_same_age(shoot_layer)){
      std::cout << "WARNING: different age for shoot path" << std::endl;
    }
    if(!is_same_age(dead_layer)){
      std::cout << "WARNING: different age for dead path" << std::endl;
    }

    //std::cout << "=== layer " << layer << " from the tip" << std::endl;
    for(int i = 0; i < nCyl-nShoot; i++){
      if(cyls[i].is_deleted)
	continue;
      if(shoot_layer[i] || dead_layer[i]){
	if(!processed_cyl[i])
	  processed_cyl[i] = true;
	else
	  std::cout << "***** Warning: already processed cyl!" << std::endl;
	//Now, sa_final and sa_current are accumulated from the top layer
	//i.e. sa_final*Lambda^2 - sa_current is the SA_delta. NOTE: sa_final is zero
	// for the dead-ends.
	A[2] += sa_final[i] * cyls[i].length * RHO;
	//Need to grow difference (SA_fol - SA_ss) for the cyl's above
	sa_current[i] -= sa_fol[i] - sa_ss[i];
	//Modify SA_fol and SA_ss to those of the current cyl, pass down the layers
	if(cyls[i].age > SSSTART)
	  sa_ss[i] = SS * sap_area(i);
	sa_fol[i] = sap_area_fol(i);
	//No need to grow (SA_fol - SA_ss) of the current cyl
	sa_current[i] += sa_fol[i] - sa_ss[i];
	//Form the term
	A[0] -= sa_current[i] * cyls[i].length * RHO;

	//Form the next layer
	if(cyls[i].parent > -1){
	  if(cyls[cyls[i].parent].extension == i){
	    sa_final[cyls[i].parent] += 1.0 * sa_final[i];
	    sa_current[cyls[i].parent] += 1.0 * sa_current[i];
	    sa_fol[cyls[i].parent] += 1.0 * sa_fol[i];
	    sa_ss[cyls[i].parent] += 1.0 * sa_ss[i];
	  }
	  else{
	    sa_final[cyls[i].parent] += PIPECOEFF * sa_final[i];
	    sa_current[cyls[i].parent] += PIPECOEFF * sa_current[i];
	    sa_fol[cyls[i].parent] += PIPECOEFF * sa_fol[i];
	    sa_ss[cyls[i].parent] += PIPECOEFF * sa_ss[i];
	  }
	  //Remove the processed cyl from the layer stack
	  if(shoot_layer[i]){
	    shoot_tmp.push_back(cyls[i].parent);
	  }
	  else if(dead_layer[i] && !shoot_layer[i]){
	    dead_tmp.push_back(cyls[i].parent);
	  }
	  shoot_layer[i] = false;
	  dead_layer[i] = false;
	}
	else{//No parent exist, base cyl
	  shoot_layer[i] = false;
	  dead_layer[i] = false;
	}
      }
    }
    // if(is_empty_array(shoot_layer))
    //   std::cout << "SHOOT IS EMPTY\n";
    // if(is_empty_array(dead_layer))
    //   std::cout << "DEAD IS EMPTY\n";
    for(int i = 0; i < shoot_tmp.size(); i++){
      shoot_layer[shoot_tmp[i]] = true;
    }
    shoot_tmp.resize(0);
    for(int i = 0; i < dead_tmp.size(); i++){
      dead_layer[dead_tmp[i]] = true;
    }
    dead_tmp.resize(0);
    for(int i = 0; i < dead_postponed.size(); i++){
      dead_layer[dead_postponed[i]] = true;
    }
    dead_postponed.resize(0);
    for(int i = 0; i < shoot_postponed.size(); i++){
      shoot_layer[shoot_postponed[i]] = true;
    }
    shoot_postponed.resize(0);
  }

  return 0;
}


bool is_same_age(std::vector<bool>& cyl_idx)
{//Check whether cyls indicated by cyl_idx have the same year
  int age = -1;
  bool first_run = true;
  for(int i = 0; i < cyl_idx.size(); i++){
    if(cyl_idx[i] && first_run){
      age = cyls[i].age;
      first_run = false;
    }
    else if(cyl_idx[i] && !first_run){
      if(cyls[i].age == age)
	continue;
      else
	return false;
    }
  }
  return true;
}

int memb_in_layer(std::vector<bool>& idx)
{
  std::cout << "Layer members: " << std::endl;
  for(int i = 0; i < idx.size(); i++){
    if(idx[i]){
      std::cout << i << ", ";
    }
  }
  std::cout << std::endl;
  return 0;
}

int postpone_layer_members(int layer,\
			   std::vector<bool>& dead_layer,std::vector<int>& postponed)
{//Function for the PIPE_MODEL. Postpones for the next iteration the layer members,
 //that do not match the age.
  //1. Find the age as layer+1, since the tip-layer has age=1.
  //2. Return an array of members of the dead layer, that do not equal the age and make
  //corresponding dead_layer entries false.
  
  int age = layer+1;
  // for(int i = 0; i < shoot_layer.size(); i++){
  //   if(shoot_layer[i]){
  //     age = cyls[i].age;
  //     break;
  //   }
  // }
  // //std::cout << "Age " << age << std::endl;

  //std::cout << "Postponed " << std::endl;
  for(int i = 0; i < dead_layer.size(); i++){
    if(dead_layer[i]){
      if(cyls[i].age != age){
	//std::cout << i << std::endl;
	postponed.push_back(i);
	dead_layer[i] = false;
      }
    }
  }
  return 0;
}

int print_deleted_cyls(int nCyl)
{
  bool first_flag = true;
  for(int i = 0; i < nCyl; i++){
    if(cyls[i].is_deleted){
      if(first_flag){
	std::cout << "Deleted cyls:\n";
	first_flag = false;
      }
      std::cout << i << ", ";
    }
  }
  std::cout << std::endl;

  return 0;
}
