
int propagate_shadow(float I[], float ip[],std::vector<CylData> &cyls,\
		     const int nCyl, VoxelBox vb, int Years, 	      \
		     Array3D &shadow,int &fatal_error)
{//Calculate the shadow values for the voxel space and return the Light Conditions
 //coefficient IP for each cylinder.
  int i,j;
  //NOTE: 3rd dimension of shadow is UP direction, i.e. vector Y.
  V3f cyl_point;
  //We will not use the SORTING by the height of the cyl's (maybe needed?)
  //Calculate the absolute value of the radiation falling unto the cyl
  float lightExp[nCyl];
  float shadow_sum;
  V3f xz_projection;
  float cos_value;
  int xv,yv,zv;//voxel indices along X-, Y- and Z-span of the VoxelBox
  float axis_fraction = 1.0/((float)NPOINCYL+1.0);
  //Calculate shadow propagation from each of the cyl's
  for(i = 0; i < nCyl; i++){
    if(cyls[i].is_deleted)
      continue;
    cyl_point = cyls[i].start + 0.5*cyls[i].length*cyls[i].axis;
    //Define the voxels over X,Y and Z spans. 0-based numbering.
    xv = ceil( (cyl_point.x - vb.offset.x) / vb.vdx ) - 1;
    if(xv < 0 || xv >= shadow.m_width){
      std::cerr << "FATAL ERROR: X span: " << xv << std::endl;
      std::cout << "cyl " << i << ": " << cyls[i] << std::endl;
      fatal_error = 1;
      return 1;
    }
    yv = ceil( (cyl_point.y - vb.offset.y) / vb.vdy ) - 1;
    if(yv < 0 || yv >= shadow.m_height){//This is UP direction of the tree.
      std::cerr << "FATAL ERROR: Y span: " << yv << std::endl;
      std::cout << "cyl " << i << ": " << cyls[i] << std::endl;
      fatal_error = 1;
      return 1;
    }
    zv = ceil( (cyl_point.z - vb.offset.z) / vb.vdz ) - 1;
    if(zv < 0 || zv >= (shadow.m_data.size()/(shadow.m_width*shadow.m_height))){
      std::cerr << "FATAL ERROR: Z span: " << zv << std::endl;
      std::cout << "cyl " << i << ": " << cyls[i] << std::endl;
      fatal_error = 1;
      return 1;
    }
    //== PYRAMID ===
    int qmax = ceil((float)SHLEN / vb.vdy);
    //std::cout << "qmax = " << qmax << std::endl;
    for(int iy = (yv-1); iy >= (yv-qmax); iy--){//Iterate over the layers down, Y-span
      for(int ix = (xv-(yv-iy)); ix <= (xv+(yv-iy)); ix++){//Iterate over X-span
	for(int iz = (zv-(yv-iy)); iz <= (zv+(yv-iy)); iz++){//Iterate over Z-span
	  //If any of the indices is outside the range we skip
	  if( (ix < 0) || (ix >= GSX) || (iy < 0) || \
	      (iy >= GSY) || (iz < 0) || (iz >= GSZ) )
	    continue;
	  // Exponential
	  //shadow(ix,iz,iy) += SH * pow(EXPBASE,-(yv-iy));
	  // Quadratic
	  // shadow(ix,iz,iy) += (float)SH *				\
	  //   ( -(float)(yv-iy) * (float)(yv-iy) / ((float)qmax*(float)qmax) + 1.0 );
	  // Linear
	  shadow(ix,iz,iy) += (float)SH * ( -(float)(yv-iy)/((float)qmax) + 1.0 );
	}
      }
    }
  }
  //***********************************************************
  //RECIPIENT CYL's
  //Calculate the IP coefficients
  for(i = 0; i < nCyl; i++){
    shadow_sum = 0.0;
    if(cyls[i].is_deleted)
      continue;
    // Calculate the shadow for the cyl
    for(j = 0; j < NPOINCYL; j++){
      cyl_point = cyls[i].start + (j+1)*axis_fraction*cyls[i].length*cyls[i].axis;
      xv = ceil( (cyl_point.x - vb.offset.x) / vb.vdx ) - 1;
      yv = ceil( (cyl_point.y - vb.offset.y) / vb.vdy ) - 1;
      zv = ceil( (cyl_point.z - vb.offset.z) / vb.vdz ) - 1;
      shadow_sum += shadow(xv,zv,yv)/(float)NPOINCYL;//average shadow
      //std::cout << shadow(xv,zv,yv) << std::endl;
    }
    //std::cout << "Avg: " << shadow_sum << std::endl;
    //Calculate the light exposure of the cyl
    // xz_projection = V3f(cyls[i].axis.x,0.0,cyls[i].axis.z);
    // if( xz_projection.Length() < TOL )
    //   cos_value = 0.0;
    // else{
    //   xz_projection.Normalize();
    //   cos_value = xz_projection * cyls[i].axis;
    // }
    if((cyls[i].Wf > TOL) && (QRAD > shadow_sum)){
      // I[i] = ((float)QRAD - shadow_sum) *		  \
      // 	(2*(cyls[i].radius+RF)*cyls[i].length*cos_value + \
      // 	 M_PI*(cyls[i].radius+RF)*(cyls[i].radius+RF)*sqrt(fabs(1.0-cos_value*cos_value)));
      // ip[i] = I[i]/\
      // 	( (float)QRAD *							\
      // 	  (2*(cyls[i].radius+RF)*cyls[i].length*cos_value +		\
      // 	   M_PI*(cyls[i].radius+RF)*(cyls[i].radius+RF)*sqrt(1.0-cos_value*cos_value)) );
      I[i] = (float)QRAD - shadow_sum;
      ip[i] = I[i]/(float)QRAD;
    }
    else{
      I[i] = 0.0;
      ip[i] = 0.0;
    }
    //std::cout << I[i] << ", " << ip[i] << std::endl;
  }

  return 0;
}

int shadow_preset(Array3D & shadow, VoxelBox vb, V3f sta, std::vector<CylData> & cyls, int nCyl)
{
  shadow.Reset(0.0);//Reset to 0.0 the shadow voxel box

  //********************************************
  //Below is the any pattern for the shadow box.
  //********************************************

  float sh_dist = 4.0;//dist at which the environment shades
  int xv,yv,zv;
  V3f cyl_point;
  // Average shadow field exerted by the surrounding
  // for(int i=0; i<nCyl; i++){
  //   if(cyls[i].is_deleted || cyls[i].order == 0)
  //     continue;
  //   if( ( cyls[i].end - sta ).Length() > (sh_dist - (float)SHLEN ) ){
  //     //cyl_point = cyls[i].start + 0.5*cyls[i].length*cyls[i].axis;
  //     cyl_point = cyls[i].end;
  //     xv = ceil( (cyl_point.x - vb.offset.x) / vb.vdx ) - 1;
  //     if(xv < 0 || xv >= shadow.m_width)
  // 	continue;
  //     yv = ceil( (cyl_point.y - vb.offset.y) / vb.vdy ) - 1;
  //     if(yv < 0 || yv >= shadow.m_height)
  // 	continue;
  //     zv = ceil( (cyl_point.z - vb.offset.z) / vb.vdz ) - 1;
  //     if(zv < 0 || zv >= (shadow.m_data.size()/(shadow.m_width*shadow.m_height)))
  // 	continue;
  //     shadow(xv,zv,yv) = (float)SH;
  //   }
  // }

  
  //Bright cube
  // for(int i = 0; i < GSY; i++){
  //   for(int j = 0; j < GSX; j++){
  //     for(int k = 0; k < GSZ; k++){
  // 	if( k >= 60 && k <= 130 && j >= 60 && j <= 130)
  // 	  continue;
  // 	shadow(j,k,i) = 1000*QRAD;
  //     }
  //   }
  // }
  //External layer of the voxel box is dark
  // int thick = 150;
  // for(int i = 0; i < GSY; i++){
  //   for(int j = 0; j < GSX; j++){
  //     for(int k = 0; k < GSZ; k++){
  // 	if( j < thick || j >= (GSX-thick) )
  // 	  shadow(j,k,i) = QRAD;
  // 	else{
  // 	  if(k < thick || k >= (GSZ-thick))
  // 	    shadow(j,k,i) = QRAD;
  // 	}
  //     }
  //   }
  // }
  //Centered pyramidal penumbra
  // float highest = 0.0;
  // float distal_x = 0.0;
  // float distal_z = 0.0;
  // for(int i=0;i<nCyl;i++){
  //   if(cyls[i].is_deleted)
  //     continue;
  //   if(highest < cyls[i].end.y){
  //     highest = cyls[i].end.y;
  //   }
  //   if(distal_x < fabs(cyls[i].end.x))
  //     distal_x = fabs(cyls[i].end.x);
  //   if(distal_z < fabs(cyls[i].end.z))
  //     distal_z = fabs(cyls[i].end.z);
  // }
  // int yv = ceil( (highest - vb.offset.y) / vb.vdy ) - 1;
  // int dxv = ceil( (fabs(distal_x-sta.x)) / vb.vdx ) - 1;
  // int dzv = ceil( (fabs(distal_z-sta.z)) / vb.vdz ) - 1;
  // //(0,0,0) point base or other sta
  // int yb = ceil( (sta.y - vb.offset.y) / vb.vdy ) - 1;
  // int xb = ceil( (sta.x - vb.offset.x) / vb.vdx ) - 1;
  // int zb = ceil( (sta.z - vb.offset.z) / vb.vdz ) - 1;
  
  // //std::cout << yv << ", " << dxv << ", " << dzv << ", " << yb << std::endl;
  // // float pivot_x[] = {0.5,-0.1,0.1,-0.1};
  // // float pivot_z[] = {0.5,0.5,-0.5,-0.5};
  // // std::cout << "Iterate over Y: " << yv-1 << ":" << yb << std::endl;
  // // std::cout << "Iterate over X: " << xb-dxv << ":" << dxv+xb << std::endl;
  // // std::cout << "Iterate over Z: " << zb-dzv << ":" << dzv+zb << std::endl;
  // // std::cout << yv << ", " << yb << ", " << std::endl;
  // for(int iy = yv; iy >= yb; iy--){//Iterate over the layers down, Y-span
  //   for(int ix = xb-dxv; ix <= dxv+xb; ix++){//Iterate over X-span
  //     for(int iz = zb-dzv; iz <= dzv+zb; iz++){//Iterate over Z-span
  // 	if( (ix < 0) || (ix >= GSX) || (iy < 0) || \
  // 	    (iy >= GSY) || (iz < 0) || (iz >= GSZ) )
  // 	  continue;
  // 	//shadow(ix,iz,iy) += QRAD * pow(1.1,-(iy-yb));
  // 	shadow(ix,iz,iy) += QRAD *					\
  // 	  -(float)(iy-yb) * (float)(iy-yb) / ((float)(yv-yb)*(yv-yb)) + 1.0;
  //     }
  //   }
  // }
  return 0;
}
