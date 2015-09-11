#ifndef FSTREAM_H
#include <fstream>
#define FSTREAM_H 1
#endif

void mtg_write(std::vector<CylData> &cyls, int nCyl)
{//Write the Multi-Scale Tree Graph (MTG) file output.
  std::ofstream mtg;
  std::vector<int> base_cyls(1,0);
  int cc;
  mtg.open("lignum.mtg",std::ofstream::out);
  while(base_cyls.size() != 0){
    // std::cout << "Base cyls: ";
    // for(int i=0;i<base_cyls.size();i++)
    //   std::cout << base_cyls[i] << " ";
    // std::cout << std::endl;
    cc = base_cyls[0];
    // if(cyls[cc].is_deleted){
    //   base_cyls.erase(base_cyls.begin());
    //   continue;
    // }
    if(cc != 0)
      mtg << "C" << cyls[cc].parent << "+";
    else
      mtg << "/";
    mtg << "C" << cc << " " << cyls[cc];

    //Go through the extensions
    while(cyls[cc].extension > -1){
      //Write all extension of the cyl CC
      mtg << "<C"<< cyls[cc].extension  << " " << cyls[cyls[cc].extension];
      //Keep the base cyls for other branches
      for(int i=0; i < cyls[cc].children.size(); i++){
      	base_cyls.push_back(cyls[cc].children[i]);
      }
      cc = cyls[cc].extension;
    }
    //Get children from the last iteration
    for(int i=0; i < cyls[cc].children.size(); i++){
      base_cyls.push_back(cyls[cc].children[i]);
    }

    //First cyl was processed, then remove
    base_cyls.erase(base_cyls.begin());
  }
  mtg.close();
}

