#ifndef PICKER_H
#define PICKER_H

#include <vector>
#include <map>
#include <memory>
#include <stdexcept>
#if __GNUG__
#   include <tr1/memory>
#endif

#include "cvec.h"
#include "geometrymaker.h" // has typdef Matrix4 SgRbtNode


#include "glsupport.h"
#include "uniforms.h"
#include "shadergeometry.h"


#include "ppm.h"
//#include "drawer.h"


extern vector < shared_ptr<Object> >  g_objectList;



class Picker  {
  
  //std::vector<std::tr1::shared_ptr<SgNode> > nodeStack_;

 typedef std::map<int, std::tr1::shared_ptr<Object>  > IdToRbtNodeMap;
 
 IdToRbtNodeMap idToRbtNode_;

  int idCounter_;

  bool srgbFrameBuffer_;

 // Drawer drawer_;

  
public:
  Picker();
  
  void addToMap(int id, std::tr1::shared_ptr<Object> node );

  std::tr1::shared_ptr<Object> find(int id);

  Cvec4 idToColor(unsigned int id); // from id to color

  unsigned int colorToId( const unsigned char p[]);
  
  std::tr1::shared_ptr<Object> getRbtNodeAtXY(int x, int y);

  int getidCounter();

};


#endif
