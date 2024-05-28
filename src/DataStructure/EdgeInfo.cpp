#include "EdgeInfo.h"

EdgeInfo::EdgeInfo(){
	v1 = v2 = nTri = -1;	
}

EdgeInfo::~EdgeInfo(){
	delete [] tris;  delete [] edgeInd; 
	delete [] biMeasure; delete [] bdMeasure;
}

int EdgeInfo::getSize(){
	return 0;
}
