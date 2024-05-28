// 
// {TINFOE, TINFOSINGLEM, TINFOBDN} = Range[3];
// . for each of the 3 edge, store the index and the position of this triangle in the edge datastructure
// . single measure, the area of the triangle
// . the number of bd edges, can be {0,1,2,3}
// example: 
// {{{1, 1}, {2, 1}, {3, 1}}, 
// 41569.2, 
// 3}
// 
#ifndef _TRIANGLE_INFO_H_
#define _TRIANGLE_INFO_H_

class TriangleInfo {
public:
	TriangleInfo();
	~TriangleInfo();

	int edges [3];
	short triInd [3];	
	float sgMeasure;
	char nBDedge;
	
	int getSize();
private:
};

#endif
