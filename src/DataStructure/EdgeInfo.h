//
// {EINFOT, EINFOSINGLEM, EINFOV, EINFOBD, EINFOBIM, EINFOBDM} = Range[6];
// . for each triangle using the edge, store the triangle# and the index of this edge in the triangle datastructure
// . per-tri measure, length of edge
// . the two end vertex
// . is input edge
// . bi-tri measure 
// . bd-tri measure
// example: 
// {{{1, 1}, {2, 1}, {3, 1}, {4, 1}}, 
// 277.128, 
// {1, 2}, 
// 1, 
// {2.47009, 2.75679, 2.8269, 2.85489, 2.78478, 3.07148}, 
// {}}
// 
#ifndef _EDGE_INFO_H_
#define _EDGE_INFO_H_

class EdgeInfo {
public:
	EdgeInfo();
	~EdgeInfo();

	int v1;
	int v2;
	int nTri;		
	int* tris;
	char* edgeInd;
	float sgMeasure;
	float* biMeasure;
	float* bdMeasure;
	bool isBD;
	
	int getSize();
private:
};

#endif
