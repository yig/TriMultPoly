#ifndef _DMWT_H_
#define _DMWT_H_
#include "tetgen.h"
#include "Point3.h"
#include "Vector3.h"
#include "CurveInfo.h"
#include "EdgeInfo.h"
#include "TriangleInfo.h"
#include "Configure.h"
#include "Boundary.h"
#include "Hole.h"
#include "Tile.h"
#include "SubKey.h"
#include <boost/unordered_map.hpp>
#include <boost/heap/priority_queue.hpp>
#include <utility> 
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <set>
#include <cmath>
#include <stdio.h>
#include <time.h>

#define point(v) Point3((float)vertices[v*3],(float)vertices[v*3+1],(float)vertices[v*3+2])
#define triangle(t,i) tris[t*3+i]
#define vertex(v,i) vertices[v*3+i]
#define TinE(ei,ti) edgeInfoList[ei]->tris[ti]
#define EIinE(ei,ti) edgeInfoList[ei]->edgeInd[ti]
#define EinT(ti,ei) triangleInfoList[ti]->edges[ei]
#define TIinT(ti,ei) triangleInfoList[ti]->triInd[ei]
#define isBoardE(e) edgeInfoList[e]->isBD
#define PI 3.1415926f

typedef boost::unordered_map< __int64, size_t> wekhashmap;

typedef boost::unordered_map< SubKey<MAXK>, int, MySubKeyHash<MAXK>,  MySubKeyEq<MAXK>> hashmap;

class DMWT{
public:
	DMWT(){}
	// use candidate triangle file
	DMWT(char* curvefile, char* canTfile, const bool pDP, const bool uWE, const bool fWE, const bool sObj){
		useNormal = false;
		badInput=false;
		isDeGen=false;

		readCurveFile(curvefile);
		readCanTFile(canTfile);
		useWE = uWE;
		paraDP    = pDP;
		filterWE = fWE;
		useMinMaxDihedral = false;
		saveObj = sObj;

		numoftilingtris = 0;
	}
	// without normal
	DMWT(char* curvefile, const bool uDT, const bool pDP, const bool uWE, const bool fWE, const bool sObj){
		useNormal = false;
		badInput=false;
		isDeGen=false;

		readCurveFile(curvefile);
		if (uDT){
			genDTTriCandidates();
		}else{
			genAllTriCandidates();
		}
		useWE = uWE;
		useDT = uDT;
		paraDP    = pDP;
		filterWE = fWE;
		useMinMaxDihedral = false;
		saveObj = sObj;

		numoftilingtris = 0;
	}
	// with normal
	DMWT(char* curvefile, char* normalfile, const bool uDT, const bool pDP, const bool uWE, const bool fWE, const bool sObj){
		useNormal = true;
		badInput=false;
		isDeGen=false;

		readCurveFile(curvefile);
		readNormalFile(normalfile);
		if (uDT){
			genDTTriCandidates();
		}else{
			genAllTriCandidates();
		}
		useWE = uWE;
		useDT = uDT;
		paraDP    = pDP;
		filterWE = fWE;
		useMinMaxDihedral = false;
		saveObj = sObj;

		numoftilingtris = 0;

		//cout<<"points------------"<<endl;
		//for (int i=0; i<numofvertices; i++){
		//	point(i).print();
		//}
		//cout<<"normals------------"<<endl;
		//for (int ni=0; ni<numofvertices; ni++){
		//	Vector3(normals[ni*3],normals[ni*3+1],normals[ni*3+2]).print();
		//}
	}
	~DMWT();
	

	float optCost;
	std::vector<int> optTile;

	//-------------evaluations--------------//
	long numSub;
	float timeReadIn;
	float timePreprocess;
	float timeTile;
	float timeTotal;
	float timeTetgen;


//protected:
	//-------------variables--------------//
	char* filename;
	int numofcurves;
	int numofnormals;
	int numofvertices;
	int numoftris;
	int numofedges;
	int numoftilingtris;
	int orgnumofvertices;

	int* tris;
	double* vertices;
	float* normals;
	int* vertex2curve;

	CurveInfo *  curveInfo;
	EdgeInfo ** edgeInfoList;
	TriangleInfo ** triangleInfoList;
	hashmap TriangulationHash;

	int* tiling;
	int zeroTile;
	int infinityTile;
	vector<OneTile> allOneTiles;
	vector<Tile> allTiles;
	Tile resTile;

	tetgenio tetgenIn;
	tetgenio tetgenOut;
	int startE;
	bool useBiTri;
	bool useWE;
	bool useDT;
	bool paraDP;
	bool filterWE;
	bool saveObj;
	bool useNormal;

	float weightTri;
	float weightEdge;
	float weightBiTri;
	float weightTriBd;
	bool useMinMaxDihedral;

	int maxTriPerEdge;

	int ** ehash;
	int ** ehashSize;

	int nPt;

	std::vector<int> mBEs1;
	std::vector<int> mBEs2;
	std::vector<int> mBBEEs;
	std::vector<int> mEdges;

	std::vector<int> toAvoidWE1, toAvoidWE2, toAvoidWE;
	std::vector<int> newList, newList2;

	std::vector<bool> isSuper;

	//-------------functions--------------//
	void buildList();
	void buildStrip();
	void tile();
	void setWeights(float wtri, float wedge, float wbitri, float wtribd);
	void saveTiling();
	void saveTiling(char* tilefile);
	void saveTilingObj(char* tilefile);
	void saveTiling(int OTid);
	void statistics();

	int scanTrianglesOnce();
	void readCurveFile(const char* file);
	void readNormalFile(const char* file);
	void readCanTFile(const char* file);
	void genDTTriCandidates();
	void genAllTriCandidates();

	float measureEdge(int v1, int v2);
	float measureTri(int v1, int v2, int v3);
	float measureBiTri(int v1, int v2, int p, int q);
	float measureTriBd(int v1, int v2, int v3, int ni);
	float costTri(float measure);
	float costEdge(float measure);
	float costBiTri(float measure);
	float costTriBd(float measure);
	float CostE(int e);
	float CostAllE(int t);
	float CostT(int t);
	float CostBI(int e, short fi1, short fi2);
	float CostBD(int e, int i);

	int computeTriangulation(Boundary<MAXK>& B,Hole& H);
	void retrieveTile(Boundary<MAXK>& B, Hole& H, int OTid);
	void retrieveTile(int OTid);

	int whichCase(Boundary<MAXK>& B, Hole& H, int vv, int& HBindex);
	void DivideBoundary(Boundary<MAXK>& B, int ee1, int ee2, int v1, int v2, int vv, int ffi1, int ffi2, int intvInd,vector<Boundary<MAXK>>& boundary );
	void CombineBoundary(Boundary<MAXK>& B, int ee1, int ee2, int v1, int v2, int vv, int ffi1, int ffi2,vector<Boundary<MAXK>>& boundary);
	bool isInvalidB(Boundary<MAXK>& B);
	void GetAllWeakEdge(Boundary<MAXK>& B,vector<int>& weSet);
	bool isBoardEdge(const Boundary<MAXK>& B);
	bool isMeetingSpanE(const Boundary<MAXK>& B);
	void getIntervalsB(const Boundary<MAXK>& B, Interval<MAXK>& Intvs);
	__int64 GetBitKey(const vector<int>& canItem,const vector<int>& allItem);
	void RecoverListFromBitKey(const vector<int>& list, __int64 key,vector<int>& outlist );

	//-------------Edge Protection
	std::vector<Point3> tempC;
	std::vector<Point3> tempOrgC;
	std::vector<std::vector<int>> tempAdj;
	std::vector<Vector3> tempNorm;
	std::vector<std::vector<int>> tempAdjNorm;
	std::vector<int> newEdge;
	std::vector<int> newAdj;
	std::vector<int> newNorm;
	std::vector<char> newClip;
	std::vector<std::vector<int>> badEdge;
	std::vector<double> radius;
	std::vector<double> orgradius;
	std::vector<std::vector<char>> cliped;
	bool badInput;
	bool isDeGen;
	double* pb_vertices;

	bool EdgeProtection();
	bool passTetGen();
	void perturbPts(double ptb);
	bool isProtected();
	void protectCorner();
	double getAngle(int p1, int p2, int p3);
	double getPt2LineDist(int p1, int p2, int p3);
	void splitEdge(int p1, int p2ind, const Point3 & newP);
	void splitEdge(int p1, int p2ind, const Point3 &newP, const Point3 &newOrgP);
	void insertMidPointsTetgen();
	bool getCurveAfterEP();
	bool sameOrientation(const vector<int> & newCurve, int start, int len);

	//-------------Tiles
	int newTile(__int64 wec, float c, int t, int nxt1, int nxt2);
	int appendTile(__int64 wec, float c, int t, int nxt1, int nxt2, vector<int> &curTile);
	void appendTile(__int64 wec, float c, int t, int nxt1, int nxt2, vector<OneTile> &locTileVec);
	bool appendTile(vector<OneTile> &locTileVec, int & start, int & end);
	int findMinTile(int tileind);

	//-------------help funcs--------------//
	int getBiPos(int len, int i, int j);
	bool isBDEdge(int v1, int v2);

	void getSortV(int e, int s, int& v1, int& v2);	
	int nextV(int v, int dir);
	void getSubActiveEdge(int e, int ff, int v1, int v2, int& ee1, int& ee2, short&ffi1, short&ffi2);
	bool isWithinInterval(int pv, int nv, int v, int d);

	float findMax(float val1, float val2, float val3);
	float findMax(float val1, float val2, float val3, float val4);



};

#endif
