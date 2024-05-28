#include "DMWT.h"
extern float timeNew;

//==================================Tile================================//

void DMWT::tile(){
	if (badInput) return;
	if (numofvertices == 3){
		optCost = CostAllE(0) + CostT(0);
		optTile.push_back(0);
		saveTiling();
		return;
	}
	numSub = 0;

	// set initial B and H
	Hole H;
	for (int i=1; i<numofcurves; i++){
		H.append(i);
	}
	Boundary<MAXK> B;
	B.append(startE,1,-1,1);
	TriangulationHash.reserve(numResvHash);
	allOneTiles.reserve(numResvOneTile);
	allTiles.reserve(numResvTile);
	mBEs1.reserve(numofcurves);
	mBEs2.reserve(numofcurves);
	mBBEEs.reserve(numofcurves);
	mEdges.reserve(numofcurves);

	//push a zero onetile and infinity onetile
	allOneTiles.push_back(OneTile(0, 0, -1, -1, -1));
	allOneTiles.push_back(OneTile(0, FLT_MAX, -1, -1, -1));

	//push a fake tile, zero tile and infinity tile
	Tile tmptile;
	allTiles.push_back(tmptile);
	tmptile.start = 0; tmptile.end = 0;
	allTiles.push_back(tmptile);
	tmptile.start = 1; tmptile.end = 1;
	allTiles.push_back(tmptile);

	zeroTile = 1;
	infinityTile = 2;

	int tile;
	tile = computeTriangulation(B,H);
	

	int minTile = findMinTile(tile);

	if (minTile == -1) return;
	optCost = allOneTiles[minTile].cost;

#if (_SaveTile==1)
	optTile.reserve(numResvTilling);
	retrieveTile(minTile);
	saveTiling();
#endif

}

void DMWT::retrieveTile(int OTid){
#if (_SaveTile==1)
	if (OTid == -1) return;
	if (allOneTiles[OTid].tri!=-1){
		optTile.push_back(allOneTiles[OTid].tri);
	}
	retrieveTile(allOneTiles[OTid].next1);
	retrieveTile(allOneTiles[OTid].next2);
#endif
}

int DMWT::findMinTile(int tileind){
	if (tileind==-1) return -1;

	int oti, minOti = -1;
	float minCost = FLT_MAX, curCost;
	int start=allTiles[tileind].start, end = allTiles[tileind].end;
	for (int i=0; i<end-start+1; i++){
		oti = start+i;
		curCost = allOneTiles[oti].cost;
		if (curCost<minCost){
			minCost = curCost;
			minOti = oti;
		}
	}
	return minOti;
}

/// Return the index of Tile(vector<int>) in allTiles
int DMWT::newTile(__int64 wec, float c, int t, int nxt1, int nxt2){
	allOneTiles.push_back(OneTile(wec, c, t, nxt1, nxt2));
	size_t nIndex(allTiles.size());
	allTiles.resize(nIndex+1);
	allTiles.back().start = (int)allOneTiles.size()-1;
	allTiles.back().end = (int)allOneTiles.size()-1;
	return (int)nIndex;
}

__forceinline void DMWT::appendTile(__int64 wec, float c, int t, int nxt1, int nxt2, vector<OneTile> &locTileVec){
	locTileVec.push_back(OneTile(wec, c, t, nxt1, nxt2));
}

/// Return true if OneTiles>0; else return false;
bool DMWT::appendTile(vector<OneTile> &locTileVec, int & st, int & en){
	if (locTileVec.size()<=1) 
		return false;
	int start = (int)allOneTiles.size(), end;
	if (filterWE){
		size_t locS = locTileVec.size();
		isSuper.clear();
		for (size_t i=0; i<locS; i++){
			isSuper.push_back(false);
		}
		__int64 andK;
		for (size_t i=1; i<locS; i++){
			if (isSuper[i] == true) continue;
			for (size_t j=i+1; j<locS; j++){
				andK = locTileVec[i].WECode & locTileVec[j].WECode;
				if (andK == locTileVec[j].WECode && (locTileVec[i].cost >= locTileVec[j].cost )){
					isSuper[i] = true; break;
				} else if (andK == locTileVec[i].WECode && (locTileVec[i].cost <= locTileVec[j].cost)){
					isSuper[j] = true; break;
				}
			}
		}
		for (size_t i=1; i<locS; i++){
			if (isSuper[i]) continue;
			allOneTiles.push_back(locTileVec[i]);
		}
	}else{
		allOneTiles.insert(allOneTiles.end(), locTileVec.begin()+1, locTileVec.end());
	}
	end = (int)allOneTiles.size()-1;

	st = start;
	en = end;
	if (en < st)
		return false;
	return true;
}

// return the index of OneTile in allOneTiles
__forceinline int DMWT::appendTile(__int64 wec, float c, int t, int nxt1, int nxt2, vector<int> &curTile){
	OneTile onetile(wec, c, t, nxt1, nxt2);
	allOneTiles.push_back(onetile);
	curTile.push_back((int)allOneTiles.size()-1);
	return (int)allOneTiles.size()-1;
}

/// Core dynamic algorithm
int DMWT::computeTriangulation(Boundary<MAXK>& B, Hole& H){
	SubKey<MAXK> subKey;
	B.getSubKey(subKey);
	H.getSubKey(subKey);

	int res;
	// exist in hash, return
	if((res=TriangulationHash[subKey])!=0) return res;

	// -------------------special cases--------------------
	if (B.isEmpty() && H.isEmpty()) return zeroTile;
	if (B.isEmpty() && !H.isEmpty()) return infinityTile;

	if (isBoardEdge(B)){
		if (H.isEmpty()){
			int e = B.getE(0);
			short fi = B.getF(0);
			//res = zeroTile;
			res = newTile(0, CostBD(e, fi), -1, -1, -1);
		} else {
			res = infinityTile; 
		}
		TriangulationHash[subKey] = res;
		numSub++;
		return res;
	}

	if (isMeetingSpanE(B)) {
		if (H.isEmpty()){
			int e = B.getE(0);
			short fi = B.getF(0);
			short ffi = B.getF(1);

			res     = newTile(0, CostBI(e, fi, ffi), -1, -1, -1);
		} else {
			res = infinityTile; 
		}
		TriangulationHash[subKey] = res;
		numSub++;
		return res;
	}

	// -------------------normal cases--------------------
	vector<int> allWE;
	wekhashmap WEKmap;
	if (useWE) {
		 GetAllWeakEdge(B,allWE);
	}
	vector<OneTile> locOneTiles;
	locOneTiles.push_back(OneTile());
	res = infinityTile;

	int e, v1, v2, nTri, HBindex;
	int ff, vv, ee1, ee2;
	short fi, ffi, ffi1, ffi2;
	char s, d, Case;
	int numSubBBs, numSubHHs;
	int oti1, oti2, oti;
	OneTile sTile;
	float costE;
	int * FF;
	int subRes1Start, subRes2Start, subResStart, subRes1End, subRes2End, subResEnd;

	B.getOneB(0,e,s,fi,d);
	nTri = edgeInfoList[e]->nTri;
	FF = edgeInfoList[e]->tris;
	getSortV(e,s,v1,v2);

	vector<Boundary<MAXK>> BBs;
	vector<pair<Hole,Hole>> HHs;
	vector<int> allWE1, allWE2, allWWEE;
	__int64 toAvoidK1,toAvoidK2, toAvoidK;

	for (int i=0; i<nTri; i++){
		ff  = FF[i];
		ffi = i;

		vv  = triangle(ff,0) + triangle(ff,1) + triangle(ff,2) - v1 - v2;
		getSubActiveEdge(e, ff, v1, v2, ee1, ee2, ffi1, ffi2);
		costE = CostAllE(ff);
		Case = whichCase(B, H, vv, HBindex);
		
		BBs.clear();
		HHs.clear();
		Hole HH;
		Hole H1;
		Hole H2;
		Boundary<MAXK> BB, B1, B2;
		
		toAvoidWE1.clear();
		toAvoidWE2.clear();
		toAvoidWE.clear();
		allWE1.clear();
		allWE2.clear();
		allWWEE.clear();

		switch(Case){
			case 1:
				DivideBoundary(B, ee1, ee2, v1, v2, vv, ffi1, ffi2, HBindex,BBs);
				HHs = H.partition();

				numSubBBs = BBs.size() == 4 ? 2 : 1;
				numSubHHs = (int)HHs.size();
				for (int b=0; b<numSubBBs; b++){
					B2 = BBs.back(); BBs.pop_back();
					B1 = BBs.back(); BBs.pop_back();
					B1.sort();
					B2.sort();
					if (isInvalidB(B1) || isInvalidB(B2)) {
						continue;
					}
					
					if (useWE){
						GetAllWeakEdge(B1,allWE1);
						GetAllWeakEdge(B2,allWE2);

						B1.getAllE(mBEs1);
						// boundary edges of B1 or B2 include e, nonmanifold, continue
						int be1, be2;
						for (be1=0; be1<B1.getSize(); be1++){
							if (mBEs1[be1] == e) break;
						}
						if (be1!=B1.getSize()){
							continue;
						}

						B2.getAllE(mBEs2);
						for (be2=0; be2<B2.getSize(); be2++){
							if (mBEs2[be2] == e) break;
						}
						if (be2!=B2.getSize()){
							continue;
						}
						// otherwise find all possible conflict edge between allWE1 and {allWE2, BEs2, e}, vice versa
						toAvoidWE1 = allWE2;
						toAvoidWE2 = allWE1;
						for (int be1=0; be1<B1.getSize(); be1++){
							toAvoidWE2.push_back(mBEs1[be1]);
						}
						toAvoidWE2.push_back(e);
						for (int be2=0; be2<B2.getSize(); be2++){
							toAvoidWE1.push_back(mBEs2[be2]);
						}

						toAvoidWE1.push_back(e);
						toAvoidK1 = GetBitKey(allWE1, toAvoidWE1);
						toAvoidK2 = GetBitKey(allWE2, toAvoidWE2);
					}
					for (int h=0; h<numSubHHs; h++){
						H1 = HHs[h].first;
						H2 = HHs[h].second;

						int nResultB1H1 = computeTriangulation(B1, H1);
						subRes1Start = allTiles[nResultB1H1].start;
						subRes1End = allTiles[nResultB1H1].end;

						int nResultB2H2 = computeTriangulation(B2, H2);
						subRes2Start = allTiles[nResultB2H2].start;
						subRes2End = allTiles[nResultB2H2].end;

						for (int tl1=0; tl1<subRes1End-subRes1Start+1; tl1++){
							for (int tl2=0; tl2<subRes2End-subRes2Start+1; tl2++){
								oti1=subRes1Start+tl1;
								oti2=subRes2Start+tl2;

								__int64 WEK1 = allOneTiles[oti1].WECode;
								float C1   = allOneTiles[oti1].cost;
								__int64 WEK2 = allOneTiles[oti2].WECode;
								float C2   = allOneTiles[oti2].cost;

								float newCost;
								__int64 newWEK = 0;
								bool goon = true;
								if (useWE) {
									goon = ((WEK1 & toAvoidK1) == 0) && ((WEK2 & toAvoidK2) == 0);
								}
								if (goon) {
									if (useMinMaxDihedral){
										newCost = findMax(CostBI(e, fi, ffi), C1, C2);
									} else {
										newCost = CostT(ff) + costE + CostBI(e, fi, ffi) + C1 + C2;
									}
									if (useWE){
										RecoverListFromBitKey(allWE1, WEK1,newList);
										RecoverListFromBitKey(allWE2, WEK2,newList2);
										newList.insert(newList.end(),newList2.begin(), newList2.end());
										newList.push_back(ee1);
										newList.push_back(ee2);
										newWEK = GetBitKey(allWE, newList);
									}
									size_t sTileInd = WEKmap[newWEK];
									if (sTileInd!=0){
										if(locOneTiles[sTileInd].cost > newCost){
											locOneTiles[sTileInd].cost = newCost;
											locOneTiles[sTileInd].tri = ff;
#if (_SaveTile==1)
									locOneTiles[sTileInd].next1 = oti1;
									locOneTiles[sTileInd].next2 = oti2;
#endif
										}
									} else {
										appendTile(newWEK, newCost,ff, oti1, oti2, locOneTiles);
										WEKmap[newWEK] = locOneTiles.size()-1;
									}
								}
							}
						}
					}
				}
				break;
			case 2:
				CombineBoundary(B, ee1, ee2, v1, v2, vv, ffi1, ffi2,BBs);
				HH = H;
				HH.erase(HBindex);
				numSubBBs = 2;

				for (int b=0; b<numSubBBs; b++){
					BB = BBs.back(); BBs.pop_back();
					BB.sort();
					if (useWE){
						GetAllWeakEdge(BB,allWWEE);

						BB.getAllE(mBBEEs);
						// boundary edges of B1 or B2 include e, nonmanifold, continue
						int be;
						for (be=0; be<BB.getSize(); be++){
							if (mBBEEs[be] == e) break;
						}
						if (be!=BB.getSize()){
							continue;
						}
						// otherwise find all possible conflict edge between allWE1 and {allWE2, BEs2, e}, vice versa
						toAvoidWE;
						toAvoidWE.push_back(e);
						toAvoidK = GetBitKey(allWWEE, toAvoidWE);
					}

					int nResultBBHH = computeTriangulation(BB, HH);
					subResStart = allTiles[nResultBBHH].start;
					subResEnd = allTiles[nResultBBHH].end;

					for (int tl=0; tl<subResEnd-subResStart+1; tl++){
						oti=subResStart+tl;

						__int64 WWEEK = allOneTiles[oti].WECode;
						float CC   = allOneTiles[oti].cost;
						float newCost;
						__int64 newWEK = 0;
						bool goon = true;
						if (useWE) {
							goon = ((WWEEK & toAvoidK) == 0);
						}
						if (goon) {
							if (useMinMaxDihedral){
								newCost = CostBI(e, fi, ffi);
								if (CC > newCost){
									newCost = CC;
								}
							} else {
								newCost = CostT(ff) + costE + CostBI(e, fi, ffi) + CC;
							}
							
							if (useWE){
								RecoverListFromBitKey(allWWEE, WWEEK,newList);
								newList.push_back(ee1);
								newList.push_back(ee2);
								newWEK = GetBitKey(allWE, newList);
							}
							size_t sTileInd = WEKmap[newWEK];

							if (sTileInd!=0){
								if(locOneTiles[sTileInd].cost > newCost){
									locOneTiles[sTileInd].cost = newCost;
									locOneTiles[sTileInd].tri = ff;
#if (_SaveTile==1)
									locOneTiles[sTileInd].next1 = oti;
									locOneTiles[sTileInd].next2 = -1;
#endif
								}
							} else {
								appendTile(newWEK, newCost,ff, oti, -1, locOneTiles);
								WEKmap[newWEK] = locOneTiles.size()-1;
							}
						}
					}
				}
				break;
			default:
				break;
		}
	}

	if (appendTile(locOneTiles, resTile.start, resTile.end)) {
		allTiles.push_back(resTile);
		res = (int)allTiles.size() - 1;
	}
	TriangulationHash[subKey] = res;
	numSub++;
	return res;
}


int DMWT::whichCase(Boundary<MAXK> &B, Hole& H, int vv, int& HBindex){
	int ci = vertex2curve[vv], c;
	// -------------Case II: lie on H----------------
	int Hind = H.find(ci);
	if (Hind!=-1){
		HBindex = Hind;
		return 2;
	}
	// -------------Case I: lie on B-----------------
	Interval<MAXK> intvs;
	getIntervalsB(B, intvs);
	int pv, nv;
	char d;
	for (int i=0; i<B.getSize(); i++){
		c = vertex2curve[intvs.getPV(i)];
		if (c != ci) 
			continue;
		intvs.getOneIntv(i, pv, nv, d);
		if (isWithinInterval(pv, nv, vv, d)){
			HBindex = i;
			return 1;
		} 
	}
	// -------Invalid Case: neither on B or H--------
	HBindex = -1;
	return -1;
}


void DMWT::DivideBoundary(Boundary<MAXK>& B, int ee1, int ee2, int v1, int v2, int vv, int ffi1, int ffi2, int intvInd,vector<Boundary<MAXK>>& newBs ){
	if (B.isEmpty()) {
		newBs.push_back(Boundary<MAXK>());
		newBs.push_back(Boundary<MAXK>());
		return;
	}
	// get new edge sign 
	int ss1, ss2;
	ss1 = vv == edgeInfoList[ee1]->v2 ? 1 : -1;
	ss2 = vv == edgeInfoList[ee2]->v1 ? 1 : -1;
	// get new direction
	char d1, di, dd1, ddi, si, sj;
	int ei, ej, vi1, vi2, vj1, vj2;
	int sizeB = B.getSize();
	d1 = B.getD(0);
	di = B.getD(intvInd);
	ei = B.getE(intvInd);
	si = B.getS(intvInd);
	ej = B.getE((intvInd+1)%sizeB);
	sj = B.getS((intvInd+1)%sizeB);

	getSortV(ei, si, vi1, vi2);
	getSortV(ej, sj, vj1, vj2);
	dd1 = (vv == vj1) ? 0 : di;
	ddi = (vv == vi2) ? 0 : di;
	// set up new boundaries
	if (vi2==vj1 && vv==vj1){
		Boundary<MAXK> B1, B2, B3, B4;
		B1.append(ee1,ss1,ffi1,di);
		B1.append(B, intvInd+1, sizeB-intvInd-1);
		B2.append(ee2, ss2, ffi2, d1);
		B2.append(B, 1, intvInd);
		B2.setD(intvInd, 0);

		B3.append(ee1, ss1, ffi1, 0);
		B3.append(B, intvInd+1, sizeB-intvInd-1);
		B4.append(ee2, ss2, ffi2, d1);
		B4.append(B, 1, intvInd);

		newBs.push_back(B1);
		newBs.push_back(B2);
		newBs.push_back(B3);
		newBs.push_back(B4);
	} else {
		Boundary<MAXK> B1,B2;
		B1.append(ee1, ss1, ffi1, dd1);
		B1.append(B, intvInd+1, sizeB-intvInd-1);
		B2.append(ee2, ss2, ffi2, d1);
		B2.append(B, 1, intvInd);
		B2.setD(intvInd, ddi);

		newBs.push_back(B1);
		newBs.push_back(B2);
	}
	
	return;
}


void DMWT::CombineBoundary(Boundary<MAXK> &B, int ee1, int ee2, int v1, int v2, int vv, int ffi1, int ffi2,vector<Boundary<MAXK>>& newBs){
	newBs.resize(2);
	if (B.isEmpty()) {	
		return;
	}
	int sizeB = B.getSize();
	// get new edge sign 
	int ss1, ss2;
	ss1 = vv == edgeInfoList[ee1]->v2 ? 1 : -1;
	ss2 = vv == edgeInfoList[ee2]->v1 ? 1 : -1;
	// get new direction
	char d1 = B.getD(0);
	// set up new boundaries 
	Boundary<MAXK> B1,B2;
	B1.append(ee1, ss1, ffi1, 1);
	B1.append(ee2, ss2, ffi2, d1);
	B1.append(B, 1, sizeB-1);
	B2.append(ee1, ss1, ffi1, -1);
	B2.append(ee2, ss2, ffi2, d1);
	B2.append(B, 1, sizeB-1);
	
	newBs[0] = B1;
	newBs[1] = B2;
	return;
}

// get a list of the weak edge.  
void DMWT::GetAllWeakEdge(Boundary<MAXK>& B,vector<int>& weSet){
	if(B.isEmpty()) 
		return;
	int e;
	int vi [2] = {0,0};
	int vj [2] = {0,0};
	B.getAllE(mEdges);
	for (int p=0; p<B.getSize(); p++){
		vi[0] = edgeInfoList[mEdges[p]]->v1;
		vi[1] = edgeInfoList[mEdges[p]]->v2;
		for (int q=p+1; q<B.getSize(); q++){
			vj[0] = edgeInfoList[mEdges[q]]->v1;
			vj[1] = edgeInfoList[mEdges[q]]->v2;
			for (int i=0; i<2; i++){
				for (int j=0; j<2; j++){
					e = ehash[vi[i]][vj[j]];
					if (e!=-1 && !isBoardE(e)){
						int b;
						for (b=0; b<B.getSize(); b++){
							if (e==mEdges[b]) break;
						}
						if (b==B.getSize()){
							weSet.push_back(e);
						}
					}
				}
			}
		}
	}
	return;
}

// get the bitkey of allItem in canItem. the key has the same length of canItem, if canItem[i] exist in allItem, key[i] is set to 1
__int64 DMWT::GetBitKey(const vector<int>& canItem,const vector<int>& allItem){
	__int64 key = 0, mask = 1, one = 1;
	int e, pos=0, size=(int)canItem.size()-1;
	size_t canS=canItem.size();
	size_t allS=allItem.size();
	size_t i,j;
	for (i=0; i<canS; i++){
		e = canItem[i];
		for (j=0; j<allS; j++){
			if (e == allItem[j]){
				break;
			}
		}
		if (j!=allS){
			mask = one << (size-pos);
			key = key | mask;
		}
		pos++;
	}
	return key;
}

// recover the sub list from list, according to the key
void DMWT::RecoverListFromBitKey(const vector<int>& list, __int64 key,vector<int>& sublist ){
	sublist.clear();
	int size = (int)list.size();
	__int64 mask = 1, exist;
	int item, pos = 0;
	for (int i=size-1; i>=0; i--){
		item = list[i];
		exist = (key >> pos ) & mask;
		if (exist) {
			sublist.push_back(item);
		}
		pos++;
	}
}

//==================================Preprocess============================//

// verify whether edge(v1,v2) is a bd edge
// v1<v2
__forceinline bool DMWT::isBDEdge(int v1, int v2){
	int c1 = vertex2curve[v1];
	int c2 = vertex2curve[v2];
	int start = curveInfo[c1].start;
	int end = curveInfo[c1].end;
	return ((c1==c2)&&((v1==start&&v2==end)||v2==v1+1));
}

// scan tris onece, 
// sort the vertex index in tris, from small to large
// build ehash, ehash[v1][v2] = e, v1<v2
// build ehashsize, the # of tri for each edge
// return # of edges
int DMWT::scanTrianglesOnce(){
	int edgenum=0; int min, max, sum = 0, mid, v, v1, v2;
	for(int t=0; t<numoftris; t++){
		// sort vertices to min-mid-max
		min = INT_MAX; max = -1; sum = 0;
		for(int i=0; i<3; i++){
			v = triangle(t,i);
			min = min < v ? min : v;
			max = max > v ? max : v;
			sum += v;
		}
		mid = sum - min - max;
		// rerangle vertice in tris
		triangle(t,0) = min; triangle(t,1) = mid; triangle(t,2) = max;
		// build ehash
		for(int i=0; i<3; i++){
			if (i==2) { v1 = triangle(t,0); v2 = triangle(t,2); }
			else { v1 = triangle(t,i); v2 = triangle(t,i+1); }
			if(ehash[v1][v2]==-1) { ehash[v1][v2]=edgenum; ehash[v2][v1]=edgenum; edgenum++;}
			ehashSize[v1][v2]++;
		}
	}
	return edgenum;
}

void DMWT::buildList(){
	if (badInput) return;
	//initialize ehash
	ehash = new int* [numofvertices];
	ehashSize = new int* [numofvertices];
	for (int i=0; i<numofvertices; i++){
		ehash[i] = new int[numofvertices];
		ehashSize[i] = new int[numofvertices];
		for (int j=0; j<numofvertices; j++){
			ehash[i][j] = -1; 
			ehashSize[i][j] = 0;
		}
	}
	triangleInfoList = new TriangleInfo * [numoftris];
	for(int i=0; i<numoftris; i++){
		triangleInfoList[i] = new TriangleInfo();
	}
	//scan triangle list once to assign index of edges
	numofedges = scanTrianglesOnce();
	//create edgeInfoList & triangleInfoList
	startE = ehash[0][1];
	edgeInfoList = new EdgeInfo * [numofedges];
	for(int i=0; i<numofedges; i++){
		edgeInfoList[i] = new EdgeInfo();
	}

	int v1, v2, ei=-1, nTri, posTri, nBDedge;
	//scan all triangles, setup triangleInfoList and most of edgeInfoList except BiTri information
	//initialize all edgeInfoList
	for(int t=0; t<numoftris; t++){
		nBDedge = 0;
		for(int i=0; i<3; i++){
			if (i==2) { v1 = triangle(t,0); v2 = triangle(t,2); }
			else { v1 = triangle(t,i); v2 = triangle(t,i+1); }
			ei = ehash[v1][v2];
			EdgeInfo * curE = edgeInfoList[ei];
			//after process one edge, set ehashSize[v1][v2] to 
			//the next slot for inserting a triangle
			if(curE->nTri==-1){
				nTri = ehashSize[v1][v2];
				curE->v1 = v1;
				curE->v2 = v2;
				curE->nTri = nTri;
				curE->isBD = isBDEdge(v1,v2);
				curE->sgMeasure = measureEdge(v1,v2);
				curE->tris = new int[nTri];
				curE->edgeInd  = new char[nTri];
				ehashSize[v1][v2] = 0;
			}
			if(curE->isBD) nBDedge++;
			posTri = ehashSize[v1][v2]; 
			TinE(ei, posTri) = t;
			EIinE(ei, posTri) = i;
			EinT(t,i) = ei;
			TIinT(t,i) = posTri;
			ehashSize[v1][v2] ++;
		}
		triangleInfoList[t]->nBDedge = nBDedge;
		triangleInfoList[t]->sgMeasure = measureTri(triangle(t,0), triangle(t,1), triangle(t,2));
	}
	maxTriPerEdge = 0;
	for (int ei=0; ei<numofedges; ei++){
			EdgeInfo * curE = edgeInfoList[ei];
			nTri = curE->nTri;
			if (nTri > maxTriPerEdge){
				maxTriPerEdge = nTri;
			}
	}
	// calculate bi-tri metric
	int nBiTri, t1, t2, count=0, vsum, p, q, nor;
	if (useBiTri) {
		for (int ei=0; ei<numofedges; ei++){
			EdgeInfo * curE = edgeInfoList[ei];
			v1 = curE->v1; v2 = curE->v2; vsum = v1 + v2;
			nTri = curE->nTri;
			nBiTri = nTri * (nTri-1) / 2;
			curE->biMeasure = new float[nBiTri];
			count = 0;
			for (int i=0; i<nTri-1; i++){
				for (int j=i+1; j<nTri; j++){
					t1 = curE->tris[i];
					t2 = curE->tris[j];
					p = triangle(t1,0) + triangle(t1,1) + triangle(t1,2) - vsum;
					q = triangle(t2,0) + triangle(t2,1) + triangle(t2,2) - vsum;
					curE->biMeasure[count] = measureBiTri(v1, v2, p, q);
					count++;
				}
			}
		}
	}
	// calculate bd-tri metric
	if (useNormal){
		for (int ei=0; ei<numofedges; ei++){
			EdgeInfo * curE = edgeInfoList[ei];
			v1 = curE->v1; v2 = curE->v2; vsum = v1 + v2;
			if (v2-v1 == 1){
				nor = v1;
			} else {
				nor     = v1;
				int tmp = v1;
				v1      = v2;
				v2      = tmp;
			}
			nTri = curE->nTri;
			curE->bdMeasure = new float[nTri];

			for (int i=0; i<nTri; i++){
				t1 = curE->tris[i];
				p = triangle(t1,0) + triangle(t1,1) + triangle(t1,2) - vsum;
				if (curE->isBD){
					curE->bdMeasure[i] = measureTriBd(v1, v2, p, nor);
				} else {
					curE->bdMeasure[i] = 0.0;
				}
			}
		}
	}

	for (int i=0; i<numofvertices; i++){
		// delete [] ehash[i];
		delete [] ehashSize[i];
	}
	// delete [] ehash;
	delete [] ehashSize; 
}


//==================================File IO===================================//

void DMWT::readCurveFile(const char* file){
	// extract the name of the curve
	filename = new char[300]; 
	strcpy_s(filename,300,file);
	char* dot=strrchr(filename,'.'); *dot = '\0';
	
	// read curves
	std::ifstream reader(file, std::ifstream::in);
	if ( !reader.good()) exit(1);
	reader >> numofcurves;
	reader >> numofvertices;
	orgnumofvertices = numofvertices;
	curveInfo    = new CurveInfo[numofcurves];
	vertex2curve = new int[numofvertices];
	vertices = new double[3 * numofvertices];
	tempC.reserve(numofvertices);
	tempOrgC.reserve(numofvertices);
	tempAdj.reserve(numofvertices);
	radius.reserve(numofvertices);
	orgradius.reserve(numofvertices);
	cliped.reserve(numofvertices);
	newEdge.push_back(0); newEdge.push_back(0);
	newAdj.push_back(0); newAdj.push_back(0);
	newClip.push_back(0); newClip.push_back(0);
	int pt = 0;
	double x,y,z;
	for (int i=0; i<numofcurves; i++){
		reader >> nPt;
		curveInfo[i].start = pt;
		curveInfo[i].end   = pt + nPt - 1;
		curveInfo[i].len = nPt;
		for (int j=0; j<nPt; j++) {
			reader >> x >> y >> z;
			tempC.push_back(Point3(x,y,z));
			tempOrgC.push_back(Point3(x,y,z));
			newAdj[0] = pt+1; newAdj[1] = pt-1; 
			tempAdj.push_back(newAdj);
			vertex2curve[pt] = i;
			vertices[pt*3] = x;
			vertices[pt*3+1] = y;
			vertices[pt*3+2] = z;
			pt++;
		}
		//tempAdj[vi][0] pt index larger than vi; [1] smaller than vi
		tempAdj[curveInfo[i].start][1] = curveInfo[i].end;
		tempAdj[curveInfo[i].end][0] = curveInfo[i].start;
	}
	reader.close();
	if (numofvertices < 3 ){
		exit(1);
	}
}

void DMWT::readNormalFile(const char* file){

	// read curve normals
	std::ifstream reader(file, std::ifstream::in);
	if ( !reader.good()) exit(1);
	reader >> numofcurves;
	reader >> numofnormals;
	newNorm.push_back(0); newNorm.push_back(0);
	//MING: clear the mem after EP. 
	tempNorm.reserve(numofnormals);
	tempAdjNorm.reserve(numofvertices);
	int pt = 0;
	double x,y,z;
	for (int i=0; i<numofcurves; i++){
		reader >> nPt;
		for (int j=0; j<nPt; j++) {
			reader >> x >> y >> z;
			tempNorm.push_back(Vector3(x,y,z));
			newNorm[0] = pt; newNorm[1] = pt-1;
			tempAdjNorm.push_back(newNorm);
			pt++;
		}
		//tempAdjNorm[vi][0] itself; [1] index smaller than vi
		tempAdjNorm[curveInfo[i].start][1] = curveInfo[i].end;
	}
	reader.close();

}

void DMWT::genDTTriCandidates(){
	if (numofvertices == 3){
		numoftris = 1;
		tris = new int[numoftris*3];
		tris[0] = 0;
		tris[1] = 1;
		tris[2] = 2;
		return;
	}
	if (!EdgeProtection()){
		return;
	}
	//set tetgenio point set input
	tetgenIn.numberofpoints = numofvertices;
	if (isDeGen){
		tetgenIn.pointlist = pb_vertices;
	}else{
		tetgenIn.pointlist = vertices;
	}
	//call tetgen
	tetrahedralize("f", &tetgenIn, &tetgenOut);
	numoftris = tetgenOut.numberoftrifaces;
	tris = tetgenOut.trifacelist;
#if (_SaveTetTri==1)
	tetgenOut.save_faces(filename);
#endif
}
void DMWT::genAllTriCandidates(){
	numoftris = numofvertices*(numofvertices-1)*(numofvertices-2)/6;
	tris = new int[numoftris*3];
	int counter=0;
	for(int i=0; i<numofvertices-2; i++){
		for(int j=i+1; j<numofvertices-1; j++){
			for(int k=j+1; k<numofvertices; k++){
				tris[counter*3]=i;
				tris[counter*3+1]=j;
				tris[counter*3+2]=k;
				counter++;
			}
		}
	}
}

/// read in .tri file, fill following datastructure:
/// int numoftris;
/// int* tris;
void DMWT::readCanTFile(const char* file){
	// read candidate triangles 
	std::ifstream reader(file, std::ifstream::in);
	if ( !reader.good()) exit(1);
	reader >> numoftris;
	tris = new int[numoftris*3];
	for (int i=0; i<numoftris; i++){
		reader >> tris[i*3] >>tris[i*3+1] >>tris[i*3+2];
	}
	reader.close();
}

void DMWT::saveTiling(){
	char * tilefile = new char[300];
	strcpy_s(tilefile,300,filename);
	if (useDT) {
		strcat_s(tilefile,300,"_dyn");
	} else {
		strcat_s(tilefile,300,"_all");
	}
	if (saveObj){
		strcat_s(tilefile,300,".obj");
		saveTilingObj(tilefile);
	}
	delete [] tilefile;
}

// save tiling
void DMWT::saveTiling(char* tilefile){
	int n = (int) optTile.size();
	std::ofstream writer(tilefile, std::ofstream::out);
	if (!writer.good()) exit(1);
	writer << "{";
	for(int i=0; i<n-1; i++){
		writer << optTile[i]+1 <<",";
	}
	if (n==0) writer << "}\n";
	else writer << optTile[n-1]+1 << "}\n";
	writer.close();
}
void DMWT::saveTilingObj(char* tilefile){
	int n = (int) optTile.size();
	std::ofstream writer(tilefile, std::ofstream::out);
	if (!writer.good()) exit(1);
	writer << "# OBJ File Generated by MMWT\n";
	writer << "# Vertices: " << numofvertices << "\n";
	writer << "# Faces: " << n << "\n";
	// write vertices
	for(int i=0; i<numofvertices; i++){
		writer << "v "<< vertex(i,0) << " " << vertex(i,1) << " " << vertex(i,2) <<"\n";
	}
	// write faces
	int t = 0;
	for (int i=0; i<n; i++){
		t = optTile[i];
		writer << "f " << triangle(t,0)+1 << " " << triangle(t,1)+1 << " " << triangle(t,2)+1 <<"\n";
	}
	writer.close();
}
//==================================Weight Functions============================//

void DMWT::setWeights(float wtri, float wedge, float wbitri, float wtribd){
	weightTri = wtri;
	weightEdge = wedge;
	weightBiTri = wbitri;
	weightTriBd = wtribd;
	if((weightBiTri==0.0f)&&(weightTriBd==0.0f)) {
		useBiTri = false;
		useMinMaxDihedral = false;
	} else {
		useBiTri = true;
		if (weightTriBd!=0.0f){
			weightBiTri=1.0f;
			useMinMaxDihedral = true;
		} else {
			useMinMaxDihedral = false;
		}
	}
}

__forceinline float DMWT::costTri(float measure){
	return weightTri*measure;
}

__forceinline float DMWT::costEdge(float measure){
	return weightEdge*measure;
}

__forceinline float DMWT::costBiTri(float measure){
	return weightBiTri*measure*measure;
}

__forceinline float DMWT::costTriBd(float measure){
	return weightTriBd*measure*measure;
}

__forceinline float DMWT::measureEdge(int v1, int v2){
	Point3 p1 = point(v1); Point3 p2 = point(v2);
	return (p2-p1).length();
}

__forceinline float DMWT::measureTri(int v1, int v2, int v3){
	Point3 p1 = point(v1); Point3 p2 = point(v2); 
	Point3 p3 = point(v3);
	return ((p2-p1)^(p3-p2)).length()/2;
}

__forceinline float DMWT::measureBiTri(int v1, int v2, int p, int q){
	Point3 p1 = point(v1); Point3 p2 = point(v2); 
	Point3 pp = point(p);  Point3 pq = point(q);
	Vector3 n1 = (p2-pp)^(p1-p2); n1.normalize();
	Vector3 n2 = (p2-p1)^(pq-p2); n2.normalize();
	float cosvalue = n1*n2; 
	cosvalue = cosvalue < 1.0 ? cosvalue : 1.0f-FLT_EPSILON;
	cosvalue = cosvalue > -1.0 ? cosvalue : -1.0f+FLT_EPSILON;
	return acos(cosvalue);
}

float DMWT::measureTriBd(int v1, int v2, int v3, int ni){
	Point3 p1 = point(v1); Point3 p2 = point(v2); Point3 p3 = point(v3);
	Vector3 n = Vector3(normals[ni*3],normals[ni*3+1],normals[ni*3+2]);
	Vector3 nt = (p2-p1)^(p3-p2); nt.normalize();
	float cosvalue = nt*n; 
	cosvalue = cosvalue < 1.0 ? cosvalue : 1.0f-FLT_EPSILON;
	cosvalue = cosvalue > -1.0 ? cosvalue : -1.0f+FLT_EPSILON;
	return acos(cosvalue);
}

__forceinline float DMWT::CostE(int e){
	float measure = measureEdge(edgeInfoList[e]->v1, edgeInfoList[e]->v2);
	return costEdge(measure);
}

__forceinline float DMWT::CostAllE(int t){
	return CostE(EinT(t,0)) + CostE(EinT(t,1)) + CostE(EinT(t,2));
}

__forceinline float DMWT::CostT(int t){
	float measure = measureTri(triangle(t,0), triangle(t,1), triangle(t,2));
	return costTri(measure);
}

__forceinline float DMWT::CostBI(int e, short fi1, short fi2){
	if(!useBiTri) return 0.0f;
	if(fi1==-1||fi2==-1) return 0.0f;
	if(fi1==fi2) return weightBiTri*PI;
	int ti1, ti2;
	if (fi1<fi2) {
		ti1 = fi1;
		ti2 = fi2;
	} else {
		ti1 = fi2;
		ti2 = fi1;
	}
	int pos = getBiPos(edgeInfoList[e]->nTri, ti1, ti2);
	return weightBiTri*edgeInfoList[e]->biMeasure[pos];
}

__forceinline float DMWT::CostBD(int e, int i){
	if (useNormal)
		return weightBiTri*edgeInfoList[e]->bdMeasure[i];
	return 0.0;
}

__forceinline int DMWT::getBiPos(int len, int i, int j){
	return (int)((i+1) * (2*len-(i+1)-1) / 2 + (j+1) - len - 1);
}

//==================================Help Func============================//
/// Check whether boundary B is a board edge
__forceinline bool DMWT::isBoardEdge(const Boundary<MAXK>& B){
	if (B.getSize()!=1) return false;
	int v1, v2, e;
	char s, d;
	e = B.getE(0);
	s = B.getS(0);
	d = B.getD(0);
	getSortV(e, s, v1, v2);
	if (nextV(v2,d)==v1 && isBoardE(e))
		return true;
	return false;
}

/// Check whether boundary B is formed by two meeting spanning edges
__forceinline bool DMWT::isMeetingSpanE(const Boundary<MAXK>& B){
	if (B.getSize()!=2) return false;
	int e1, e2;
	char s1, s2, d1, d2;
	e1   = B.getE(0);
	s1   = B.getS(0);
	d1   = B.getD(0);
	e2   = B.getE(1);
	s2   = B.getS(1);
	d2   = B.getD(1);
	if (d1==0 && d2==0 && e1==e2 && s1==-s2 && !isBoardE(e1))
		return true;
	return false;
}

/// Get the two end points of e in order of s
__forceinline void DMWT::getSortV(int e, int s, int& v1, int& v2){
	if (s == 1) {
		v1 = edgeInfoList[e]->v1;
		v2 = edgeInfoList[e]->v2;
	} else {
		v1 = edgeInfoList[e]->v2;
		v2 = edgeInfoList[e]->v1;
	}
}

/// Get the next vertex of v in direction d
__forceinline int DMWT::nextV(int v, int dir){
	int ci = vertex2curve[v];
	int lo = curveInfo[ci].start;
	int hi = curveInfo[ci].end;
	int nv = v + dir;
	if (nv > hi)
		return lo;
	if (nv < lo)
		return hi;
	return nv;
}

/// Get the next 2 span edges. e must be an edge of ff
__forceinline void DMWT::getSubActiveEdge(int e, int ff, int v1, int v2, int& ee1, int& ee2, short&ffi1, short&ffi2){
	int e1, e2, fi1, fi2;
	e1 = EinT(ff,0);
	fi1 = TIinT(ff,0);
	e2  = EinT(ff,1);
	fi2 = TIinT(ff,1);
	if (e1 == e){
		e1  = EinT(ff,2);
		fi1 = TIinT(ff,2);
	} else if (e2 == e) {
		e2  = EinT(ff,2);
		fi2 = TIinT(ff,2);
	}
	if (edgeInfoList[e1]->v1==v1||edgeInfoList[e1]->v2==v1){
		ee1  = e1;
		ee2  = e2;
		ffi1 = fi1;
		ffi2 = fi2;
	} else {
		ee1  = e2;
		ee2  = e1;
		ffi1 = fi2;
		ffi2 = fi1;
	}
}

/// check whether v is in the interval of pv->(d)->nv
/// v can be equal to pv or nv
/// pv, nv must be on the same curve
__forceinline bool DMWT::isWithinInterval(int pv, int nv, int v, int d){
	int ci = vertex2curve[pv];
	int start = curveInfo[ci].start;
	int end   = curveInfo[ci].end;
	switch(d)
	{
		case 0:
			if(pv==nv&&pv==v) return true;
			break;
		case 1:
			if(pv<nv && v>=pv && v<=nv) return true;
			if(pv>=nv && ((v>=start && v<=nv)||(v>=pv && v<=end))) return true;
			break;
		case -1:
			if(pv>nv && v>=nv && v<=pv) return true;
			if(pv<=nv && ((v>=start && v<=pv)||(v>=nv && v<=end))) return true;
			break;
		default:
			break;
	}
	return false;
}

/// Get the intervals of Boundary
__forceinline void DMWT::getIntervalsB(const Boundary<MAXK>& B, Interval<MAXK>& Intvs){
	if (B.getSize() <= 0) {
		return;
	}

	int e = B.getE(0);
	char s = B.getS(0);
	char d;
	int size = (int) B.getSize();
	int cv1, cv2, nv1, nv2;

	getSortV(e, s, cv1, cv2);
	for (int i = 0; i<size-1; i++){
		getSortV(B.getE(i+1), B.getS(i+1), nv1, nv2);
		d = B.getD(i);
		Intvs.append(cv2,nv1,d);
		cv1  = nv1;
		cv2  = nv2;
	}
	getSortV(B.getE(0), B.getS(0), nv1, nv2);
	d = B.getD(size-1);
	Intvs.append(cv2,nv1,d);
}

/// Find the max value among the 3
__forceinline float DMWT::findMax(float val1, float val2, float val3){
	if (val1 > val2){
		if (val1 > val3)
			return val1;
		return val3;
	} else {
		if (val2 > val3)
			return val2;
		return val3;
	}
}
/// Find the max value among the 4
__forceinline float DMWT::findMax(float val1, float val2, float val3, float val4){
	float val = findMax(val1,val3,val3);
	if (val > val4)
		return val;
	return val4;
}

/// Check whether Boudary B is valid
__forceinline bool DMWT::isInvalidB(Boundary<MAXK>& B){
	if (B.isEmpty()) return false;
	int e = B.getE(0);
	if (isBoardE(e) && !isBoardEdge(B)) 
		return true;
	return false;
}

//==================================Edge Protection============================//
/// edge protection : implementation of "Constrained Delaunay Tetrahedralizations and
/// Provably Good Boundary Recovery"
bool DMWT::EdgeProtection(){
	// perturbation on vertex if plain case is detected
	int perturbNum = 0;
	while (!passTetGen()){
		isDeGen = true;
		perturbPts(plainPTB);
		perturbNum++;
		if (perturbNum>500){
			badInput = true;
			cout<<"MWT: plain input, >500 perturbation"<<endl;
			return false;
		}
	}

	// edge protection for multiple curves
	if(!isProtected()){
		if (!badInput){
			protectCorner();
			insertMidPointsTetgen();
		}
	}

	if (badInput){
		return false;
	}
	getCurveAfterEP();

	#if SAVE_NEWCURVE
	saveCurve(filename,vertices,numofpoints);
	#endif

	return true;
}

/// For getCurveAfterEP(), to check whether the orientation changes
/// if so, normals should be fliped
bool DMWT::sameOrientation(const vector<int> & newCurve, int start, int len){
	int pos1 = -1,pos2 = -1,pos=start;
	while((pos1==-1||pos2==-1)&& pos<=len){
		if(newCurve[pos]==start + 1) pos1 = pos;
		if(newCurve[pos]==start + 2) pos2 = pos;
		pos++;
	}
	return pos2>pos1;
}

/// Recover ordered curve points after edge protection
bool DMWT::getCurveAfterEP(){
	int start,pt=-1, cur, pre, next;
	vector<int> newCurve;
	delete [] vertex2curve;
	vertex2curve = new int[numofvertices];
	for (int i=0; i<numofcurves; i++){
		start = curveInfo[i].start;
		curveInfo[i].start = (int)newCurve.size();
		newCurve.push_back(start);
		vertex2curve[++pt]=i;
		pre = start; cur = tempAdj[start][0];
		while(cur!=start){
			newCurve.push_back(cur);
			vertex2curve[++pt]=i;
			next = tempAdj[cur][0]+tempAdj[cur][1]-pre;
			pre = cur; cur = next;
		}
		curveInfo[i].end = (int)newCurve.size()-1;
		curveInfo[i].len = curveInfo[i].end - curveInfo[i].start + 1;
	}
	numofvertices = (int)newCurve.size();
	delete [] vertices;
	vertices = new double[3 * numofvertices];
	/*have the risk of flipping the orientation of the curve*/
	if (isDeGen){
		pb_vertices = new double[3 * numofvertices];
		for (int i=0; i<numofvertices; i++) {
			vertices[i*3+0] = tempOrgC[newCurve[i]][0];
			vertices[i*3+1] = tempOrgC[newCurve[i]][1];
			vertices[i*3+2] = tempOrgC[newCurve[i]][2];
			pb_vertices[i*3+0] = tempC[newCurve[i]][0];
			pb_vertices[i*3+1] = tempC[newCurve[i]][1];
			pb_vertices[i*3+2] = tempC[newCurve[i]][2];
		}
	} else {
		for (int i=0; i<numofvertices; i++) {
			vertices[i*3+0] = tempC[newCurve[i]][0];
			vertices[i*3+1] = tempC[newCurve[i]][1];
			vertices[i*3+2] = tempC[newCurve[i]][2];
		}
	}
	if(useNormal){
		normals = new float[3 * numofvertices];
		Vector3 norm;
		//bool sameOri;
		for (int ci=0; ci<numofcurves; ci++){
			/*sameOri = sameOrientation(newCurve, curveInfo[ci].start, curveInfo[ci].len);*/
			/*if(sameOri){
				for (int i=curveInfo[ci].start; i<=curveInfo[ci].end; i++) {
					norm = tempNorm[newCurve[i]];
					normals[i*3+0] = (float)norm[0];
					normals[i*3+1] = (float)norm[1];
					normals[i*3+2] = (float)norm[2];
				}
			}else{
				for (int i=curveInfo[ci].start; i<=curveInfo[ci].end; i++) {
					norm = tempNorm[newCurve[i]];
					normals[i*3+0] = (float)-norm[0];
					normals[i*3+1] = (float)-norm[1];
					normals[i*3+2] = (float)-norm[2];
				}
			}*/
			for (int i=curveInfo[ci].start; i<=curveInfo[ci].end; i++) {
				norm = tempNorm[newCurve[i]];
				normals[i*3+0] = (float)norm[0];
				normals[i*3+1] = (float)norm[1];
				normals[i*3+2] = (float)norm[2];
			}
		}
	}
	return true;
}

/// Insert mid points on edges that are not edge-protected
void DMWT::insertMidPointsTetgen(){
	int besize, p, q, j; Point3 p1, p2, newP, orgp1,orgp2, newOrgP;
	while(!isProtected()){
		if (badInput) break;
		besize = (int)badEdge.size();
		for (int i=0; i<besize; i++){
			p = badEdge[i][0]; j = badEdge[i][1]; 
			p1 = tempC[p]; q = tempAdj[p][j];
			p2 = tempC[q]; Vector3 vec = p2-Point3();
			newP = p1 + vec; newP*=0.5;

			if (isDeGen){
				orgp1 = tempOrgC[p]; 
				orgp2 = tempOrgC[q]; 
				Vector3 orgvec = orgp2-Point3();
				newOrgP = orgp1 + orgvec; newOrgP*=0.5;
				splitEdge(p,j,newP,newOrgP);
			} else {
				splitEdge(p,j,newP);
			}
		}
	}
}
/// Protect acute corners
void DMWT::protectCorner(){
	radius.clear();
	cliped.clear();
	vector<int> acuteList;
	
	for(int i=0; i<numofvertices; i++){
		newClip[0] = 0; newClip[1] = 0;
		cliped.push_back(newClip);
		radius.push_back(FLT_MAX);
	}
	if (isDeGen){
		orgradius.clear();
		for(int i=0; i<numofvertices; i++){
			orgradius.push_back(FLT_MAX);
		}
	}
	//gather acute corners
	for (int i=0; i<numofvertices; i++){
		if(getAngle(i,tempAdj[i][0],tempAdj[i][1])<hfPI){
			acuteList.push_back(i);
			cliped[i][0] = 1; cliped[i][1] = 1;
		}
	}
	//computer distance to other edges (circular protection)
	int acsize = (int)acuteList.size();
	int p1,p2,p3; int start, end; double dist;
	for (int i=0; i<acsize; i++){
		p1 = acuteList[i];
		for (int ci=0; ci<numofcurves; ci++){
			start=curveInfo[ci].start;
			end=curveInfo[ci].end;
			for (int j=start; j<=end; j++){
				p2 = j; p3 = j+1;
				//MING: not necessary to ensure the order
				if(j==end){ p2=start; p3=end; }
				if(p1==p2||p1==p3)	continue;
				dist = getPt2LineDist(p1,p2,p3);
				radius[p1] = min(dist,radius[p1]);
			}
		}
	}
	//computer cut positons min(dist,S,L/3,2D/3)
	for (int i=0; i<acsize; i++){
		double M = FLT_MAX;
		double L = FLT_MAX;
		double D = FLT_MAX;
		double S = FLT_MAX;

		double orgM = FLT_MAX;
		double orgL = FLT_MAX;
		double orgD = FLT_MAX;
		double orgS = FLT_MAX;
		int pos;
		p1 = acuteList[i];
		for (int j=0; j<2; j++){
			p2 = tempAdj[p1][j];

			M = (tempC[p1]-tempC[p2]).length();
			S = min(S,M);
			if(cliped[p1][j])	L=min(L,M);
			pos = 0;
			if(tempAdj[p2][0]!=p1) pos = 1;
			if(cliped[p2][pos])	D = min(D,M);

			if(isDeGen){
				orgM = (tempOrgC[p1]-tempOrgC[p2]).length();
				orgS = min(orgS,orgM);
				if(cliped[p1][j])	orgL=min(orgL,orgM);
				if(cliped[p2][pos])	orgD = min(orgD,orgM);
			}
		}
		radius[p1] = min(min(radius[p1],S),min(L/3.0,D*2.0/3.0));
		if (isDeGen){
			orgradius[p1] = min(min(orgradius[p1],orgS),min(orgL/3.0,orgD*2.0/3.0));
		}
	}
	//split cliped edges
	for (int i=0; i<acsize; i++){
		p1 = acuteList[i];
		for (int j=0; j<2; j++){
			p2 = tempAdj[p1][j];

			Point3 newP(tempC[p1]);
			Vector3 dir(tempC[p2]-tempC[p1]);
			dir.normalize();
			newP += dir*(radius[p1]*0.9);

			if (isDeGen){
				Point3 newOrgP(tempOrgC[p1]);
				Vector3 Orgdir(tempOrgC[p2]-tempOrgC[p1]);
				Orgdir.normalize();
				newOrgP += Orgdir*(orgradius[p1]*0.9);
				splitEdge(p1,j,newP,newOrgP);
			} else {
				splitEdge(p1,j,newP);
			}
		}
	}
}

/// Split Edge with newP
void DMWT::splitEdge(int p1, int p2ind, const Point3 & newP){
	//add new point
	int p2 = tempAdj[p1][p2ind];
	tempC.push_back(newP); numofvertices++;
	newAdj[0] = p1; newAdj[1] = p2;
	tempAdj.push_back(newAdj);
	if(useNormal){
		//and its adj norms
		int newNind = tempAdjNorm[p1][p2ind];
		Vector3 newN = tempNorm[newNind];
		tempNorm.push_back(newN);
		newNorm[0] = newNind; newNorm[1] = newNind;
		tempAdjNorm.push_back(newNorm);
	}
	tempAdj[p1][p2ind] = numofvertices-1;
	int pos = 0;
	if(tempAdj[p2][0]!=p1) pos = 1;
	tempAdj[p2][pos] = numofvertices-1;
}

/// Split Edge with newP in perturbed point set, with newOrgP in orginal point set
void DMWT::splitEdge(int p1, int p2ind, const Point3 &newP, const Point3 &newOrgP){
	//add new point
	int p2 = tempAdj[p1][p2ind];
	tempC.push_back(newP); numofvertices++;
	tempOrgC.push_back(newOrgP);
	newAdj[0] = p1; newAdj[1] = p2;
	tempAdj.push_back(newAdj);
	if(useNormal){
		//and its adj norms
		int newNind = tempAdjNorm[p1][p2ind];
		Vector3 newN = tempNorm[newNind];
		tempNorm.push_back(newN);
		newNorm[0] = newNind; newNorm[1] = newNind;
		tempAdjNorm.push_back(newNorm);
	}
	tempAdj[p1][p2ind] = numofvertices-1;
	int pos = 0;
	if(tempAdj[p2][0]!=p1) pos = 1;
	tempAdj[p2][pos] = numofvertices-1;
}
double DMWT::getAngle(int p1, int p2, int p3){
	double res = 0.0;
	Point3 v1 = tempC[p1];
	Point3 v2 = tempC[p2];
	Point3 v3 = tempC[p3];
	Vector3 vec1 = v1-v2; Vector3 vec2 = v1-v3;
	vec1.normalize(); vec2.normalize();
	res = vec1*vec2;
	res = res < -1.0 ? -1.0 : res;
	res = res > 1.0 ? 1.0 : res;
	return acos(res);
}
double DMWT::getPt2LineDist(int p1, int p2, int p3){
	double res = 0.0;
	Point3 v1 = tempC[p1];
	Point3 v2 = tempC[p2];
	Point3 v3 = tempC[p3];
	Vector3 vec12 = v2-v1;
	Vector3 vec23 = v3-v2;
	Vector3 vec13 = v3-v1;
	if(vec12*vec23>0)	return vec12.length();
	if(vec13*vec23<0)	return vec13.length();
	return sqrt(vec13*vec13-pow(1.0*(vec13*vec23),2.0)/(vec23*vec23));
}

/// Check whether the curve is alread edge-protected
bool DMWT::isProtected(){
	badEdge.clear();
	char ** used = new char*[numofvertices];
	for(int i=0; i<numofvertices; i++){
		used[i] = new char[2];
		used[i][0] = 0; used[i][1] = 0;
	}
	//set tetgenio point set input
	tetgenio tetIn;
	tetgenio tetOut;
	tetIn.numberofpoints = numofvertices;
	double * org_pts = new double[numofvertices*3];
	for (int i=0; i<numofvertices; i++){
		org_pts[i*3+0] = tempC[i][0];
		org_pts[i*3+1] = tempC[i][1];
		org_pts[i*3+2] = tempC[i][2];
	}
	tetIn.pointlist = org_pts;
	//call tetgen
	try{
		tetrahedralize("f", &tetIn, &tetOut);
	} catch (int e){
		cout<<"MWT: tetgen problem "<<e<<endl;
		badInput = true;
		delete [] used;
		return true;
	}
	int tmpTn = tetOut.numberoftrifaces;
	int * tmpTris = tetOut.trifacelist;
	int tri[3]; int p1, p2; int pos;
	for (int i=0; i<tmpTn; i++){
		tri[0] = tmpTris[i*3+0];
		tri[1] = tmpTris[i*3+1];
		tri[2] = tmpTris[i*3+2];
		for (int x=0; x<=1; x++){
			for (int y=x+1; y<=2; y++){
				p1 = tri[x]; p2 = tri[y];
				pos = -1;
				if(tempAdj[p1][0]==p2) pos = 0;
				if(tempAdj[p1][1]==p2) pos = 1;
				if(pos!=-1) used[p1][pos] = 1;
				p1 = tri[y]; p2 = tri[x];
				pos = -1;
				if(tempAdj[p1][0]==p2) pos = 0;
				if(tempAdj[p1][1]==p2) pos = 1;
				if(pos!=-1) used[p1][pos] = 1;
			}
		}
	}
	//save badedges for latter splitting
	for(int i=0; i<numofvertices; i++){
		for (int j=0; j<2; j++) {
			if((!used[i][j])&&i<tempAdj[i][j]){
				newEdge[0] = i; newEdge[1] = j;
				badEdge.push_back(newEdge);
			}
		}
	}
	int besize = (int)badEdge.size();
	if(besize>BADEDGE_LIMIT){
		badInput = true;
		delete [] used;
		return true;
	}
	delete [] used;
	if(besize>0) return false;
	return true;
}

/// Perturbe plain curve for tetgen
void DMWT::perturbPts(double ptb){
	for (unsigned int i=0; i<tempC.size(); i++){
		tempC[i].pertube(ptb*(i%5));
	}
}

/// Check whether input curve is spacial enough to pass tetgen
bool DMWT::passTetGen(){
	double* tempV = new double[tempC.size()*3];
	for (unsigned int i=0; i<tempC.size(); i++){
		tempV[i*3+0] = tempC[i][0];
		tempV[i*3+1] = tempC[i][1];
		tempV[i*3+2] = tempC[i][2];
	}
	tetgenio tetIn;
	tetgenio tetOut;
	tetIn.numberofpoints = (int)tempC.size();
	tetIn.pointlist = tempV;

	//call tetgen
	try{
		tetrahedralize("f", &tetIn, &tetOut);
	} catch (int e){
		cout<<"DMWT: tetgen error "<<e<<endl;
		return false;
	}

	int tmpTn = tetOut.numberoftrifaces;
	if (tmpTn == 0) 
		return false;
	return true;
}
//==================================Statistics============================//

void DMWT::statistics(){
}

