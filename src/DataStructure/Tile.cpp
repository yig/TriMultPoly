//#include "Tile.h"

//typedef boost::unordered_map< __int64, OneTile*> tilemap;
// ----------------OneTile
/*OneTile::OneTile(){
	tileNum++;
	WECode = 0; cost = FLT_MAX; tri = -1; //next1 = -1; next2 = -1;
}
OneTile::OneTile(__int64 wec, float c, int t, int nxt1, int nxt2){
	tileNum++;
	WECode = wec; cost = c; tri = t; //next1 = nxt1; next2 = nxt2;
}
OneTile::~OneTile(){}
*/
// ----------------Tile
//Tile::Tile(){}
//Tile::~Tile(){}
/*
int Tile::min(){
	if (onetiles.empty()) return -1;
	float minCost = FLT_MAX;
	/*OneTile * minTile = NULL;
	for (OneTile * curTile = head; curTile!=NULL; curTile=curTile->next){
		if (curTile->cost < minCost) {
			minCost = curTile->cost;
			minTile = curTile;
		}
	}*/
/*
	OneTile * minTile = NULL;
	OneTile * curTile = head;
	for (int i=0; i<size; i++){
		if (curTile->cost < minCost) {
			minCost = curTile->cost;
			minTile = curTile;
		}
		curTile = curTile->next;
	}
	return minTile;
}
*/