// 
// Tile class
// 
#ifndef _Tile_H_
#define _Tile_H_
#include <float.h>
#include <stdlib.h>
#include <vector>
#include <unordered_map>

using namespace std;
extern int tileNum;

class OneTile {
public:
	int64_t WECode;
	float cost;
	int tri;
#if (_SaveTile==1)
	int next1;
	int next2;
#endif

	OneTile(){
		//tileNum++;
		WECode = 0; cost = FLT_MAX; tri = -1;
#if (_SaveTile==1)
		next1 = -1; next2 = -1;
#endif
	}
	OneTile(int64_t wec, float c, int t, int nxt1, int nxt2){
		//tileNum++;
		WECode = wec; cost = c; tri = t; 
#if (_SaveTile==1)
		next1 = nxt1; next2 = nxt2;
#endif
	}
	~OneTile(){}
};
/*
class Tile {
public:
	
	Tile();
	~Tile();

	vector<int> onetiles;
	//int min();

private:
};
*/
struct Tile {
	int start;
	int end;
	__forceinline void assign(int s, int e){
		start = s; end = e;
	}
};

/*class OneTile {
public:
	unsigned long WECode;
	float cost;
	int* tris;
	int nT;
	OneTile *next;
	OneTile();
	OneTile(int wec, float c, int* t, int n);
	~OneTile();
};

class Tile {
public:
	
	Tile();
	Tile(int wec, float c, int* t, int n);
	~Tile();

	int size;
	OneTile *head;
	OneTile *tail;

	OneTile * append(int wec, float c, int* t, int n);
	void append(OneTile * onetile);
	OneTile * min();
	Tile KeepOnlyTheMins();
	int getSize();
private:
};*/

#endif
