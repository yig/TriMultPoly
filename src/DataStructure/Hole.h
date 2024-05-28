// 
// Hole class
// 
#ifndef _HOLE_H_
#define _HOLE_H_

#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <time.h>
#include "SubKey.h"
#include "Configure.h"

class Hole {
public:
	std::vector<int> curves;

	Hole(){
	}
	Hole(Hole const & other){
		curves = other.curves;
	}
	Hole(std::vector<int> v){
		curves = v;
	}
	~Hole(){}
	Hole & operator=(Hole & H);
	int size(){
		return (int)curves.size();
	}
	int find(int ci);
	std::vector<std::pair<Hole,Hole>> partition();

	int countBit(__int64 bit, int n);
	void getSubKey(SubKey<MAXK> & sk);

	__forceinline void append(int item){
		curves.push_back(item);
	}
	__forceinline void erase(int hInd){
		if (hInd<0) return;
		curves.erase(curves.begin()+hInd);
	}
	__forceinline bool isEmpty(){
		return curves.empty();
	}
};


#endif
