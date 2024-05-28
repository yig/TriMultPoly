// 
// Boundary class
// 
#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "SubKey.h"
#include "Configure.h"
#include <boost/lexical_cast.hpp>
#include <string>
#include <sstream>
#include <time.h>
#include <vector>

#define Infinity 2000000000

template <int k>
class Boundary {
public:
	
	Boundary(){ size=0; }
	~Boundary(){}

	//void clear(){}
	void getAllE(std::vector<int>& edges){
		edges.clear(); 
		if (size==0) 
			return;

		for (int i=0; i<size; i++){
			edges.push_back(getE(i));
		}
	}
	__forceinline int  getE(int i) const{
		return ((OneBs[i] >> 11) & 0x1FFFFF) - 1;
	}
	__forceinline char getS(int i) const{
		return ((OneBs[i] >> 10) & 0x1) == 0? -1 : 1;
	}
	__forceinline short getF(int i) const{
		return ((OneBs[i] >> 2) & 0xFF) - 1;
	}
	__forceinline char getD(int i) const{
		return (OneBs[i] & 0x3) - 1;
	}

	__forceinline void setS(int i, char s){
		char S = s == -1? 0 : 1;
		OneBs[i] = (OneBs[i] & 0xFFFFFBFF) | (S << 10);
	}
	__forceinline void flipS(int i){
		OneBs[i] = OneBs[i] ^ 0x00000400;
	}
	__forceinline void setD(int i, char d){
		OneBs[i] = (OneBs[i] & 0xFFFFFFFC) | (d+1);
	}

	__forceinline void  setOneB(int i, int e, char s, short f, char d){
		int S = s == -1? 0 : 1;
		OneBs[i] = ((e+1) << 11) | (S << 10) | ((f+1) << 2) | (d+1);
	}
	__forceinline int  getOneB(int i){
		return	OneBs[i];
	}
	__forceinline void  append(int e, char s, short f, char d){
		setOneB(size, e,s,f,d);
		size++;
	}
	__forceinline void append(Boundary<k> &orgB, int ind, int num){
		if (num <=0) return;
		int oSize = orgB.getSize();
		if (num > oSize) return;
		for (int i=0; i<num; i++){
			OneBs[size] = orgB.getOneB((ind+i)%oSize);
			size++;
		}
	}
	__forceinline void getOneB(int i, int &e, char &s, short &f, char &d) const{
		e = getE(i);
		s = getS(i);
		f = getF(i);
		d = getD(i);
	}

	__forceinline bool isEmpty() const {
		return size==0;
	}
	__forceinline void sort(){
		if (size==0) {
			return;
		}
		int minE = Infinity;
		int minI = -1;
		int curE;
		char curS, curD;
		for (int i=0; i<size; i++){
			curE = getE(i);
			if (curE < minE){
				minE = curE;
				minI = i;
			}
		}
		curS = getS(minI);
		if (curS > 0){
			if (minI !=0 ){
				PosShift(minI);
			}
		} else {			
			NegShift(minI);
			char headd = getD(0);
			for (int i=0; i<size-1; i++){
				flipS(i);
				curD = getD(i+1);
				setD(i,-curD);
			}
			flipS(size-1);
			setD(size-1,-headd);
		}
	}

	__forceinline void Reverse(int b, int e){
		for (; b<e; b++, e--){
			int tmp = OneBs[e];
			OneBs[e] = OneBs[b];
			OneBs[b] = tmp;
		}
	}

	__forceinline void PosShift(int i){
		Reverse(i, size-1);
		Reverse(0, i-1);
		Reverse(0, size-1);
	}
	__forceinline void NegShift(int i){
		Reverse(0,i);
		Reverse(i+1, size-1);
	}

	__forceinline void getSubKey(SubKey<k> &sk) const {
		for (int i=0; i<size; i++){
			sk.key[i] = OneBs[i];
		}
		for (int i=size; i<k; i++){
			sk.key[i] = 0;
		}
	}

	__forceinline char getSize() const {
		return size;
	}

private:
	int OneBs[k];
	char size;	
};

template <int k>
class Interval {
public:
	
	Interval(){ size=0; }
	~Interval(){}

	__forceinline int getPV(int i) const{
		return ((intvs[i] >> 17) & 0x1FFFFF);
	}
	__forceinline int getNV(int i) const{
		return ((intvs[i] >> 2) & 0x7FFF);
	}
	__forceinline char getD(int i) const{
		return (intvs[i] & 0x3) - 1;
	}

	__forceinline void  setOneIntv(int i, int pv, int nv, char d){
		intvs[i] = (pv << 17) | (nv << 2) | (d+1);
	}
	__forceinline void  append(int pv, int nv, char d){
		setOneIntv(size, pv, nv, d);
		size++;
	}
	__forceinline void getOneIntv(int i, int &pv, int &nv, char &d) const{
		pv = getPV(i);
		nv = getNV(i);
		d = getD(i);
	}

	__forceinline bool isEmpty() const {
		return size==0;
	}

	__forceinline char getSize() const {
		return size;
	}

private:
	int intvs[k];
	char size;
};


#endif
