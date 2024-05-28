#include "Boundary.h"
extern float timeSprintf;
extern float timeNew;
// ----------------OneB
/*OneB::OneB(int ee, char ss, short ff, char dd){
	e = ee; s = ss; f = ff; d = dd; next = NULL; pre = NULL;
}
OneB::~OneB(){}

// ----------------OneIntv
OneIntv::OneIntv(){}
OneIntv::OneIntv(int pV, int nV, char dd, char cc){
	pv = pV;
	nv = nV;
	d  = dd;
	c  = cc;
}
OneIntv::~OneIntv(){}
*/
/*template <int k>
void Boundary<k>::clear(){
}8/

// ----------------Boundary
/*
void Boundary::clear(){
	OneB * curB;
	OneB * nxtB = head;
	for (int i=0; i<size; i++){
		curB = nxtB;
		nxtB = curB->next;
		delete curB;
		//curB = nxtB;
		//nxtB = curB->next;
	}
}

void Boundary::append(int e, char s, short f, char d){
	if (size==0) {
		head = new OneB(e,s,f,d);
		head->next = head;
		head->pre = head;
		tail = head;
		size = 1;
	} else {
		OneB * newB = new OneB(e,s,f,d);
		newB->next = head;
		newB->pre = tail;
		tail->next = newB;
		head->pre = newB;
		tail = newB;
		size++;
	}
}*/
//get a uniform boundary
/*template <int k>
void Boundary::sort(){
	int minE = Infinity;
	int minI = -1;
	OneB * curB = head;
	OneB * startB;
	OneB * tmpB;
	int curE, curS;
	for (int i=0; i<size; i++){
		curE = curB->e;
		if (curE < minE){
			minE = curE;
			startB = curB;
		}
		curB = curB->next;
	}
	curS = startB->s;
	if (startB->s > 0){
		if (startB != head) {
			head = startB;
			tail = startB->pre;
		}
	} else {
		head = startB;
		tail = startB->next;
		curB = head;
		int headd = head->d;
		for (int i=0; i<size; i++){
			tmpB       = curB->next;
			curB->next = curB->pre;
			curB->pre  = tmpB;
			curB->s    = -curB->s;
			curB->d    = i==size-1? -headd : -curB->next->d;
			curB       = curB->next;
		}
	}
}*/
/*
// is empty boundary
bool Boundary::isEmpty()const{
	return size==0;
}
// return the point of #ind boundary
OneB * Boundary::getB(int ind){
	if (ind>=size) return NULL;
	OneB * curB = head;
	for (int i=0; i<ind; i++){
		curB = curB->next;
	}
	return curB;
}*/
// return a list of edges
/*template <int k>
void Boundary<k>::getAllE(std::vector<int>& edges)
{
	edges.clear(); 
	if (size==0) 
		return;
	clock_t start = clock(),end;

	
	edges.reserve(size);

	//int * allE = new int[size];
	
	for (int i=0; i<size; i++){
		edges.push_back(getE(i));
	}

	end = clock();
	timeNew += (float)(end-start)/CLOCKS_PER_SEC; 
	//return allE;

}*/
/*
// return a list of edges
int * Boundary::getAllE(){
	if (size==0) return NULL;
	clock_t start = clock(),end;

	int * allE = new int[size];
	OneB * curB = head;
	for (int i=0; i<size; i++){
		allE[i] = curB->e;
		curB = curB->next;
	}

	end = clock();
	timeNew += (float)(end-start)/CLOCKS_PER_SEC; 
	return allE;

}
// append partial of another boundary to current B
void Boundary::append(OneB * start, int num){
	if (num<=0) return;
	OneB * curB = start;
	for (int i=0; i<num; i++){
		append(curB->e, curB->s, curB->f, curB->d);
		curB = curB->next;
	}
}
*/
/*template <int k>
void Boundary::getSubKey(SubKey<k> &sk){
	

	clock_t start = clock(),end;
	for (int i=0; i<k; i++){
		sk.key[i] = OneBs[i];
	}
	end = clock();
	timeSprintf += (float)(end-start)/CLOCKS_PER_SEC; 
}*/