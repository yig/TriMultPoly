#include "Hole.h"

Hole & Hole::operator=(Hole & H){
	if (&H != this){
		curves = H.curves;
	}
	return *this;
}

int Hole::find(int ci){
	int ind = -1;
	int n = (int)curves.size();
	for (int i=0; i<n; i++){
		if (curves[i]==ci){
			ind = i;
			break;
		}
	}
	return ind;
}

void Hole::getSubKey(SubKey<MAXK> & sk){
	int n = (int)curves.size();
	int newK = 0;

	for (int i=0; i<n; i++){
		newK = newK ^ (1 << curves[i]);	
	}

	sk.key[MAXK] = newK;
}

std::vector<std::pair<Hole,Hole>> Hole::partition(){
	std::vector<std::pair<Hole,Hole>> res;
	const int n = (int)curves.size();
	int mask=1;
	unsigned long total = (unsigned long)pow(2.0f, n);
	for (unsigned long bit=0; bit<total; bit++){
		Hole H1, H2;
		for (int pos=0; pos<n; pos++){
			if (((bit >> pos) & mask) == 1) {
				H1.curves.push_back(curves[pos]);
			} else {
				H2.curves.push_back(curves[pos]);
			}
		}
		std::pair<Hole,Hole> sub(H1, H2);
		res.push_back(sub);
	}
	return res;
}

int Hole::countBit(int64_t bit, int n){
	int count = 0;
	int mask=1;
	for (int pos=0; pos<n; pos++){
		if (((bit >> pos) & mask) == 1) {
			count++;
		}
	}
	return count;
}



