// 
// Implement simple linked list
// 
#ifndef _SUB_KEY_H_
#define _SUB_KEY_H_

template <int k>
class SubKey{
public:
	int key[k+1];
	SubKey(){
		/*for (int i=0; i<k+1; i++){
			key[i] = 0;
		}*/
	}
	~SubKey(){}
	void clear(){
		memset (key, 0, k+1);
	}
};

template <int k>
struct MySubKeyHash
{
	__forceinline size_t operator() (const SubKey<k> &x) const {

		const char *key = (const char *)&(x.key[0]);
		size_t nHash = 0;
		for (size_t i = 0; i < (sizeof(int)*(k+1)); i++) {
			nHash += key[i];
			nHash += (nHash << 10);
			nHash ^= (nHash >> 6);
		}
		nHash += (nHash << 3);
		nHash ^= (nHash >> 11);
		nHash += (nHash << 15);
		return nHash;

		/*int sum=0;
		for (int i=0; i<k+1; i++){
			sum +=x.key[i];
		}
		return sum;*/
	}
};

template <int k>
struct MySubKeyEq
{
	__forceinline bool operator() (const SubKey<k> &x, const SubKey<k> &y) const {
		for (int i=0; i<k+1; i++){
			if (x.key[i] != y.key[i])
				return false;
		}
		return true;
	}
};

#endif
