#include "DMWT.h"
#include <string>
#include <time.h>
// #include "windows.h"
// #include "psapi.h"
// #pragma comment(lib, "psapi.lib") 
using namespace std;
float timeNew;
int tileNum;

void pause(){
	cout<<"pause......"<<endl;
	int tmp;
	cin >> tmp;
	cout<<"pause end."<<endl;
}

// for traigulation cmd line
// usage: TMP.exe <curveName> <useDT> <useMinSet> <areaWeight> <edgeWeight> <dihedralWeight> <useMinMaxDihedral> <saveObj> <useNormal>
// example: TMP.exe monkey 1 1 0 0 1 0 1 0
void suf_cmd( int argc, char *argv[] ){

	if (argc!=10) {
		cout<<"\nUsage: \nTMP.exe <curveName> <useDT> <useMinSet> <areaWeight> <edgeWeight> <dihedralWeight> <useMinMaxDihedral> <saveObj> <useNormal>"<<endl;
		cout<<"\nExample: \nTMP.exe monkey 1 1 0 0 1 0 1 0"<<endl;
		return;
	}

	std::string testfile = argv[1];
	std::string curvefile;
	std::string normalfile;
	std::string canTfile;
	curvefile = testfile + ".curve";
	normalfile = testfile + ".normal";
	canTfile = testfile + ".tri";

	bool useDT=true;
	const char* tagDT = argv[2];
	useDT=(strcmp(tagDT,"1")==0);

	bool filterWE=false;
	const char* tagFWE = argv[3];
	filterWE=(strcmp(tagFWE,"1")==0);

	float weightTri = stof(argv[4]);
	float weightEdge = stof(argv[5]);
	float weightBiTri = stof(argv[6]);
	float weightTriBd = stof(argv[7]);

	bool saveObj=false;
	const char* tagObj = argv[8];
	saveObj=(strcmp(tagObj,"1")==0);

	bool useNormal=false;
	const char* tagNormal = argv[9];
	useNormal=(strcmp(tagNormal,"1")==0);

	DMWT * myDMWT;
	if (useNormal){
		myDMWT = new DMWT(curvefile.c_str(), normalfile.c_str(), useDT, 0, 1, filterWE, saveObj);
	} else {
		myDMWT = new DMWT(curvefile.c_str(), useDT,  0, 1, filterWE, saveObj);
	}
	myDMWT->setWeights(weightTri, weightEdge, weightBiTri, weightTriBd);

	// ---------build lists--------
	clock_t start = clock(),end;
	myDMWT->buildList();
	end = clock();
	float timeList = (float)(end-start)/CLOCKS_PER_SEC; 

	// ---------tile--------
	start = clock();
	myDMWT->tile();
	end = clock();
	float timeTile = (float)(end-start)/CLOCKS_PER_SEC;

	// ---------statistics--------
	cout<<"============================="<<endl;
	cout<<"Curve Name: "<<testfile<<endl;
	cout<<"----------------------------"<<endl;
	cout<<"<N> numCurves: "<<myDMWT->numofcurves<<endl;
	cout<<"<N> numPts: "<<myDMWT->orgnumofvertices<<endl;
	cout<<"<T> time: "<<timeList+timeTile<<" (secs)\n";
// 	HANDLE hProc = GetCurrentProcess();
// 	PROCESS_MEMORY_COUNTERS_EX info;
// 	info.cb = sizeof(info);
// 	BOOL okay = GetProcessMemoryInfo(hProc, (PROCESS_MEMORY_COUNTERS*)&info, info.cb);
// 	size_t memusage = info.WorkingSetSize;
// 	cout<<"<M> memory: "<<memusage/1024<<" KB"<<endl;
	cout<<"----------------------------"<<endl;
	cout<<"cost: "<<myDMWT->optCost<<endl;
	cout<<"numSub: "<<myDMWT->numSub<<endl;
	cout<<"============================="<<endl;
}

// for edge protection cmd line
// usage: TMP.exe -ep <curveName> <useNormal>
// example: TMP.exe monkey 0
void ep_cmd( int argc, char *argv[] ) {

	if (argc!=4) {
		cout<<"\nUsage: \nTMP.exe <curveName> <useNormal>"<<endl;
		cout<<"\nExample: \nTMP.exe monkey 0 "<<endl;
		return;
	}

	std::string testfile = argv[2];
	std::string curvefile;
	std::string normalfile;
	curvefile = testfile + ".curve";
	normalfile = testfile + ".normal";

	bool useDT=true;
	bool filterWE=false;
	bool saveObj=false;
	bool useNormal=false;
	const char* tagNormal = argv[3];
	useNormal=(strcmp(tagNormal,"1")==0);

	float weightTri = 1.0;
	float weightEdge = 0.0;
	float weightBiTri = 0.0;
	float weightTriBd = 0.0;

	DMWT * myDMWT;
	if (useNormal){
		myDMWT = new DMWT(curvefile.c_str(), normalfile.c_str(), useDT, 0, 1, filterWE, saveObj);
	} else {
		myDMWT = new DMWT(curvefile.c_str(), useDT,  0, 1, filterWE, saveObj);
	}
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

int main( int argc, char *argv[] ){

	timeNew = 0.0;
	tileNum = 0;

	if (cmdOptionExists(argv, argv+argc, "-ep")) {
		cout << "Start edge protection..."<<endl;
		ep_cmd(argc, argv);
	} else {
		suf_cmd(argc, argv);
	}
	return 0;
}
