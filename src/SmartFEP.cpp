#include "header.h"
#include "analyze.h"
#include "placement.h"
#include "partition.h"
#include ".\szx\iplace.h"
#include "fileChecker.h"

int main() {
	//string instance = "motorbike";
	string instance = "propellerTip";
	//string instance = "flange";
	//string instance = "bullet";
	//string initialAlgorithm = "random";
	string initialAlgorithm = "projection";

	string basePath = ".\\data\\" + instance + "\\";
	string iterationsSaveFolder;
	int runIdx = 1;
	for (; ; ++runIdx) {
		iterationsSaveFolder = basePath + initialAlgorithm + "\\run" + to_string(runIdx);
		if (checkFolderExist(iterationsSaveFolder)) continue; else break;
	}
	system(("md " + iterationsSaveFolder).c_str());
	
	string finput = basePath + instance + ".in";
	string fouput = iterationsSaveFolder + "\\" + instance + ".out";
	string iterationsSaveFile = iterationsSaveFolder + "\\" + instance + "Iterations.txt";

	char *problem_input = (char *) finput.c_str();
	char *solution_output = (char *) fouput.c_str();

#ifdef LOCAL
	freopen(problem_input, "r", stdin);
	//freopen(solution_output, "w", stdout);
#endif // LOCAL
	///problem parameter
	int dimension;
	vector<int>cubeShape;
	vector<int>fabricShape;
	double alpha;
	double beta;
	vector<double>heatMap;

	///for placement
	int inv_sampling_step;
	int tileCnt, adapterCnt;
	vector<vector<int>>cnts;
	vector<vector<int>>pris_info;

	Analyze::readProblem(dimension, cubeShape, fabricShape, alpha, beta, heatMap);      //returns

	clock_t start = clock();

	Partition::solve(dimension, cubeShape, fabricShape, alpha, beta, heatMap,  //inputs
		inv_sampling_step, pris_info);                            //returns

#ifdef LOCAL
	printf("Partitioing OK, Spend Time: %.3lfs\n", (double)(clock() - start) / CLOCKS_PER_SEC);
#endif // LOCAL

	Analyze::analyze_prisms(dimension, cubeShape, fabricShape, alpha, beta, heatMap,    //inputs
		inv_sampling_step, pris_info, tileCnt, adapterCnt, cnts);   //returns

#ifdef LOCAL
	db("analyze_prisms OK!");
#endif // LOCAL

	int prePlacementCnt = 1;
	bool startWithCenter = false;
	int maxIteration = 100000;
	int searchRadiu = 5;
	unsigned randSeed = 0;
	bool localSearch = true;
	int visibleInterval = 1;
	double timeLimited = 1800.;

	vector<pair<int, int>>init_place;

	Vec<Pos> jobPosHints(calcJobPosHints(dimension, pris_info, tileCnt, adapterCnt, cnts));
	Vec<Pos> workerOfJobs(solveQapByLocalSearch(cnts, { fabricShape[0], fabricShape[1] }, jobPosHints));

	for (size_t j = 0; j < tileCnt + adapterCnt; ++j) {
		init_place.push_back({ workerOfJobs[j].x , workerOfJobs[j].y });
	}

	cnts.resize(tileCnt + adapterCnt);
	Placement::Solution bestSol = Placement::solve(inv_sampling_step, cubeShape, fabricShape, tileCnt, adapterCnt, cnts, pris_info,
		init_place,
		prePlacementCnt,
		startWithCenter,
		maxIteration,
		searchRadiu,
		randSeed,
		localSearch,
		visibleInterval,
		timeLimited,
		initialAlgorithm,
		iterationsSaveFile);

#ifdef LOCAL
	Placement::save_solution(solution_output, bestSol.place_of_tile);
#endif // LOCAL

#ifndef LOCAL
	Placement::outputSolution(bestSol.place_of_tile);
#endif // LOCAL

	return 0;
}

