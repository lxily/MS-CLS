//#include <iostream>
//#include <fstream>
#include <algorithm>
#include <cmath>

#include "Typedef.h"
#include "ConsecutiveIdSet.h"
#include "Random.h"
#include "Timer.h"
#include "geometry.h"


using namespace std;
using namespace hust;
using namespace hust::util;


void reprint(Str s) { cerr << s; fill(s.begin(), s.end(), '\b'); cerr << s; }


struct Pos {
    ID x;
    ID y;
	ID z;

	Pos() {}
	Pos(ID _x) :x(_x), y(0), z(0) {}
	Pos(ID _x, ID _y):x(_x), y(_y), z(0) {}
	Pos(ID _x, ID _y, ID _z) : x(_x), y(_y), z(_z) {}

    friend bool operator==(const Pos &l, const Pos &r) { return (l.x == r.x) && (l.y == r.y) && (l.z == r.z); }
    friend bool operator!=(const Pos &l, const Pos &r) { return (l.x != r.x) || (l.y != r.y) || (l.z != r.z); }
};
enum { DistanceAmp = 1000 };

using CalcDistance = Func<Int(const Pos &w0, const Pos &w1)>;
using OnWorker = Func<bool(const Pos &w)>;

struct SpiralScanner {
    Pos o;
    ID r;
    Pos d;

    SpiralScanner(const Pos &origin, ID initRadius = 1) : o(origin), r(initRadius), d({ r, 0 }) {}

    bool next(OnWorker onWorker) {
        if (onWorker({ o.x + d.x, o.y + d.y })) { return false; }
        if (onWorker({ o.x - d.x, o.y - d.y })) { return false; }
        if (onWorker({ o.x + d.y, o.y - d.x })) { return false; }
        if (onWorker({ o.x - d.y, o.y + d.x })) { return false; }
        if ((--d.x) > 0) { ++d.y; } else { d = { ++r, 0 }; }
        return true;
    }

    void scan(OnWorker onWorker) { while (next(onWorker)) {} }
};

void normalizeJobPosHints(const Pos &workerMat, Vec<Pos> &jobPosHints) {
    Pos leftBottom = jobPosHints.front();
    Pos rightTop = jobPosHints.back();
    for (auto p = jobPosHints.begin(); p != jobPosHints.end(); ++p) {
        math::updateMin(leftBottom.x, p->x);
        math::updateMin(leftBottom.y, p->y);
        math::updateMax(rightTop.x, p->x);
        math::updateMax(rightTop.y, p->y);
    }
    Pos jobMat = { rightTop.x - leftBottom.x + 1, rightTop.y - leftBottom.y + 1 };
    Real scaleX = static_cast<Real>(workerMat.x) / jobMat.x;
    Real scaleY = static_cast<Real>(workerMat.y) / jobMat.y;
    for (auto j = jobPosHints.begin(); j != jobPosHints.end(); ++j) {
        (j->x -= leftBottom.x) *= scaleX;
        (j->y -= leftBottom.y) *= scaleY;
    }
}

void padJobPosHints(const Pos &workerMat, Vec<Pos> &jobPosHints, Random &r) {
    ID workerNum = workerMat.x * workerMat.y;
    ID originalJobNum = static_cast<ID>(jobPosHints.size());
    jobPosHints.resize(workerNum);
    for (ID j = originalJobNum; j < workerNum; ++j) {
        jobPosHints[j] = { r.pick(workerMat.x), r.pick(workerMat.y) };
    }
}

Vec<Pos> solveQapByLocalSearch(Vec<Vec<ID>> &jobAdjList, const Pos &workerMat,
    CalcDistance normDistance, Vec<Pos> &jobPosHints, CalcDistance distance) {
    Random r;

    ID workerNum = workerMat.x * workerMat.y;
    ID originalJobNum = static_cast<ID>(jobAdjList.size());
    ID jobNum = workerNum;
    jobAdjList.resize(jobNum);
    
    //normalizeJobPosHints(workerMat, jobPosHints);
    padJobPosHints(workerMat, jobPosHints, r);

    Vec<Pos> workerOfJobs(jobNum);
    Arr2D<ID, ID> jobForWorkers(workerMat.x, workerMat.y, InvalidId);
    auto isWorker = [&](const Pos &w) { return (w.x >= 0) && (w.x < workerMat.x) && (w.y >= 0) && (w.y < workerMat.y); };
    auto isFreeWorker = [&](const Pos &w) { return jobForWorkers[w.x][w.y] == InvalidId; };
    auto assignJob = [&](ID j, const Pos &w) {
        workerOfJobs[j] = w;
        jobForWorkers[w.x][w.y] = j;
    };
    auto tryAssignJob = [&](ID j, const Pos &w) {
        if (isFreeWorker(w)) { assignJob(j, w); return true; }
        return false;
    };
    auto jobForWorker = [&](const Pos &w) { return jobForWorkers[w.x][w.y]; };
    auto swapJobs = [&](ID j0, ID j1) {
        Pos newWorker(workerOfJobs[j1]);
        assignJob(j1, workerOfJobs[j0]);
        assignJob(j0, newWorker);
    };

    ConsecutiveIdSet freeJobs(jobNum);
    for (ID j = 0; j < jobNum; ++j) { freeJobs.insert(j); }

    auto initByJob = [&]() {
        while (!freeJobs.empty()) {
            ID j = freeJobs.itemAt(r.pick(freeJobs.size()));

            if (!tryAssignJob(j, jobPosHints[j])) {
                SpiralScanner ss(jobPosHints[j]);
                ss.scan([&](const Pos &w) {
                    if (!isWorker(w)) { return false; }
                    return tryAssignJob(j, w);
                });
            }

            freeJobs.eraseItem(j);
        }
    };
    //initByJob();

    auto initByWorker = [&]() {
        for (Pos w = { 0, 0 }; w.x < workerMat.x; ++w.x) {
            for (w.y = 0; w.y < workerMat.y; ++w.y) {
                ID bestJob = InvalidId;
                Int bestDist = SafeLimit<Int>::max;
                ID maxStep = max(freeJobs.size() / (1 << 11), 1);
                for (ID j = 0; j < freeJobs.size();) {
                    ID job = freeJobs.itemAt(j);
                    Int dist = distance(w, jobPosHints[job]);
                    if (dist <= 0) { bestJob = job; break; }
                    if (math::updateMin(bestDist, dist)) { bestJob = job; }
                    (++j) += r.pick(maxStep);
                }
                if (bestJob == InvalidId) { break; }
                assignJob(bestJob, w);
                freeJobs.eraseItem(bestJob);
            }

#ifdef LOCAL
			reprint(to_string(w.x) + "    " + to_string(freeJobs.size()) + "  ");
#endif // LOCAL
        }
        cerr << endl;
    };
    //initByWorker();

	Func<void(ID, ID, ID, ID, Vec<ID>)> initByRecursiveTwoPartition = [&](ID left, ID bottom, ID right, ID top, Vec<ID> jobs) {
		//db2(right - left , top - bottom);
		enum { Left, Right, CornerNum };

		ID partJobNum = static_cast<ID>(jobs.size());
		if (partJobNum == 1) {
			assignJob(jobs[0], { left, bottom });
			return;
		}
		else if (jobs.empty()) {
			return;
		}

		int minX = inf, minY = inf, minZ = inf;
		int maxX = -inf, maxY = -inf, maxZ = inf;
		for (auto &job : jobs) {
			math::updateMin(minX, jobPosHints[job].x);
			math::updateMax(maxX, jobPosHints[job].x);
			math::updateMin(minY, jobPosHints[job].y);
			math::updateMax(maxY, jobPosHints[job].y);
			math::updateMin(minZ, jobPosHints[job].z);
			math::updateMax(maxZ, jobPosHints[job].z);
		}

		Arr<Pos, CornerNum> corners = { {
			{ minX, minY, minZ },
			{ maxX, minY, maxZ },
		} }; // TODO[szx][9]: { subMat.x - 1, subMat.y - 1 }

		Pos leftBottomPos, rightTopPos;
		if (right - left > top - bottom) {
			/*		 mid
			|---------------------|
			|					  |
			|---------------------|
			*/
			double leftRatio = 0.5;
			int midx = left + int((right - left) * leftRatio);
			math::updateMax(midx, left + 1);
			math::updateMin(midx, right - 1);
			leftBottomPos = { midx, bottom };
			rightTopPos = { midx, top };
			corners[Left] = { minX, (minY + maxY) / 2, 0};
			corners[Right] = { maxX, (minY + maxY) / 2, 0 };
		}
		else {
			double leftRatio = 0.5;
			int midy = bottom + int((top - bottom) * leftRatio);
			math::updateMax(midy, bottom + 1);
			math::updateMin(midy, top - 1);
			leftBottomPos = { left, midy };
			rightTopPos = { right, midy };
			corners[Left] = { (minX + maxX) / 2, minY, 0 };
			corners[Right] = { (minX + maxX) / 2, maxY, 0 };
		}

		Arr<ID, CornerNum> cornerWorkerNums;
		cornerWorkerNums[Left] = (rightTopPos.x - left) * (rightTopPos.y - bottom);
		cornerWorkerNums[Right] = (right - leftBottomPos.x) * (top - leftBottomPos.y);

		//db2(cornerWorkerNums[Left], cornerWorkerNums[Right]);
		//db2(cornerWorkerNums[Left] + cornerWorkerNums[Right], partJobNum);

		Vec<Int> distanceToCorners(partJobNum);
		Arr<Vec<ID>, CornerNum> jobRanks;
		for (ID c = 0; c < CornerNum; ++c) {
			if (cornerWorkerNums[c] <= 0) { continue; }
			for (ID j = 0; j < partJobNum; ++j) {
				distanceToCorners[j] = distance(jobPosHints[jobs[j]], corners[c]);
			}
			jobRanks[c].resize(partJobNum); // OPTIMIZE[szx][9]: shuffle?
			for (ID j = 0; j < partJobNum; ++j) { jobRanks[c][j] = j; }
			sort(jobRanks[c].begin(), jobRanks[c].end(), [&](ID l, ID r) { return distanceToCorners[l] > distanceToCorners[r]; });
		}

		Vec<bool> isFreeJobs(partJobNum, true);
		Arr<Vec<ID>, CornerNum> jobAtCorners;
		for (ID c = 0; c < CornerNum; ++c) { jobAtCorners[c].reserve((partJobNum / CornerNum) + 1); }
		while (partJobNum > 0) {
			for (ID c = 0; c < CornerNum; ++c) {
				if (jobAtCorners[c].size() >= cornerWorkerNums[c]) { continue; }
				while (partJobNum > 0) {
					ID j = jobRanks[c].back();
					jobRanks[c].pop_back();
					if (!isFreeJobs[j]) { continue; }
					jobAtCorners[c].push_back(jobs[j]);
					isFreeJobs[j] = false;
					--partJobNum;
					break;
				}
			}
		}

		fatalif(cornerWorkerNums[Left] != jobAtCorners[Left].size() || 
				cornerWorkerNums[Right] != jobAtCorners[Right].size(), "ERROR on assign!");

		initByRecursiveTwoPartition(left, bottom, rightTopPos.x, rightTopPos.y, jobAtCorners[Left]);
		initByRecursiveTwoPartition(leftBottomPos.x, leftBottomPos.y, right, top, jobAtCorners[Right]);
	};

    Func<void(ID, ID, ID, ID, Vec<ID>)> initByRecursivePartition = [&](ID left, ID bottom, ID right, ID top, Vec<ID> jobs) {
        enum { LeftBottom, RightBottom, LeftTop, RightTop, CornerNum };

        ID partJobNum = static_cast<ID>(jobs.size());
        if (partJobNum == 1) {
            assignJob(jobs[0], { left, bottom });
            return;
        } else if (jobs.empty()) {
            return;
        }

		int minX = inf, minY = inf, minZ = inf;
		int maxX = -inf, maxY = -inf, maxZ = inf;
		for (auto &job : jobs) {
			math::updateMin(minX, jobPosHints[job].x);
			math::updateMax(maxX, jobPosHints[job].x);
			math::updateMin(minY, jobPosHints[job].y);
			math::updateMax(maxY, jobPosHints[job].y);
			math::updateMin(minZ, jobPosHints[job].z);
			math::updateMax(maxZ, jobPosHints[job].z);
		}

		Pos midp = { (minX + maxX) / 2,(minY + maxY) / 2,(minZ + maxZ) / 2 };
		Arr<Pos, CornerNum> corners = { {
			{ minX, minY, minZ },	//论文结果
			{ maxX, minY, maxZ },
			{ minX, maxY, minZ },
			{ maxX, maxY, maxZ },
			//{ minX, minY, minZ },	//最远四个
			//{ maxX, minY, maxZ },
			//{ minX, maxY, maxZ },
			//{ maxX, maxY, minZ },
			//{ minX, minY, minZ },	//底面四个
			//{ maxX, minY, minZ },
			//{ minX, maxY, minZ },
			//{ maxX, maxY, minZ },
		} }; // TODO[szx][9]: { subMat.x - 1, subMat.y - 1 }

        Pos mid = { (left + right) / 2, (bottom + top) / 2 };
        Arr<ID, CornerNum> cornerWorkerNums;
        cornerWorkerNums[LeftBottom] = (mid.x - left) * (mid.y - bottom);
        cornerWorkerNums[RightBottom] = (right - mid.x) * (mid.y - bottom);
        cornerWorkerNums[LeftTop] = (mid.x - left) * (top - mid.y);
        cornerWorkerNums[RightTop] = (right - mid.x) * (top - mid.y);

        Vec<Int> distanceToCorners(partJobNum);
		Vec<Int> distanceToOCorners(partJobNum);
        Arr<Vec<ID>, CornerNum> jobRanks;
        for (ID c = 0; c < CornerNum; ++c) {
            if (cornerWorkerNums[c] <= 0) { continue; }
            for (ID j = 0; j < partJobNum; ++j) { 
				distanceToCorners[j] = distance(jobPosHints[jobs[j]], corners[c]) 
									 + abs(jobPosHints[jobs[j]].z - corners[c].z)
					; 
				distanceToOCorners[j] = inf;
				for (ID _c = 0; _c < CornerNum; ++_c) {
					if (_c == c) continue;
					math::updateMin(distanceToOCorners[j], 
						distance(jobPosHints[jobs[j]], corners[_c])
						+ abs(jobPosHints[jobs[j]].z - corners[_c].z)
					);
				}
			}

            jobRanks[c].resize(partJobNum); // OPTIMIZE[szx][9]: shuffle?
            for (ID j = 0; j < partJobNum; ++j) { jobRanks[c][j] = j; }
            sort(jobRanks[c].begin(), jobRanks[c].end(), [&](ID l, ID r) {
				Int dl = distanceToCorners[l];
				Int dr = distanceToCorners[r];
				return dl > dr || dl == dr && distanceToOCorners[l] < distanceToOCorners[r];
			});
        }

        Vec<bool> isFreeJobs(partJobNum, true);
        Arr<Vec<ID>, CornerNum> jobAtCorners;
        for (ID c = 0; c < CornerNum; ++c) { jobAtCorners[c].reserve((partJobNum / CornerNum) + 1); }
        while (partJobNum > 0) {
            for (ID c = 0; c < CornerNum; ++c) {
                if (jobAtCorners[c].size() >= cornerWorkerNums[c]) { continue; }
                while (partJobNum > 0) {
                    ID j = jobRanks[c].back();
                    jobRanks[c].pop_back();
                    if (!isFreeJobs[j]) { continue; }
                    jobAtCorners[c].push_back(jobs[j]);
                    isFreeJobs[j] = false;
                    --partJobNum;
                    break;
                }

            }
        }

        initByRecursivePartition(left, bottom, mid.x, mid.y, jobAtCorners[LeftBottom]);
        initByRecursivePartition(mid.x, bottom, right, mid.y, jobAtCorners[RightBottom]);
        initByRecursivePartition(left, mid.y, mid.x, top, jobAtCorners[LeftTop]);
        initByRecursivePartition(mid.x, mid.y, right, top, jobAtCorners[RightTop]);
    };

    Vec<ID> jobs(jobNum);
    for (ID j = 0; j < jobNum; ++j) { jobs[j] = j; }
    initByRecursivePartition(0, 0, workerMat.x, workerMat.y, jobs);
	//initByRecursiveTwoPartition(0, 0, workerMat.x, workerMat.y, jobs);

    auto distanceSum = [&](ID job) {
        Int d = 0;
        for (auto adjJob = jobAdjList[job].begin(); adjJob != jobAdjList[job].end(); ++adjJob) {
            d += normDistance(workerOfJobs[job], workerOfJobs[*adjJob]);
        }
        return d;
    };
    auto distanceNewSum = [&](ID job, const Pos &worker) {
        Int d = 0;
        for (auto adjJob = jobAdjList[job].begin(); adjJob != jobAdjList[job].end(); ++adjJob) {
            const Pos &workerOfAdjJob((worker != workerOfJobs[*adjJob]) ? workerOfJobs[*adjJob] : workerOfJobs[job]);
            d += normDistance(worker, workerOfAdjJob);
        }
        return d;
    };

    auto localSearch = [&](ID radius) {
        ID doubleR = 2 * radius;
        ID diameter = doubleR + 1;
        
        Vec<Int> jobDistances(jobNum);
        auto refreshDistance = [&](ID job) { jobDistances[job] = distanceSum(job); };
        auto refreshDistances = [&](ID job) {
            refreshDistance(job);
            for (auto adjJob = jobAdjList[job].begin(); adjJob != jobAdjList[job].end(); ++adjJob) {
                refreshDistance(*adjJob); // OPTIMIZE[szx][3]: only the wire which connects `job` is changed.
            } // OPTIMIZE[szx][3]: (jobDistances[*adjJob] -= normDistance(oldWorker, workerOfJobs[*adjJob])) += normDistance(newWorker, workerOfJobs[*adjJob])
        };
        for (ID j = 0; j < jobNum; ++j) { refreshDistance(j); }

        Vec<Arr2D<Int, ID>> moveDelta(jobNum, Arr2D<Int, ID>(diameter, diameter));
        auto refreshMoveDeltas = [&](ID job) {
            Pos newWorker = { workerOfJobs[job].x - radius };
            for (ID x = 0; x < diameter; ++x, ++newWorker.x) { // OPTIMIZE[szx][1]: dx is the same when moving along y-axis.
                newWorker.y = workerOfJobs[job].y - radius;
                for (ID y = 0; y < diameter; ++y, ++newWorker.y) {
                    if (!isWorker(newWorker)) { moveDelta[job][x][y] = SafeLimit<Int>::max; continue; }
                    if (newWorker == workerOfJobs[job]) { moveDelta[job][x][y] = SafeLimit<Int>::max; continue; }
                    moveDelta[job][x][y] = distanceNewSum(job, newWorker) - jobDistances[job];
                }
            }
        };
        for (ID j = 0; j < jobNum; ++j) { refreshMoveDeltas(j); }

        auto refreshSwapDeltas = [&](ID job) {
            refreshMoveDeltas(job);
            for (auto adjJob = jobAdjList[job].begin(); adjJob != jobAdjList[job].end(); ++adjJob) {
                refreshMoveDeltas(*adjJob);
            }
        };

        struct Move {
            Int objDelta;
            ID jobToSwap;
        };

        Int obj = 0;
        for (auto j = jobDistances.begin(); j != jobDistances.end(); ++j) { obj += *j; }
        obj /= 2;
        #ifdef LOCAL
        cerr << obj << " (InitObj)" << endl;
        #endif

        Timer timer(25min);
        Sampling1 s(r);
        for (ID iteration = 0; !timer.isTimeOut(); ++iteration) {
            Move bestMove = { SafeLimit<Int>::max };

            ID j = r.pick(jobNum);

            // find the best move.
            Pos newWorker = { workerOfJobs[j].x - radius };
            for (ID x = 0; x < diameter; ++x, ++newWorker.x) { // OPTIMIZE[szx][1]: dx is the same when moving along y-axis.
                if (newWorker.x < 0) { continue; }
                if (newWorker.x >= workerMat.x) { break; }
                newWorker.y = workerOfJobs[j].y - radius;
                for (ID y = 0; y < diameter; ++y, ++newWorker.y) {
                    if (newWorker.y < 0) { continue; }
                    if (newWorker.y >= workerMat.y) { break; }
                    ID jobToSwap = jobForWorker(newWorker);
                    Int delta = moveDelta[j][x][y] + moveDelta[jobToSwap][doubleR - x][doubleR - y];
                    if (s.isMinimal(bestMove.objDelta, delta)) { bestMove = { delta, jobToSwap }; }
                }
            }

            // make the move.
            swapJobs(j, bestMove.jobToSwap);
            obj += bestMove.objDelta;

            // refresh cache.
            refreshDistances(j);
            refreshDistances(bestMove.jobToSwap);
            refreshSwapDeltas(j); // refresh swap delta after all distances have been refreshed since there are dependencies.
            refreshSwapDeltas(bestMove.jobToSwap);

            #ifdef LOCAL
            if (iteration % (1 << 12) == 0) { reprint(to_string(obj) + "    " + to_string(bestMove.objDelta) + "  "); }
            #endif
            // check objective value.
            //Int cobj = 0;
            //for (auto j = jobDistances.begin(); j != jobDistances.end(); ++j) { cobj += *j; }
            //cobj /= 2;
            //if (cobj != obj) {
            //    db(cobj - obj);
            //}
        }
    };
    //localSearch(3);

    return workerOfJobs;
}


Vec<Pos> solveQapByLocalSearch(Vec<Vec<ID>> &jobAdjList, const Pos &workerMat, Vec<Pos> &jobPosHints) {
    Int maxDistance = workerMat.x + workerMat.y + 2;
    Vec<Int> normDistances(maxDistance);
    for (Int i = 0; i < maxDistance; ++i) {
        normDistances[i] = Int(i * sqrt(i) * DistanceAmp);
    }
    auto distance = [](const Pos &w0, const Pos &w1) {
        return abs(w0.x - w1.x) + abs(w0.y - w1.y);
    };
    auto normDistance = [&](const Pos &w0, const Pos &w1) {
        return normDistances[distance(w0, w1)];
    };

    return solveQapByLocalSearch(jobAdjList, workerMat, normDistance, jobPosHints, distance);
}

enum { XY, XZ, YZ };
enum {
	Resolution,
	OriginX,
	OriginY,
	OriginZ,
};

Pos projection(ID x, ID y, ID z, int dir) {
	switch (dir){
		case XY:
			return Pos(x, y, z);
		case XZ:
			return Pos(x, z, y);
		case YZ:
			return Pos(y, z, x);
		default:
			fatalif(1, "The projection direction should be XY or XZ or YZ");
	}
}

int maxAreaProjectonDirection(const std::vector<std::vector<int>> &prisms_info) {
	int direction = XY;
	double maxArea = 0;
	vector<pair<int, int>>pts;

	ID max_z = -1;
	for (int dir = XY; dir <= YZ; ++dir) {
		ID mz = -1;
		for (size_t p = 0; p < prisms_info.size(); ++p) { // OPT[szx][0]: project the largest facet to the fabric.
			if (prisms_info[p][Resolution] == -1) { continue; }
			const ID &x = prisms_info[p][OriginX];
			const ID &y = prisms_info[p][OriginY];
			const ID &z = prisms_info[p][OriginZ];
			Pos pos = projection(x, y, z, dir);
			pts.push_back({ pos.x, pos.y });

			mz = max(mz, pos.z);
		}
		double area = areaOfConvexhull(pts);
		
		if (maxArea < area) {
			maxArea = area;
			direction = dir;

			max_z = mz;
		}
		pts.clear();

#ifdef LOCAL
		//db2(dir, area);
#endif // LOCAL
	}

	db2(direction, max_z);

	return direction;
}

Vec<Pos> calcJobPosHints(int dim,
    const std::vector<std::vector<int>> &prisms_info,
    int tile_cnt, int adapter_cnt,
    const std::vector<std::vector<int>> &connections) {
    
	auto posDiscretization = [&](Vec<Pos> &jobPos) {
		vector<int>xs, ys, zs;
		for (auto &p : jobPos) {
			xs.emplace_back(p.x);
			ys.emplace_back(p.y);
			zs.emplace_back(p.z);
		}
		sort(xs.begin(), xs.end()); xs.erase(unique(xs.begin(), xs.end()), xs.end());
		sort(ys.begin(), ys.end()); ys.erase(unique(ys.begin(), ys.end()), ys.end());
		sort(zs.begin(), zs.end()); zs.erase(unique(zs.begin(), zs.end()), zs.end());

		for (auto &p : jobPos) {
			p.x = lower_bound(xs.begin(), xs.end(), p.x) - xs.begin();
			p.y = lower_bound(ys.begin(), ys.end(), p.y) - ys.begin();
			p.z = lower_bound(zs.begin(), zs.end(), p.z) - zs.begin();
		}
	};

	auto queryJobPosRange = [&](const Vec<Pos> &jobPos, bool print, ID &xl, ID &xh, ID &yl, ID &yh, ID &zl, ID &zh) {
		xl = inf, xh = -inf;
		yl = inf, yh = -inf;
		zl = inf, zh = -inf;
		for (const auto &job : jobPos) {
			math::updateMin(xl, job.x);
			math::updateMax(xh, job.x);
			math::updateMin(yl, job.y);
			math::updateMax(yh, job.y);
			math::updateMin(zl, job.z);
			math::updateMax(zh, job.z);
		}

		if(print)db3(xh - xl, yh - yl, zh - zl);
	};

	/*
		bullet -> XZ or YZ
		propellerTip -> XZ
	*/
	int projectionDir = maxAreaProjectonDirection(prisms_info);

	db(projectionDir);

    int nodeNum = tile_cnt + adapter_cnt;

    Vec<int> idMap(prisms_info.size(), -1); // `idMap[t]` means the `t`_th non-empty tile is which prism.

    //ofstream ofsp("1.prism.nb");
    //ofsp << "ListPointPlot3D[{";

    Vec<Pos> jobPosHints(nodeNum);
    int tileNum = 0;
    for (size_t p = 0; p < prisms_info.size(); ++p) { // OPT[szx][0]: project the largest facet to the fabric.
        if (prisms_info[p][Resolution] == -1) { continue; }
		ID radius = 0;// (1 << prisms_info[p][Resolution]) * 5;
		const ID &x = prisms_info[p][OriginX] + radius;
		const ID &y = prisms_info[p][OriginY] + radius;
		const ID &z = prisms_info[p][OriginZ] + radius;
		jobPosHints[tileNum] = projection(x, y, z, projectionDir);

        idMap[tileNum++] = int(p);
        //ofsp << "{" << prisms_info[p][OriginX] << "," << prisms_info[p][OriginY] << "," << prisms_info[p][OriginZ] << "}," << endl;
    }

	ID minX = inf, maxX = -inf;
	ID minY = inf, maxY = -inf;
	ID minZ = inf, maxZ = -inf;

	queryJobPosRange(jobPosHints, true, minX, maxX, minY, maxY, minZ, maxZ);

	posDiscretization(jobPosHints);
	queryJobPosRange(jobPosHints, true, minX, maxX, minY, maxY, minZ, maxZ);

	//Rotation
	Pos mid = { (minX + maxX) / 2, (minY + maxY) / 2, 0 };
	double rad = 0;
	//rad = acos(-1) / 2;
	for (int i = 0; i < tileNum; ++i) {
		Pos tmp = jobPosHints[i];
		jobPosHints[i].x = (tmp.x - mid.x)*cos(-rad) - (tmp.y - mid.y)*sin(-rad) + mid.x;
		jobPosHints[i].y = (tmp.x - mid.x)*sin(-rad) + (tmp.y - mid.y)*cos(-rad) + mid.y;
	}
	db3(minX, minY, minZ);
	db3(maxX, maxY, maxZ);


	//Translation
	queryJobPosRange(jobPosHints, true, minX, maxX, minY, maxY, minZ, maxZ);

	for (int i = 0; i < tileNum; ++i) {
		jobPosHints[i].x -= minX;
		jobPosHints[i].y -= minY;
		jobPosHints[i].z -= minZ;
	}
	maxX -= minX; maxY -= minY; maxZ -= minZ;
	minX = minY = minZ = 0;
	db3(maxX - minX, maxY - minY, maxZ - minZ);

    //ofsp.seekp(-1, ios::end);
    //ofsp << "}]";
    //ofsp.close();

    for (int a = tile_cnt; a < nodeNum; ++a) {
        jobPosHints[a] = { 0, 0, 0 };
        // calcualte the coord of the adapter.
        std::set<int> diffX;
        std::set<int> diffY;
		std::set<int> diffZ;
        for (auto t = connections[a].begin(); t != connections[a].end(); ++t) {
            if (diffX.find(jobPosHints[*t].x) == diffX.end()) {
                jobPosHints[a].x += jobPosHints[*t].x;
                diffX.insert(jobPosHints[*t].x);
            }
            if (diffY.find(jobPosHints[*t].y) == diffY.end()) {
                jobPosHints[a].y += jobPosHints[*t].y;
                diffY.insert(jobPosHints[*t].y);
            }
			if (diffZ.find(jobPosHints[*t].z) == diffZ.end()) {
				jobPosHints[a].z += jobPosHints[*t].z;
				diffZ.insert(jobPosHints[*t].z);
			}
        }
		jobPosHints[a].x /= (ID)diffX.size();
        jobPosHints[a].y /= (ID)diffY.size();
		jobPosHints[a].z /= (ID)diffZ.size();
    }

	int toX = 1, toY = 1;

	double ratioY = 1;
	//ratioY = maxX * 1060. / 800 / maxY;

	db(ratioY);

	for (ID i = 0; i < nodeNum; ++i) {
		jobPosHints[i].x = (jobPosHints[i].x * maxZ + jobPosHints[i].z * toX);
		jobPosHints[i].y = (jobPosHints[i].y * maxZ + jobPosHints[i].z * toY) * ratioY;
	}

	posDiscretization(jobPosHints);
	queryJobPosRange(jobPosHints, true, minX, maxX, minY, maxY, minZ, maxZ);

	db3(maxX - minX, maxY - minY, 1.0 * (maxY - minY) / (maxX - minX));

    return jobPosHints;
}


void testLocalSearch() {
    //   z
    //   ^
    //   |      4
    //   |    / |
    //   5  1 - 3
    //   | /   /
    //   0 - 2 ----> x
    Vec<Vec<ID>> jobAdjList({
        { 1, 2, 5 },
        { 0, 3, 5 },
        { 0, 3, 4 },
        { 1, 2 },
        { 2 },
        { 0, 1 },
    });
    Vec<Pos> jobPosHints({
        { 0, 0 },
        { 0, 1 },
        { 1, 0 },
        { 1, 1 },
        { 1, 1 },
        { 0, 0 },
    });

    Pos workerMat = { 2, 3 };

    solveQapByLocalSearch(jobAdjList, workerMat, jobPosHints);
}
