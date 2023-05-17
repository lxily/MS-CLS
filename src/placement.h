namespace Placement {
	#define MAXPECNTS 2000000
	#define ll long long
	#define inf 0x3f3f3f3f
	#define pii pair<int,int>
	#define mkp(a,b) make_pair(a,b)
	#define dequal(x, y) (abs(x - y) < DOUBLE_EPSILON)

	#define db(x) cerr<<#x<<" = "<<(x)<<endl
	#define db2(x,y) cerr<<#x<<" = "<<(x)<<", "<<#y<<" = "<<(y)<<endl
	#define db3(x,y,z) cerr<<#x<<" = "<<(x)<<", "<<#y<<" = "<<(y)<<", "<<#z<<" = "<<(z)<<endl
	
	/***************Problem parameter*****************/

	vector<vector<int>> connections;    //connect graph
	vector<vector<int>> prisms_info;    //prisms info
	vector<int>cubeShape;
	int inverse_sampling_step;
	int tile_cnt = 0, adapter_cnt = 0;        //
	int fabric_XSize = 0, fabric_YSize = 0;  //size of fabric

	ofstream saveIterationData;

	/*************************************************/

	/**************For Initial solution***************/

	bool is_used_tile[MAXPECNTS];
	bool is_used_place[MAXPECNTS];
	vector<double>value1D5;
	vector<vector<int>>fabric;
	const double theoreticalOpt[7] = {
		0., 1., 2., 3., 4., 6.8284271247, 9.6568542495
	};
	const int dir[4][2] = { {1, 0}, {-1, 0}, {0, 1}, {0, -1} };

	inline int get_idx_of_fabric(int x, int y) {
		return fabric[x][y];
	}

	struct Point {
		int x, y;
		Point(int _x = 0, int _y = 0) :x(_x), y(_y) {}
		bool operator == (const Point &p) const {
			return x == p.x && y == p.y;
		}
		bool operator != (const Point &p) const {
			return x != p.x || y != p.y;
		}

		inline int dist(const Point &p) {
			return abs(x - p.x) + abs(y - p.y);
		}

	};
	bool operator < (const Point &A, const Point &B) {
		return A.x - B.x ? A.x < B.x : A.y < B.y;
	}

	inline Point get_point_by_idx(int idx) {
		return Point(idx % fabric_XSize, idx / fabric_XSize);
	}

	struct Buckets {
		vector<vector<int>>_vals;
		vector<pair<int, int>>pos;

		Buckets() {
			_vals = {};
			pos = {};
		}

		void clear() {
			vector<vector<int>>().swap(_vals);
			vector<pair<int, int>>().swap(pos);
		}

		void insert(int b, int v) {
			if (_vals.size() <= b) {
				_vals.resize(b + 1);
			}
			if (pos.size() <= v) {
				pos.resize(v + 1);
			}
			pos[v] = make_pair(b, (int)_vals[b].size());
			_vals[b].push_back(v);
		}

		void remove(int v) {
			pair<int, int>p = pos[v];

			vector<int> &V = _vals[p.first];

			int sz = (int)V.size();
			if (p.second != sz - 1) {
				V[p.second] = V[sz - 1];
				pos[V[sz - 1]] = p;
			}
			V.pop_back();
		}
	};

	struct Solution {
		const static int alpha = 1;
		const static int beta = 0;

		double score;

		vector<int>tile_in_place;
		vector<Point>place_of_tile;
		vector<double>realTimeCost;
		vector<int>tabuList[2];

		Buckets costBucket;

		void initCostBucket() {
			double avgCost = score / connections.size();
			for (int i = 0; i < (int)connections.size(); ++i) {
				costBucket.insert(realTimeCost[i] >= avgCost, i);
			}
		}

		void insertToBucket(int tile_id) {
			double avgCost = score / (int)connections.size();
			costBucket.insert(
				realTimeCost[tile_id] >= avgCost * 2 ? 3 :
				realTimeCost[tile_id] >= avgCost ? 2 :
				realTimeCost[tile_id] >= avgCost / 2 ? 1 : 0, tile_id);
		}

		void removeFromBucket(int tile_id) {
			costBucket.remove(tile_id);
		}

		vector<int>sampleTilesFromBucket(vector<int>sampleSize) {
			vector<int>ret;
			unordered_map<int, bool>vis;
			for (int i = 0; i < 2; ++i) {
				int sz = (int)costBucket._vals[i].size();
				for (int j = 0; j < sampleSize[i]; ++j) {
					int id = costBucket._vals[i][rand() % sz];
					if (vis.count(id))continue;
					vis[id] = true;
					ret.push_back(id);
				}
			}
			return ret;
		}

		Solution() {
			score = inf;
			//costBucket.clear();
			tile_in_place.resize(MAXPECNTS, -1);
			place_of_tile.resize(connections.size());
			realTimeCost.resize(connections.size(), 0.);
		}

		void initRealTimeCost() {
			for (int i = 0; i < (int)connections.size(); ++i) {
				for (int j = 0; j < (int)connections[i].size(); ++j) {
					int adj_id = connections[i][j];
					realTimeCost[i] += value1D5[place_of_tile[i].dist(place_of_tile[adj_id])];
				}
			}
		}

		double getProbability(int tile_id) {
			double avgCost = score / connections.size();
			return 1. / (1. + exp(- (realTimeCost[tile_id] - avgCost) / avgCost / 2));
		}

		double sumOfCost() {
			double ret = 0;
			for (int i = 0; i < (int)connections.size(); ++i) {
				ret += realTimeCost[i];
			}
			return ret;
		}

		void insert(int tile_id, Point plc) {
			int place_id = get_idx_of_fabric(plc.x, plc.y);
			fatalif(is_used_tile[tile_id], "Tile %d have been placed!", tile_id);
			fatalif(is_used_place[place_id], "Position (%d, %d) already occupied!", plc.x, plc.y);

			tile_in_place[place_id] = tile_id;
			place_of_tile[tile_id] = plc;

			is_used_tile[tile_id] = true;
			is_used_place[place_id] = true;
		}

		void remove(int tile_id) {
			fatalif(!is_used_tile[tile_id], "Tile %d have not been placed!", tile_id);
			Point plc = place_of_tile[tile_id];
			int place_id = get_idx_of_fabric(plc.x, plc.y);

			tile_in_place[place_id] = -1;
			is_used_tile[tile_id] = false;
			is_used_place[place_id] = false;
		}

		inline Point get_place(int tile_id) {
			fatalif(!is_used_tile[tile_id], "Tile %d have not been placed!", tile_id);
			return place_of_tile[tile_id];
		}

		inline int get_tile_id(int x, int y) {
			return tile_in_place[get_idx_of_fabric(x, y)];
		}

		inline int get_tile_id(Point plc) {
			return tile_in_place[get_idx_of_fabric(plc.x, plc.y)];
		}

		inline int get_tile_id(int place_id) {
			return tile_in_place[place_id];
		}

		inline double getRealTimeCost(int tile_id) {
			return tile_id == -1 ? 0 : realTimeCost[tile_id];
		}

		void swap(int tile_id1, int tile_id2) {
			Point plc1 = place_of_tile[tile_id1];
			Point plc2 = place_of_tile[tile_id2];
			int place_id1 = get_idx_of_fabric(plc1.x, plc1.y);
			int place_id2 = get_idx_of_fabric(plc2.x, plc2.y);

			std::swap(place_of_tile[tile_id1], place_of_tile[tile_id2]);
			std::swap(tile_in_place[place_id1], tile_in_place[place_id2]);

			double new_cost_id1 = 0;
			for (int i = 0; i < (int)connections[tile_id1].size(); ++i) {
				int adj_id = connections[tile_id1][i];
				double newLength = value1D5[plc2.dist(get_place(adj_id))];
				new_cost_id1 += newLength;
				if (adj_id != tile_id1 && adj_id != tile_id2) {
					realTimeCost[adj_id] += newLength - value1D5[get_place(adj_id).dist(plc1)];
				}
			}
			//removeFromBucket(tile_id1);
			realTimeCost[tile_id1] = new_cost_id1;
			//insertToBucket(tile_id1);
			

			double new_cost_id2 = 0;
			for (int i = 0; i < (int)connections[tile_id2].size(); ++i) {
				int adj_id = connections[tile_id2][i];
				double newLength = value1D5[plc1.dist(get_place(adj_id))];
				new_cost_id2 += newLength;
				if (adj_id != tile_id1 && adj_id != tile_id2) {
					realTimeCost[adj_id] += newLength - value1D5[get_place(adj_id).dist(plc2)];
				}
			}
			//removeFromBucket(tile_id2);
			realTimeCost[tile_id2] = new_cost_id2;
			//insertToBucket(tile_id2);
		}

		void init_tabu_list() {
			tabuList[0].resize(tile_cnt + adapter_cnt, 0);
			tabuList[1].resize(tile_cnt + adapter_cnt, 0);
		}

		inline void set_tabu(int list_id, int tile_id, int iter) {//list_id == 1 -> worse
			if(tile_id >=0)tabuList[list_id][tile_id] = iter + rand() % alpha + beta;
		}

		inline bool is_tabu(int list_id, int tile_id, int iter) {
			return tile_id < 0? false: tabuList[list_id][tile_id] > iter;
		}

		void reset_tabu() {
			tabuList[0].clear();
			tabuList[1].clear();
			tabuList[0].resize(tile_cnt + adapter_cnt, 0);
			tabuList[1].resize(tile_cnt + adapter_cnt, 0);
		}
	};

	struct EdgeSpacePoint {
		vector<Point>pts;
		vector<bool>inque;
		vector<int>idxinpts;

		EdgeSpacePoint() {
			pts.clear();
			inque.resize(MAXPECNTS, false);
			idxinpts.resize(MAXPECNTS, -1);
		}

		void clear() {
			pts.clear();
		}

		void insert(Point _p) {
			if (is_inque(_p.x, _p.y))return;
			int pid = get_idx_of_fabric(_p.x, _p.y);
			idxinpts[pid] = (int)pts.size();
			pts.push_back(_p);
			inque[pid] = true;
		}

		void remove(Point _p) {

			int pid = get_idx_of_fabric(_p.x, _p.y);
			int pidx = idxinpts[pid];
			if (!inque[pid])return;

			if (pidx != (int)pts.size() - 1) {
				Point &lp = pts[pts.size() - 1];
				int lpid = get_idx_of_fabric(lp.x, lp.y);
				pts[pidx] = lp;
				idxinpts[lpid] = pidx;
			}

			pts.pop_back();
			idxinpts[pid] = -1;
			inque[pid] = false;
		}

		bool is_inque(int x, int y) {
			return inque[get_idx_of_fabric(x, y)];
		}

		Point get_min_cost_point(Solution &sol, int tile_id) {
			//vector<Point>best_points;

			double minCost = (double)inf;
			Point best_Point = Point(fabric_XSize, fabric_YSize);
			for (int i = 0; i < (int)pts.size(); ++i) {
				double cost = 0;
				Point p = pts[i];
				fatalif(is_used_place[get_idx_of_fabric(pts[i].x, pts[i].y)], "(%d, %d) Used!\n", pts[i].x, pts[i].y);
				for (int j = 0; j < (int)connections[tile_id].size(); ++j) {
					int adj_id = connections[tile_id][j];
					if (is_used_tile[adj_id]) {
						cost += value1D5[sol.get_place(adj_id).dist(p)];
					}
				}
				if (minCost > cost || (dequal(minCost, cost) && p < best_Point)) {
					minCost = cost;
					best_Point = p;
				}

				/*if (minCost > cost) {
					minCost = cost;
					best_points.clear();
					best_points.push_back(p);
				}
				else if (dequal(minCost, cost)) {
					best_points.push_back(p);
				}*/
			}

			return best_Point;
		}
	}espt;//最近的候选放置点集

	struct MyPrism {
		int resolution;
		int dim;
		int pre_tiles_cnt;
		int origin[3];
		int shape[3];

		MyPrism() {
		}

		void print_info() {
			if (dim == 2) {
				printf("Origin: [%d, %d]\n", origin[0], origin[1]);
				printf("Shape: [%d, %d]\n", shape[0], shape[1]);

			}
			else {
				printf("Origin: [%d, %d, %d]\n", origin[0], origin[1], origin[2]);
				printf("Shape: [%d, %d, %d]\n", shape[0], shape[1], shape[2]);
			}
		}

		bool operator < (const MyPrism &p) const {
			for (int i = dim - 1; i >= 0; --i) {
				if (origin[i] != p.origin[i]) {
					return origin[i] < p.origin[i];
				}
			}
			if (dim == 2) {
				fatalif(1, "There are two identical origin coordinates [%d, %d]!", origin[0], origin[1]);
			}
			else {
				fatalif(1, "There are two identical origin coordinates [%d, %d, %d]!", origin[0], origin[1], origin[2]);
			}

			return false;
		}

		inline int get_tiles()const {
			return dim == 2 ? shape[0] * shape[1] : shape[0] * shape[1] * shape[2];
		}

		inline int get_tile_id(int x, int y, int z) {
			if (dim == 2) {
				return pre_tiles_cnt + shape[0] * y + x;
			}
			else {
				return pre_tiles_cnt + shape[0] * shape[1] * z + shape[0] * y + x;
			}
		}
	};
	vector<MyPrism>prisms;

	void get_prisms() {

		int pre_tiles = 0;
		for (int i = 0; i < (int)prisms_info.size(); ++i) {
			if (prisms_info[i][0] < 0)continue;

			MyPrism p;
			p.resolution = prisms_info[i][0];
			p.dim = prisms_info[i].size() == 5 ? 2 : 3;
			p.pre_tiles_cnt = pre_tiles;
			for (int j = 0; j < p.dim; ++j) {
				p.origin[j] = prisms_info[i][j + 1];
			}
			for (int j = 0; j < p.dim; ++j) {
				p.shape[j] = prisms_info[i][j + p.dim + 1];
			}
			pre_tiles += p.get_tiles();
			prisms.push_back(p);
		}

		//sort(prisms.begin(), prisms.end());
	}

	void place_one_tile(const int &tile_id, const Point &plc, Solution &sol) {
		sol.insert(tile_id, plc);
		if (espt.is_inque(plc.x, plc.y)) {
			espt.remove(plc);
		}

		int x_min = max(0, plc.x - 2);
		int x_max = min(fabric_XSize, plc.x + 3);
		int y_min = max(0, plc.y - 2);
		int y_max = min(fabric_YSize, plc.y + 3);
		for (int x = x_min; x < x_max; ++x) {
			for (int y = y_min; y < y_max; ++y) {
				if (!espt.is_inque(x, y) && !is_used_place[get_idx_of_fabric(x, y)]) {
					espt.insert(Point(x, y));
				}
			}
		}
	}

	void clear_one_tile(const int &tile_id, Solution &sol) {
		espt.insert(sol.get_place(tile_id));
		sol.remove(tile_id);
	}

	vector<vector<int>>divide_graph() {
		vector<vector<int>>ret;

		vector<bool>vis(MAXPECNTS, false);
		for (int i = 0; i < (int)connections.size(); ++i) {
			if (vis[i])continue;
			queue<int>Q;
			Q.push(i);
			vis[i] = true;
			vector<int>sub_gra;
			while (!Q.empty()) {
				int u = Q.front();
				Q.pop();
				sub_gra.push_back(u);
				for (auto &v : connections[u]) {
					if (vis[v])continue;
					vis[v] = true;
					Q.push(v);
				}
			}
			ret.push_back(sub_gra);
		}
		return ret;
	}

	vector<pii> get_dist_info(int start_tile_id) {
		vector<pii>ret;
		queue<pii>Q;
		vector<bool>vis(MAXPECNTS, false);

		ret.push_back(mkp(0, start_tile_id));
		Q.push(mkp(0, start_tile_id));
		vis[start_tile_id] = true;

		while (!Q.empty()) {
			pii out = Q.front();
			Q.pop();

			int u = out.second;
			int d = out.first;

			for (int i = 0; i < (int)connections[u].size(); ++i) {
				int v = connections[u][i];
				if (vis[v])continue;
				vis[v] = true;
				Q.push(mkp(d + 1, v));
				ret.push_back(mkp(d + 1, v));
			}
		}

		sort(ret.begin(), ret.end());

		return ret;
	}

	vector<Point>get_position(int pre_cnt, bool startWithCenter = false) {
		int use_fabX = (int)sqrt((int)connections.size());
		use_fabX = min(use_fabX, fabric_XSize - 1);

		int use_fabY = (int)ceil((int)connections.size() / use_fabX);
		if (use_fabX * use_fabY < (int)connections.size()) {
			if (use_fabX < fabric_XSize - 1) ++use_fabX;
			if (use_fabY < fabric_YSize - 1) ++use_fabY;
		}

		use_fabY = min(use_fabY, fabric_YSize - 1);

#ifdef LOCAL
		db2(use_fabX, use_fabY);
#endif //LOCAL
		vector<Point>ret;

		int x = (int)sqrt(pre_cnt);
		int y = (int)ceil(1.0 * pre_cnt / x);

		int x_delta = use_fabX / x;
		int y_delta = use_fabY / y;
		for (int i = 0; i < x; ++i) {
			for (int j = 0; j < y; ++j) {
				int px = i * x_delta;
				int py = j * y_delta;
				if (startWithCenter) {
					px += x_delta / 2;
					py += y_delta / 2;
				}
				ret.push_back(Point(px, py));
			}
		}

		sort(ret.begin(), ret.end(), [](Point &p1, Point &p2) {
			int d1 = p1.x + p1.y;
			int d2 = p2.x + p2.y;
			return d1 - d2 ? d1 < d2 : p1 < p2;
		});

		return ret;
	}

	vector<pair<int, Point>> get_tile_place_pair(int pre_cnt = 1, bool startWithCenter = false) {
		vector<pii>dist_info = get_dist_info(prisms[0].pre_tiles_cnt);

		vector<int>preTiles;
		int tileIdx = prisms[0].pre_tiles_cnt;
		preTiles.push_back(tileIdx);
		int delta = pre_cnt - 1 ? (int)dist_info.size() / (pre_cnt - 1) : inf;
		for (int i = 0; i < pre_cnt - 1; ++i) {
			int idx = min((int)dist_info.size() - 1, tileIdx + delta);
			if (idx != tileIdx) {
				preTiles.push_back(dist_info[tileIdx = idx].second);
			}
		}

		vector<Point>prePlace = get_position(pre_cnt, startWithCenter);

		vector<pair<int, Point>>ret;


		pre_cnt = (int)min(prePlace.size(), preTiles.size());
		for (int i = 0; i < pre_cnt; ++i) {
			ret.push_back({ preTiles[i], prePlace[i] });

#ifdef LOCAL
			//printf("Tile %d place at Point (%d, %d).\n", preTiles[i], prePlace[i].x, prePlace[i].y);
#endif // LOCAL
		}

		return ret;
	}

	Solution get_init_solution(const vector<pair<int, int>> &init_place, int pre_placement_cnt = 1, bool startWithCenter = false, string initialAgorithm = "projection") {
		
		
#ifdef LOCAL
		/*vector<vector<int>>graph_info = divide_graph();
		printf("Graph Infomation: ");
		for (auto &V : graph_info) { printf("%d ", (int)V.size()); }putchar('\n');*/
#endif // LOCAL

		queue<int>tileWaitingQue;
		Solution initSol;

		vector<pair<int, Point>>pre_deal;

		if (initialAgorithm == "projection") {
			for (int i = 0; i < (int)connections.size(); ++i) {
				place_one_tile(i, Point(init_place[i].first, init_place[i].second), initSol);
			}
		}
		else if (initialAgorithm == "greedy") {
			pre_deal = get_tile_place_pair(pre_placement_cnt, startWithCenter);
			for (int i = 0; i < (int)pre_deal.size(); ++i) {
				tileWaitingQue.push(pre_deal[i].first);
				place_one_tile(pre_deal[i].first, pre_deal[i].second, initSol);
			}
			while (!tileWaitingQue.empty()) {
				int tile_id = tileWaitingQue.front();
				tileWaitingQue.pop();

				for (int i = 0; i < (int)connections[tile_id].size(); ++i) {

					int nxt_tile_id = connections[tile_id][i];
					if (is_used_tile[nxt_tile_id])continue;

					Point pos = espt.get_min_cost_point(initSol, nxt_tile_id);

					place_one_tile(nxt_tile_id, pos, initSol);
					tileWaitingQue.push(nxt_tile_id);
				}
			}

			for (int i = 0; i < (int)connections.size(); ++i) {
				if (!is_used_tile[i]) {
					Point pos = espt.get_min_cost_point(initSol, i);
					place_one_tile(i, pos, initSol);
				}
			}
		}
		else if (initialAgorithm == "random") {
			int totalPEs = fabric_XSize * fabric_YSize;
			vector<int>order(totalPEs, 0);
			for (int i = 1; i < totalPEs; ++i) {
				order[i] = i;
				unsigned randPos = rand() * rand();
				randPos %= (i + 1);
				swap(order[randPos], order[i]);
			}
			for (unsigned i = 0; i < connections.size(); ++i) {
				int x = order[i] / fabric_YSize;
				int y = order[i] % fabric_YSize;
				place_one_tile(i, Point(x, y), initSol);
			}
		}
		else {
			fatalif("1", "%s not implemented!", initialAgorithm.c_str());
		}

#ifdef LOCAL
		printf("Init placement OK! (Use pre_deal = %d)\n", (int)pre_deal.size());
#endif // LOCAL
		return initSol;
	}

	/*************************************************/

	/***************For Save and Debug****************/

	void save_solution(char *save_road, const vector<Point> &psol) {
		ofstream outFile;
		outFile.open(save_road);

		outFile << "# inverse_sampling_step\n";
		outFile << inverse_sampling_step << '\n';

		outFile << "# Prisms\n";
		for (int i = 0; i < (int)prisms_info.size(); ++i) {
			for (int j = 0; j < (int)prisms_info[i].size(); ++j) {
				if (j)outFile << ' ';
				outFile << prisms_info[i][j];
			}
			outFile << '\n';
		}

		outFile << "compute_map:";
		for (int i = 0; i < tile_cnt; ++i) {
			outFile << '\n' << psol[i].x << ' ' << psol[i].y;
		}

		outFile << "\nadapter_map:";
		for (int i = tile_cnt; i < (int)psol.size(); ++i) {
			int reput = (int)connections[i].size() - 1;
			for (int j = 0; j < reput; ++j) {
				outFile << '\n' << psol[i].x << ' ' << psol[i].y;
			}
		}
		outFile.close();
	}

	void outputSolution(const vector<Point> &psol) {
		printf("# inverse_sampling_step\n");
		printf("%d\n", inverse_sampling_step);
		printf("# Prisms\n");

		for (int i = 0; i < (int)prisms_info.size(); ++i) {
			for (int j = 0; j < (int)prisms_info[i].size(); ++j) {
				if (j)putchar(' ');
				printf("%d", prisms_info[i][j]);
			}
			putchar('\n');
		}

		printf("compute_map:");
		for (int i = 0; i < tile_cnt; ++i) {
			printf("\n%d %d", psol[i].x, psol[i].y);
		}

		printf("\nadapter_map:");
		for (int i = tile_cnt; i < (int)psol.size(); ++i) {
			int reput = (int)connections[i].size() - 1;
			for (int j = 0; j < reput; ++j) {
				printf("\n%d %d", psol[i].x, psol[i].y);
			}
		}
	}

	void save_graph_and_fabric(char *save_road) {
		unordered_map<int, unordered_map<int, bool>>rec;

		ofstream outFile;
		outFile.open(save_road);

		int jobNum = (int)connections.size();
		int conCnt = 0;
		for (int i = 0; i < (int)connections.size(); ++i) {
			conCnt += (int)connections[i].size();
		}
		outFile << jobNum << ' ' << conCnt / 2 << '\n';
		for (int i = 0; i < jobNum; ++i) {
			for (int j = 0; j < (int)connections[i].size(); ++j) {
				if (rec[i][connections[i][j]] || rec[connections[i][j]][i])continue;
				outFile << i << ' ' << connections[i][j] << '\n';
				rec[i][connections[i][j]] = rec[connections[i][j]][i] = true;
			}
		}

		int workerNum = fabric_XSize * fabric_YSize;
		outFile << workerNum << '\n';
		for (int i = 0; i < workerNum; ++i) {
			for (int j = 0; j < workerNum; ++j) {
				Point p1 = get_point_by_idx(i);
				Point p2 = get_point_by_idx(j);
				double d = p1.dist(p2);
				outFile << d * sqrt(d) << ' ';
			}
			outFile << '\n';
		}

		outFile.close();
	}

	void print_connection_info() {
		printf("\nConnections Info: \n");
		for (int i = 0; i < (int)connections.size(); i += 1) {
			printf("%d -> : ", i);
			for (int j = 0; j < (int)connections[i].size(); j += 1) {
				printf("%2d ", (int)connections[i][j]);
			}
			printf("\n");
		}
	}

	void print_prisms_info() {
		printf("\nPrisms Info: \n");
		for (int i = 0; i < (int)prisms_info.size(); i += 1) {
			printf("Prism ID : %d, Shape: -> [", i);
			for (int j = 0; j < (int)prisms_info[i].size(); j += 1) {
				printf(j ? " ,%d" : "%d", prisms_info[i][j]);
			}
			printf("]\n");
		}
	}

	void print_placement(const vector<Point> &psol) {
		printf("compute_map:\n");
		for (int i = 0; i < tile_cnt; ++i) {
			printf("%d %d\n", psol[i].x, psol[i].y);
		}
		printf("adapter_map:\n");
		for (int i = tile_cnt; i < (int)psol.size(); ++i) {
			printf("%d %d\n", psol[i].x, psol[i].y);
		}
	}

	void analyze_solution(Solution &sol) {
		printf("\n=========================Analyze solution=========================\n\n");

		int optCnt = 0;
		double biasSum = 0;
		double biasAvg = 0;
		double degreeAvg = 0;
		double theoreticallyOpt = 0;

		for (int i = 0; i < (int)connections.size(); ++i) {
			double bias = fabs(theoreticalOpt[connections[i].size()] - sol.getRealTimeCost(i));
			optCnt += (bias < DOUBLE_EPSILON);
			biasSum += bias;
			degreeAvg += connections[i].size();
			theoreticallyOpt += theoreticalOpt[connections[i].size()];
		}
		biasAvg = biasSum / connections.size();
		degreeAvg /= connections.size();

		printf("Theoretically optimal: %.6lf, You got: %.6lf\n", pow(theoreticallyOpt, 2. / 3), pow(sol.score, 2. / 3));
		printf("Lower_bound: %.6lf, Normalized Connectivity: %.6lf\n", 
			pow(100. * fabric_XSize * fabric_YSize, 2. / 3), pow(100. * fabric_XSize * fabric_YSize, 2./3) / pow(sol.score, 2. / 3));
		printf("Total tiles = %d, Average degree = %.3lf\n", (int)connections.size(), degreeAvg);
		printf("OptCnt rate = %d / %d, biasSum = %.6lf\nbiasAvg (each tile) = %.6lf, biasAvg (each wire) = %.6lf\n",
			optCnt, (int)connections.size(), biasSum, biasAvg, biasAvg / degreeAvg);

		int lowerAvgCnt = 0;
		double avgCost = sol.score / connections.size();
		for (int i = 0; i < (int)connections.size(); ++i) {
			lowerAvgCnt += sol.getRealTimeCost(i) < avgCost;
		}
		printf("Here are %d below average, %d above average.\n", lowerAvgCnt, (int)connections.size() - lowerAvgCnt);

		//printf("Cost buckets infomation: "); for (int i = 0; i < (int)sol.costBucket._vals.size(); ++i)printf("%d ",
		//	(int)sol.costBucket._vals[i].size()); putchar('\n');

		if (connections.size() < 20) {
			for (int i = 0; i < (int)connections.size(); ++i) {
				printf("For tile %d, [Opt - %.8lf, Cur - %.8lf]\n", i, theoreticalOpt[connections[i].size()], sol.getRealTimeCost(i));
			}
		}
		printf("\n===========================Analyze End============================\n\n");
	}

	/*************************************************/

	/*****************For Local Search****************/

	double value_of_tile_at(int tile_id, Point plc, Solution &sol) {
		double res = 0;
		for (int i = 0; i < (int)connections[tile_id].size(); ++i) {
			int adj_tile_id = connections[tile_id][i];
			res += value1D5[plc.dist(sol.get_place(adj_tile_id))];
		}
		return res;
	}

	double evaluate(Solution &sol) {
		double score = 0;
		for (int i = 0; i < (int)connections.size(); ++i) {
			score += value_of_tile_at(i, sol.get_place(i), sol);
		}
		return score;
	}

	double trySwapTile(int tile_id1, int tile_id2, Point plc, Solution &sol) {
		double cost = 0;
		if (tile_id2 == -1) {
			for (int i = 0; i < (int)connections[tile_id1].size(); ++i) {
				cost += value1D5[plc.dist(sol.get_place(connections[tile_id1][i]))];
			}
		}
		else {
			swap(sol.place_of_tile[tile_id1], sol.place_of_tile[tile_id2]);
			cost += value_of_tile_at(tile_id1, sol.get_place(tile_id1), sol);
			cost += value_of_tile_at(tile_id2, sol.get_place(tile_id2), sol);
			swap(sol.place_of_tile[tile_id1], sol.place_of_tile[tile_id2]);
		}
		return cost;
	}

	void swapOperation(int tile_id1, int tile_id2, Point plc, Solution &sol) {
		if (tile_id2 == -1) {  //move to space
			Point prePlc = sol.get_place(tile_id1);
			clear_one_tile(tile_id1, sol);
			place_one_tile(tile_id1, plc, sol);
			double new_cost = 0;
			for (int i = 0; i < (int)connections[tile_id1].size(); ++i) {
				int adj_id = connections[tile_id1][i];
				double newLength = value1D5[plc.dist(sol.get_place(adj_id))];
				new_cost += newLength;
				sol.realTimeCost[adj_id] += newLength - value1D5[prePlc.dist(sol.get_place(adj_id))];
			}
			//sol.removeFromBucket(tile_id1);
			sol.realTimeCost[tile_id1] = new_cost;
			//sol.insertToBucket(tile_id1);
		}
		else {
			sol.swap(tile_id1, tile_id2);
		}
	}

	void moveAround(int tile_id, int searchRadiu, Solution &sol) {

		Point plc = sol.get_place(tile_id);

		double pre_value;
		double aft_value;
		int best_swap_id = -1;
		Point best_place = Point(-1, -1);
		double best_impro = (double)-inf;  //best improve

		int ctx = plc.x;
		int cty = plc.y;

		int use_radius = searchRadiu;
		int height = 2 * use_radius + 1;

		for (int p = 0; p < height; ++p) {

			int y = cty + use_radius - p;
			if (y < 0 || y >= fabric_YSize)continue;
			int x_left = max(0, p <= use_radius ? ctx - p : ctx - 2 * use_radius + p);
			int x_right = min(fabric_XSize - 1, p <= use_radius ? ctx + p : ctx + 2 * use_radius - p);

			for (int x = x_left; x <= x_right; ++x) {

				int swap_id = sol.get_tile_id(x, y);
				if (swap_id == tile_id)continue;

				pre_value = sol.getRealTimeCost(tile_id) + sol.getRealTimeCost(swap_id);
				aft_value = trySwapTile(tile_id, swap_id, Point(x, y), sol);

				if (aft_value < pre_value) {
					swapOperation(tile_id, swap_id, Point(x, y), sol);
					sol.score += 2 * (aft_value - pre_value);
				}

				if (best_impro < pre_value - aft_value) {
					best_swap_id = swap_id;
					best_place = Point(x, y);
					best_impro = pre_value - aft_value;
				}
			}
		}
		if (best_impro < 0 && (best_swap_id != -1 || best_place.x >= 0)) {
			swapOperation(tile_id, best_swap_id, best_place, sol);
			sol.score -= 2 * best_impro;
		}
	}

	void captureNeighbors(int tile_id, int searchRadiu, Solution &sol) {

		for (const auto &adj_id : connections[tile_id]) {
			double pre_value;
			double aft_value;
			int best_swap_id = -1;
			Point best_place = Point(-1, -1);
			double best_impro = (double)-inf;  //best improve

			Point plc = sol.get_place(adj_id);

			if (plc.dist(sol.get_place(tile_id)) == 1) {
				continue;
			}

			int ctx = plc.x;
			int cty = plc.y;
			int orgx = sol.get_place(tile_id).x;
			int orgy = sol.get_place(tile_id).y;

			int use_radius = min(searchRadiu, max(abs(ctx - orgx), abs(cty - orgy)) / 2 + 1);
			int height = 2 * use_radius + 1;

			int biasX = abs(orgx - ctx) * 1 / 2;
			int biasY = abs(orgy - cty) * 1 / 2;
			ctx += orgx < ctx ? -biasX : orgx > ctx ? biasX : 0;
			cty += orgy < cty ? -biasY : orgy > cty ? biasY : 0;
			ctx = max(0, ctx); ctx = min(fabric_XSize - 1, ctx);
			cty = max(0, cty); cty = min(fabric_YSize - 1, cty);

			for (int p = 0; p < height; ++p) {

				int y = cty + use_radius - p;
				if (y < 0 || y >= fabric_YSize)continue;
				int x_left = max(0, p <= use_radius ? ctx - p : ctx - 2 * use_radius + p);
				int x_right = min(fabric_XSize - 1, p <= use_radius ? ctx + p : ctx + 2 * use_radius - p);

				for (int x = x_left; x <= x_right; ++x) {

					int swap_id = sol.get_tile_id(x, y);
					if (swap_id == tile_id || swap_id == adj_id)continue;

					pre_value = sol.getRealTimeCost(adj_id) + sol.getRealTimeCost(swap_id);
					aft_value = trySwapTile(adj_id, swap_id, Point(x, y), sol);

					if (aft_value < pre_value) {
						swapOperation(adj_id, swap_id, Point(x, y), sol);
						sol.score += 2 * (aft_value - pre_value);
					}

					if (best_impro < pre_value - aft_value) {
						best_swap_id = swap_id;
						best_place = Point(x, y);
						best_impro = pre_value - aft_value;
					}
				}
			}
			if (best_impro < 0 && plc.dist(sol.get_place(tile_id)) > 2 && (best_swap_id != -1 || best_place.x >= 0)) {
				swapOperation(adj_id, best_swap_id, best_place, sol);
				sol.score -= 2 * best_impro;
			}
		}

	}

	void closeToNeighbors(int tile_id, int searchRadiu, Solution &sol) {
		bool improved_i = false;
		double best_impro = (double)-inf;  //最优改进
		int best_swap_id = -1;
		Point best_place = Point(-1, -1);
		double pre_value;
		double aft_value;
		for (int j = 0; j < (int)connections[tile_id].size(); ++j) {
			int adj_id = connections[tile_id][j];
			Point plc = sol.get_place(adj_id);

			int ctx = plc.x;
			int cty = plc.y;
			int orgx = sol.get_place(tile_id).x;
			int orgy = sol.get_place(tile_id).y;

			int use_radius = min(searchRadiu, max(abs(ctx - orgx), abs(cty - orgy)) / 2 + 1);
			int height = 2 * use_radius + 1;

					
			int biasX = abs(orgx - ctx) / 2;
			int biasY = abs(orgy - cty) / 2;
			ctx += orgx < ctx ? -biasX : orgx > ctx ? biasX : 0;
			cty += orgy < cty ? -biasY : orgy > cty ? biasY : 0;
			ctx = max(0, ctx); ctx = min(fabric_XSize - 1, ctx);
			cty = max(0, cty); cty = min(fabric_YSize - 1, cty);

			for (int p = 0; p < height; ++p) {

				int y = cty + use_radius - p;
				if (y < 0 || y >= fabric_YSize)continue;
				int x_left = max(0, p <= use_radius ? ctx - p : ctx - 2 * use_radius + p);
				int x_right = min(fabric_XSize - 1, p <= use_radius ? ctx + p : ctx + 2 * use_radius - p);

				for (int x = x_left; x <= x_right; ++x) {

					int swap_id = sol.get_tile_id(x, y);
					if (swap_id == tile_id || swap_id == adj_id)continue;

					pre_value = sol.getRealTimeCost(tile_id) + sol.getRealTimeCost(swap_id);
					aft_value = trySwapTile(tile_id, swap_id, Point(x, y), sol);

					if (aft_value < pre_value) {
						improved_i = true;
						swapOperation(tile_id, swap_id, Point(x, y), sol);
						sol.score += 2 * (aft_value - pre_value);
					}

					if (best_impro < pre_value - aft_value) {
						best_swap_id = swap_id;
						best_place = Point(x, y);
						best_impro = pre_value - aft_value;
					}
				}
			}
		}
		if(best_impro < 0 && (best_swap_id != -1 || best_place.x >= 0)){
			swapOperation(tile_id, best_swap_id, best_place, sol);
			sol.score -= 2 * best_impro;
		}
	}

	void local_search(Solution &sol,
		Solution &best_sol,
		int max_iteration = 1000,
		int searchRadiu = 2,
		int visibleInterval = 50,
		double timeLimited = 600.) {

		searchRadiu = 2;

		double norm = pow(100. * fabric_XSize * fabric_YSize, 2. / 3);

		clock_t start = clock();
		clock_t iterTime = clock();

		sol.init_tabu_list();

		vector<pair<bool, int>>improved((int)connections.size(), {true, 0});

		double pre_score = sol.score;
		int impTimes = 0;

		for (int iter = 1; iter <= max_iteration; ++iter) {
#ifdef LOCAL
			if (iter % visibleInterval == 0) {
				printf("Iteration progress: %3d / %3d | ", iter, max_iteration);
			}
#endif // LOCAL

			//searchRadiu = int(pow(sol.score / connections.size() / 6, 2. / 3));

			int _move = 0, _close = 0, _capture = 0;

			for (int i = 0; i < (int)connections.size(); ++i) {

				if (sol.is_tabu(0, i, iter)) continue;
				
				double pre_score = sol.score;

				int actionType = rand() % 1000;
				actionType = actionType < 990 ? 0 : actionType < 995 ? 1 : 2;

				if (actionType == 0) {
					moveAround(i, searchRadiu, sol);
					_move++;
				}
				else if (actionType == 1) {
					closeToNeighbors(i, searchRadiu, sol);
					_close++;
				}
				else {
					captureNeighbors(i, searchRadiu, sol);
					_capture++;
				}

				//if (pre_score >= sol.score) {
					sol.set_tabu(0, i, iter);
				//}
			}
			
#ifdef LOCAL
			if (iter % visibleInterval == 0) {
				printf("Action type = [%d, %d, %d] | ", _move, _close, _capture);
				double endtime = (double)(clock() - iterTime) / CLOCKS_PER_SEC;
				printf(" Time: %.3lfs | Score = %.6lf | ESPT [%d]\n", endtime, pow(sol.score, 2 / 3.), searchRadiu);
				saveIterationData << iter << ' ' << (double)(clock() - start) / CLOCKS_PER_SEC << ' ' 
								  << pow(sol.score, 2. / 3) << ' ' << norm / pow(sol.score, 2 / 3.) << endl;
				iterTime = clock();
			}
#endif // LOCAL

			//if (pre_score > sol.score) {
			//	impTimes++;
			//	if (impTimes > 5) {
			//		searchRadiu--;
			//	}
			//	if (searchRadiu < 5) {
			//		searchRadiu = 5;
			//	}
			//}
			//else if(pre_score < sol.score){
			//	impTimes = 0;
			//	searchRadiu += 1;
			//	if (searchRadiu > 10)searchRadiu = 5;
			//}
			//pre_score = sol.score;

			if (best_sol.score > sol.score) {
				best_sol = sol;
			}

			if ((double)(clock() - start) / CLOCKS_PER_SEC > timeLimited) {
				return;
			}

		}
	}

	/*************************************************/

	void parameterInitialization(unsigned randSeed) {
		srand(randSeed ? randSeed : (unsigned int)time(0));

		value1D5.resize(fabric_XSize + fabric_YSize + 2, 0);
		for (int i = 0; i < (int)value1D5.size(); ++i) {
			value1D5[i] = 1.0 * i * sqrt(i);
		}

		fabric.resize(fabric_XSize, vector<int>(fabric_YSize));
		for (int y = 0; y < fabric_YSize; ++y) {
			for (int x = 0; x < fabric_XSize; ++x) {
				fabric[x][y] = y * fabric_XSize + x;
			}
		}

		memset(is_used_tile, false, sizeof is_used_tile);
		memset(is_used_place, false, sizeof is_used_place);
		espt.clear();
		get_prisms();
	}

	Solution solve(const int &inv_sampling_step,
		const vector<int> &_cubeShape,
		const vector<int> &fabricShape,
		const int &tileCnt,
		const int &adapterCnt,
		const vector<vector<int>> &cnts,
		const vector<vector<int>> &pris_info,
		const vector<pair<int, int>> &init_place,

		int prePlacementCnt = 1,
		bool startWithCenter = false,
		int maxIteration = 1000,
		int searchRadiu = 2,
		unsigned randSeed = 0,
		bool localSearch = true,
		int visibleInterval = 1,
		double timeLimited = 600.,
		string initialAgorithm = "projection",
		string iterationsSaveFile = "Null") {

		cubeShape = _cubeShape;
		fabric_XSize = fabricShape[0];
		fabric_YSize = fabricShape[1];
		prisms_info = pris_info;
		tile_cnt = tileCnt;
		adapter_cnt = adapterCnt;
		connections = cnts;
		inverse_sampling_step = inv_sampling_step;

#ifdef LOCAL
		db2(tile_cnt, adapter_cnt);
		db2(fabric_XSize, fabric_YSize);
#endif // LOCAL

		parameterInitialization(randSeed); 

		Solution initSol = get_init_solution(init_place, prePlacementCnt, startWithCenter, initialAgorithm);

		initSol.score = evaluate(initSol);
		initSol.initRealTimeCost();
		//initSol.initCostBucket();

		Solution bestSol = initSol;

#ifdef LOCAL
		analyze_solution(bestSol);
		printf("Init solution score = %.6lf\n", pow(initSol.score, 2 / 3.));
#endif // LOCAL

		saveIterationData.open(iterationsSaveFile.c_str());

		saveIterationData << "0" << " 0.0 " << pow(initSol.score, 2. / 3) << ' ' 
						  << pow(100. * fabric_XSize * fabric_YSize, 2. / 3) / pow(initSol.score, 2. / 3) << endl;

		if (localSearch) {
			local_search(initSol, bestSol, maxIteration, searchRadiu, visibleInterval, timeLimited);
		}

#ifdef LOCAL
		printf("Best solution score = %.6lf\n", pow(bestSol.score, 2 / 3.));
		printf("Best solution score [For QAP Compare] = %.6lf\n", bestSol.score / 2);
		analyze_solution(bestSol);
#endif // LOCAL

		saveIterationData.close();

		return bestSol;
	}
}
