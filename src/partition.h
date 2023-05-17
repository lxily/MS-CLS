namespace Partition {
	int dimension;
	double alpha, beta;
	vector<int>cubeShape;
	vector<int>fabricShape;
	vector<double>heatMap;

	struct Coord {
		double x, y, z;
		Coord(double _x = 0, double _y = 0, double _z = 0) :x(_x), y(_y), z(_z) {}

		bool operator < (const Coord & c)const {
			return  dequal(z, c.z) ? dequal(y, c.y) ? x < c.x : y < c.y : z < c.z;
		}
		bool operator == (const Coord & c)const {
			return dequal(x, c.x) && dequal(y, c.y) && dequal(z, c.z);
		}
		bool operator != (const Coord & c)const {
			return !dequal(x, c.x) || !dequal(y, c.y) || !dequal(z, c.z);
		}

	};

	struct Tile {
		int id;				//Globally unique
		short resolution;
		short pos[3];	//Not a coordinate point, but a grid point

		vector<int>adjTile[6];//0-left,1-right,2-front,3-back,4-down,5-up

		Tile() {
		}

		~Tile() {
			for (int i = 0; i < 6; ++i) adjTile[i].clear();
		}

		bool operator < (const Tile & t)const {
			if (fabs(resolution - t.resolution) > 1e-6) return resolution < t.resolution;
			int s1 = pos[0] + pos[1] + pos[2];
			int s2 = t.pos[0] + t.pos[1] + t.pos[2];
			//if (s1 != s2) return s1 < s2;
			if (pos[2] - t.pos[2])return pos[2] < t.pos[2];
			if (pos[1] - t.pos[1])return pos[1] < t.pos[1];
			return pos[0] < t.pos[0];
		}
	};

	bool v_insert(vector<int> &v, const int val) {
		if (v.size() > 3)return false; //Up to four connections
		for (auto &x : v) if (x == val) return false;	//already exists
		v.push_back(val);
		return true;
	}

	bool v_erase(vector<int> &v, const int val) {
		bool success = false;
		short sz = (short)v.size();
		for (short i = 0; i < sz; ++i) {
			if (v[i] == val) {
				success = true;
				if (i != sz - 1) {
					v[i] = v[sz - 1];
				}
				v.pop_back();
				break;
			}
		}
		
		return success;
	}

	bool v_find(const vector<int> &v, const int val) {
		for (auto &x : v)if (x == val)return true;
		return false;
	}

	double resolution_of_heatMap(int x, int y, int z) {
		if (dimension == 2) {
			return heatMap[(cubeShape[0] + 1) * y + x];
		}
		else {
			return heatMap[(cubeShape[1] + 1) * (cubeShape[0] + 1) * z + (cubeShape[0] + 1) * y + x];
		}
	}

	double get_heat_value_in_heatMap(Coord pos) {
		int x_low = (int)floor(pos.x);
		int x_hig = (int)ceil(pos.x);

		int y_low = (int)floor(pos.y);
		int y_hig = (int)ceil(pos.y);

		int z_low = (int)floor(pos.z);
		int z_hig = (int)ceil(pos.z);

		//z-low
		///y-low
		////x-low
		double f_111 = resolution_of_heatMap(x_low, y_low, z_low);
		////x-high
		double f_211 = resolution_of_heatMap(x_hig, y_low, z_low);
		/////x-low-value
		double f_11 = pos.x - x_low < DOUBLE_EPSILON ? f_111 : x_hig - pos.x < DOUBLE_EPSILON ? f_211 : ((pos.x - x_low) * f_211 + (x_hig - pos.x) * f_111) / (x_hig - x_low);

		///y-hig
		////x-low
		double f_121 = resolution_of_heatMap(x_low, y_hig, z_low);
		////x-high
		double f_221 = resolution_of_heatMap(x_hig, y_hig, z_low);
		/////x-high-value
		double f_21 = pos.x - x_low < DOUBLE_EPSILON ? f_121 : x_hig - pos.x < DOUBLE_EPSILON ? f_221 : ((pos.x - x_low) * f_221 + (x_hig - pos.x) * f_121) / (x_hig - x_low);

		//////z-low-value
		double f_1 = pos.y - y_low < DOUBLE_EPSILON ? f_11 : y_hig - pos.y < DOUBLE_EPSILON ? f_21 : ((pos.y - y_low) * f_21 + (y_hig - pos.y) * f_11) / (y_hig - y_low);

		if (dimension == 2)return f_1;

		//z-hig
		///y-low
		////x-low
		double f_112 = resolution_of_heatMap(x_low, y_low, z_hig);
		////x-hig
		double f_212 = resolution_of_heatMap(x_hig, y_low, z_hig);
		/////x-low_value
		double f_12 = pos.x - x_low < DOUBLE_EPSILON ? f_112 : x_hig - pos.x < DOUBLE_EPSILON ? f_212 : ((pos.x - x_low) * f_212 + (x_hig - pos.x) * f_112) / (x_hig - x_low);

		///y-hig
		////x-low
		double f_122 = resolution_of_heatMap(x_low, y_hig, z_hig);
		////x-hig
		double f_222 = resolution_of_heatMap(x_hig, y_hig, z_hig);
		/////x-high-value
		double f_22 = pos.x - x_low < DOUBLE_EPSILON ? f_122 : x_hig - pos.x < DOUBLE_EPSILON ? f_222 : ((pos.x - x_low) * f_222 + (x_hig - pos.x) * f_122) / (x_hig - x_low);

		//////z-high-value
		double f_2 = pos.y - y_low < DOUBLE_EPSILON ? f_12 : y_hig - pos.y < DOUBLE_EPSILON ? f_22 : ((pos.y - y_low) * f_22 + (y_hig - pos.y) * f_12) / (y_hig - y_low);

		return pos.z - z_low < DOUBLE_EPSILON ? f_1 : z_hig - pos.z < DOUBLE_EPSILON ? f_2 : ((pos.z - z_low) * f_2 + (z_hig - pos.z) * f_1) / (z_hig - z_low);
	}

	// Integrate using trapezoid rule over heatmap at given y, z over given x range
	double integrate_x(double x_low, double x_high, double y, double z) {
		fatalif(x_low > x_high, "Integrating x wrong way!");

		int range_low = (int)ceil(x_low);
		int range_high = (int)floor(x_high);

		// Occurs when both are in same grid space
		if (range_low > range_high) {
			Coord low = { x_low, y, z };
			Coord high = { x_high, y, z };

			double low_val = get_heat_value_in_heatMap(low);
			double high_val = get_heat_value_in_heatMap(high);

			return 0.5 * (high_val + low_val) * (x_high - x_low);
		}

		double sum = 0.0;

		for (int i = range_low; i < range_high; i += 1) {
			Coord low = { (double)i, y, z };
			Coord high = { (double)i + 1.0, y, z };

			double low_val = get_heat_value_in_heatMap(low);
			double high_val = get_heat_value_in_heatMap(high);

			// Since units of 1, can just integrate directly without delta
			sum += 0.5 * (high_val + low_val);
		}

		double low_diff = (double)range_low - x_low;
		double high_diff = x_high - (double)range_high;

		// Add endpoints if they were not exactly on heatmap
		if (low_diff > DOUBLE_EPSILON) {
			Coord low = { x_low, y, z };
			Coord high = { (double)range_low, y, z };

			double low_val = get_heat_value_in_heatMap(low);
			double high_val = get_heat_value_in_heatMap(high);

			sum += 0.5 * (high_val + low_val) * low_diff;
		}

		if (high_diff > DOUBLE_EPSILON) {
			Coord low = { (double)range_high, y, z };
			Coord high = { x_high, y, z };

			double low_val = get_heat_value_in_heatMap(low);
			double high_val = get_heat_value_in_heatMap(high);

			sum += 0.5 * (high_val + low_val) * high_diff;
		}

		return sum;
	}

	// Double integrate heatmap using trapezoid rule at given z over given x, y range
	double integrate_y(double x_low, double x_high, double y_low, double y_high, double z) {

		fatalif(y_low > y_high, "Integrating y wrong way! [%.3lf > %.3lf]", y_low, y_high);

		int range_low = (int)ceil(y_low);
		int range_high = (int)floor(y_high);

		// Occurs when both are in same grid space
		if (range_low > range_high) {
			double low_val = integrate_x(x_low, x_high, y_low, z);
			double high_val = integrate_x(x_low, x_high, y_high, z);

			return 0.5 * (high_val + low_val) * (y_high - y_low);
		}

		double sum = 0.0;

		for (int i = range_low; i < range_high; i += 1) {
			double low_val = integrate_x(x_low, x_high, i, z);
			double high_val = integrate_x(x_low, x_high, i + 1, z);

			// Trapezoidal rule (unit step)
			sum += 0.5 * (high_val + low_val);
		}

		double low_diff = (double)range_low - y_low;
		double high_diff = y_high - (double)range_high;

		// Add endpoints if they were not exactly on heatmap
		if (low_diff > DOUBLE_EPSILON) {

			double low_val = integrate_x(x_low, x_high, y_low, z);
			double high_val = integrate_x(x_low, x_high, range_low, z);

			sum += 0.5 * (high_val + low_val) * low_diff;
		}

		if (high_diff > DOUBLE_EPSILON) {

			double low_val = integrate_x(x_low, x_high, range_high, z);
			double high_val = integrate_x(x_low, x_high, y_high, z);

			sum += 0.5 * (high_val + low_val) * high_diff;
		}

		return sum;
	}

	// Triple integrate heatmap using trapezoid rule over given x, y, z range
	double integrate_z(double x_low, double x_high, double y_low, double y_high, double z_low, double z_high) {
		fatalif(z_low > z_high, "Integrating z wrong way!");

		int range_low = (int)ceil(z_low);
		int range_high = (int)floor(z_high);

		// Occurs when both are in same grid space
		if (range_low > range_high) {
			double low_val = integrate_y(x_low, x_high, y_low, y_high, z_low);
			double high_val = integrate_y(x_low, x_high, y_low, y_high, z_high);

			return 0.5 * (high_val + low_val) * (z_high - z_low);
		}

		double sum = 0.0;

		for (int i = range_low; i < range_high; i += 1) {
			double low_val = integrate_y(x_low, x_high, y_low, y_high, i);
			double high_val = integrate_y(x_low, x_high, y_low, y_high, i + 1);

			// Trapezoidal rule (unit step)
			sum += 0.5 * (high_val + low_val);
		}

		double low_diff = (double)range_low - z_low;
		double high_diff = z_high - (double)range_high;

		// Add endpoints if they were not exactly on heatmap
		if (low_diff > DOUBLE_EPSILON) {

			double low_val = integrate_y(x_low, x_high, y_low, y_high, z_low);
			double high_val = integrate_y(x_low, x_high, y_low, y_high, range_low);

			sum += 0.5 * (high_val + low_val) * low_diff;
		}

		if (high_diff > DOUBLE_EPSILON) {

			double low_val = integrate_y(x_low, x_high, y_low, y_high, range_high);
			double high_val = integrate_y(x_low, x_high, y_low, y_high, z_high);

			sum += 0.5 * (high_val + low_val) * high_diff;
		}

		return sum;
	}

	// Find max overshoot between target heatmap and solution resolution at given y, z over given x range
	double find_max_heat_value_x(double x_low, double x_high, double y, double z) {
		fatalif(x_low > x_high, "x wrong way!");

		int range_low = (int)ceil(x_low);
		int range_high = (int)floor(x_high);

		// Occurs when both are in same grid space
		if (range_low > range_high) {
			Coord low = { x_low, y, z };
			Coord high = { x_high, y, z };

			double low_val = get_heat_value_in_heatMap(low);
			double high_val = get_heat_value_in_heatMap(high);

			return max(low_val, high_val);
		}

		double max_overshoot = 0;

		for (int i = range_low; i <= range_high; i += 1) {
			Coord curr_coord = { (double)i, y, z };
			double val = get_heat_value_in_heatMap(curr_coord);

			max_overshoot = max(max_overshoot, val);
		}

		double low_diff = (double)range_low - x_low;
		double high_diff = x_high - (double)range_high;

		// Add endpoints if they were not exactly on heatmap
		if (low_diff > DOUBLE_EPSILON) {
			Coord curr_coord = { (double)x_low, y, z };

			double val = get_heat_value_in_heatMap(curr_coord);

			max_overshoot = max(max_overshoot, val);
		}

		if (high_diff > DOUBLE_EPSILON) {
			Coord curr_coord = { (double)x_high, y, z };

			double val = get_heat_value_in_heatMap(curr_coord);

			max_overshoot = max(max_overshoot, val);
		}

		return max_overshoot;
	}

	// Find max overshoot between target heatmap and solution resolution at given z over given x, y range
	double find_max_heat_value_y(double x_low, double x_high, double y_low, double y_high, double z) {

		fatalif(y_low > y_high, "y wrong way!");

		int range_low = (int)ceil(y_low);
		int range_high = (int)floor(y_high);

		// Occurs when both are in same grid space
		if (range_low > range_high) {
			double low_val = find_max_heat_value_x(x_low, x_high, y_low, z);
			double high_val = find_max_heat_value_x(x_low, x_high, y_high, z);

			return max(low_val, high_val);
		}

		double max_overshoot = 0.0;

		for (int i = range_low; i <= range_high; i += 1) {
			double val = find_max_heat_value_x(x_low, x_high, i, z);

			max_overshoot = max(max_overshoot, val);
		}

		double low_diff = (double)range_low - y_low;
		double high_diff = y_high - (double)range_high;

		// Add endpoints if they were not exactly on heatmap
		if (low_diff > DOUBLE_EPSILON) {
			double val = find_max_heat_value_x(x_low, x_high, y_low, z);

			max_overshoot = max(max_overshoot, val);
		}

		if (high_diff > DOUBLE_EPSILON) {
			double value = find_max_heat_value_x(x_low, x_high, y_high, z);

			max_overshoot = max(max_overshoot, value);
		}

		return max_overshoot;
	}

	// Find max overshoot between target heatmap and solution resolution over given x, y, z range
	double find_max_heat_value_z(double x_low, double x_high, double y_low, double y_high, double z_low, double z_high) {

		fatalif(z_low > z_high, "z wrong way!");

		int range_low = (int)ceil(z_low);
		int range_high = (int)floor(z_high);

		// Occurs when both are in same grid space
		if (range_low > range_high) {
			double low_val = find_max_heat_value_y(x_low, x_high, y_low, y_high, z_low);
			double high_val = find_max_heat_value_y(x_low, x_high, y_low, y_high, z_high);

			return max(low_val, high_val);
		}

		double max_overshoot = 0.0;

		for (int i = range_low; i <= range_high; i += 1) {
			double val = find_max_heat_value_y(x_low, x_high, y_low, y_high, i);

			max_overshoot = max(max_overshoot, val);
		}

		double low_diff = (double)range_low - z_low;
		double high_diff = z_high - (double)range_high;

		// Add endpoints if they were not exactly on heatmap
		if (low_diff > DOUBLE_EPSILON) {
			double val = find_max_heat_value_y(x_low, x_high, y_low, y_high, z_low);

			max_overshoot = max(max_overshoot, val);
		}

		if (high_diff > DOUBLE_EPSILON) {
			double val = find_max_heat_value_y(x_low, x_high, y_low, y_high, z_high);

			max_overshoot = max(max_overshoot, val);
		}

		return max_overshoot;
	}

	struct Solution {
		using t_iter = set<Tile>::iterator;

		set<Tile>tiles;
		vector<t_iter>idToIter;

		double score;
		int adapter_cnt;				//Initialize in get_init_solution()

		vector<int>sampleShape;			//Initialize in init_sample_map

		int inverse_sample_step;		//Initialize in init_sample_map()
		double max_scale_limitation;

		int fabric_limitation;

		vector<vector<double>>resAvg; 	//Initialize in init_resolution_info()
		vector<vector<double>>resMax;	//same as resMax

		vector<int>IDSPOOR;				//Initialize in get_init_solution()

		set<pair<double, int>>max_scales;

		vector<vector<int>>empty_prisms;

		Solution() {
			score = 0;
			adapter_cnt = 0;
			set<Tile>().swap(tiles);
			vector<int>().swap(IDSPOOR);
			vector<t_iter>().swap(idToIter);
			set<pair<double, int>>().swap(max_scales);
			vector<vector<double>>().swap(resAvg);
			vector<vector<double>>().swap(resMax);
			vector<vector<int>>().swap(empty_prisms);
		}

		void clear() {
			score = 0;
			adapter_cnt = 0;
			set<Tile>().swap(tiles);
			vector<int>().swap(IDSPOOR);
			vector<t_iter>().swap(idToIter);
			set<pair<double, int>>().swap(max_scales);
			vector<vector<double>>().swap(resAvg);
			vector<vector<double>>().swap(resMax);
			vector<vector<int>>().swap(empty_prisms);
		}

		void insert(const Tile &t) {
			if (t.id >= (int)idToIter.size()) {
				idToIter.resize(t.id + 1);
			}
			idToIter[t.id] = tiles.insert(t).first;
		}

		void remove(const Tile &t) {
			tiles.erase(idToIter[t.id]);
		}

		void init_sampleShape(int iss) {
			inverse_sample_step = iss;

			sampleShape = { cubeShape[0] * inverse_sample_step / 10,
							cubeShape[1] * inverse_sample_step / 10,
							cubeShape[2] * inverse_sample_step / 10 };
		}

		void init_resolution_info() {
			fatalif(!sampleShape.size(), "No sampleShape.");
			vector<vector<double>>().swap(resAvg);
			vector<vector<double>>().swap(resMax);

			int mxLength = min(sampleShape[0], min(sampleShape[1], sampleShape[2]));
			int miResolution = int(floor(log2(mxLength) + DOUBLE_EPSILON)) + 1;

			int total_grids = sampleShape[0] * sampleShape[1] * sampleShape[2];
			resAvg.resize(total_grids, vector<double>(miResolution, 0));
			resMax.resize(total_grids, vector<double>(miResolution, 0));

			for (int r = 0; r < miResolution; ++r) {
				for (int x = 0; x < (int)sampleShape[0]; ++x) {
					for (int y = 0; y < (int)sampleShape[1]; ++y) {
						for (int z = 0; z < (int)sampleShape[2]; ++z) {
							int id = sampleMapId(x, y, z);
							if (r == 0) {
								Coord leftDown = { 1.* x, 1.* y, 1.* z };
								Coord rightUp = { x + 1., y + 1., z + 1. };

								leftDown = coord_from_sampleMap_to_heatMap(leftDown);
								rightUp = coord_from_sampleMap_to_heatMap(rightUp);

								resAvg[id][0] = integrate_z(leftDown.x, rightUp.x, leftDown.y, rightUp.y, leftDown.z, rightUp.z);
								resMax[id][0] = find_max_heat_value_z(leftDown.x, rightUp.x, leftDown.y, rightUp.y, leftDown.z, rightUp.z);
							}
							else if (x + (1 << r) <= sampleShape[0] &&
								y + (1 << r) <= sampleShape[1] &&
								z + (1 << r) <= sampleShape[2]) {
								int radiu = (1 << (r - 1));

								double a111 = resAvg[id][r - 1];
								double m111 = resMax[id][r - 1];

								int id211 = sampleMapId(x + radiu, y, z);
								double a211 = resAvg[id211][r - 1];
								double m211 = resMax[id211][r - 1];

								int id121 = sampleMapId(x, y + radiu, z);
								double a121 = resAvg[id121][r - 1];
								double m121 = resMax[id121][r - 1];

								int id221 = sampleMapId(x + radiu, y + radiu, z);
								double a221 = resAvg[id221][r - 1];
								double m221 = resMax[id221][r - 1];

								int id112 = sampleMapId(x, y, z + radiu);
								double a112 = resAvg[id112][r - 1];
								double m112 = resMax[id112][r - 1];

								int id212 = sampleMapId(x + radiu, y, z + radiu);
								double a212 = resAvg[id212][r - 1];
								double m212 = resMax[id212][r - 1];

								int id122 = sampleMapId(x, y + radiu, z + radiu);
								double a122 = resAvg[id122][r - 1];
								double m122 = resMax[id122][r - 1];

								int id222 = sampleMapId(x + radiu, y + radiu, z + radiu);
								double a222 = resAvg[id222][r - 1];
								double m222 = resMax[id222][r - 1];

								resAvg[id][r] = (a111 + a211 + a121 + a221 + a112 + a212 + a122 + a222);
								resMax[id][r] = max(m111, max(m211, max(m121, max(m221, max(m112, max(m212, max(m122, m222)))))));
							}
						}
					}
				}
			}
		}

		inline int sampleMapId(int x, int y, int z) {
			return x + y * sampleShape[0] + z * sampleShape[0] * sampleShape[1];
		}

		inline double max_scale() {
			return (--max_scales.end())->first;
		}

		Coord coord_from_sampleMap_to_heatMap(Coord p) {
			double nx = p.x * 10 / inverse_sample_step;
			double ny = p.y * 10 / inverse_sample_step;
			double nz = p.z * 10 / inverse_sample_step;
			return { nx , ny, nz };
		}

		int query_resolution(const int &tile_id) {
			fatalif(idToIter[tile_id]->id < 0, "Query_resolution error.[By tile_id]");
			return idToIter[tile_id]->resolution;
		}

		double query_max_require_resolution(const short *pos, short resol) {
			/*int leftdownID = sampleMapId(pos[0], pos[1], pos[2]);
			return resMax[leftdownID][resol];*/

			Coord leftDown = { 1.* pos[0], 1.* pos[1], 1.* pos[2] };
			Coord rightUp = { 1. * pos[0] + (1 << resol), 1. * pos[1] + (1 << resol), 1.* pos[2] + (1 << resol) };

			leftDown = coord_from_sampleMap_to_heatMap(leftDown);
			rightUp = coord_from_sampleMap_to_heatMap(rightUp);

			return find_max_heat_value_z(leftDown.x, rightUp.x, leftDown.y, rightUp.y, leftDown.z, rightUp.z);
		}

		double query_max_require_resolution(int tile_id) {
			return query_max_require_resolution(idToIter[tile_id]->pos, idToIter[tile_id]->resolution);
		}

		double query_avg_require_resolution(const short *pos, short resol) {
			/*int leftdownID = sampleMapId(pos[0], pos[1], pos[2]);
			return resAvg[leftdownID][resol];*/

			Coord leftDown = { 1.* pos[0], 1.* pos[1], 1.* pos[2] };
			Coord rightUp = { 1. * pos[0] + (1 << resol), 1. * pos[1] + (1 << resol), 1.* pos[2] + (1 << resol) };

			leftDown = coord_from_sampleMap_to_heatMap(leftDown);
			rightUp = coord_from_sampleMap_to_heatMap(rightUp);

			return integrate_z(leftDown.x, rightUp.x, leftDown.y, rightUp.y, leftDown.z, rightUp.z);
		}

		double query_avg_require_resolution(int tile_id) {
			return query_avg_require_resolution(idToIter[tile_id]->pos, idToIter[tile_id]->resolution);
		}

		double calc_tile_covered_resolution(int tile_id, int increased = 0) {
			const Tile &t = *idToIter[tile_id];
			const short *p = t.pos;
			double radius = 1 << (t.resolution + increased);
			Coord leftDown = coord_from_sampleMap_to_heatMap({ (double)p[0], (double)p[1], (double)p[2] });
			Coord rightUp = coord_from_sampleMap_to_heatMap({ p[0] + radius, p[1] + radius, p[2] + radius });

			return 1. / radius * (rightUp.x - leftDown.x) * (rightUp.y - leftDown.y) * (rightUp.z - leftDown.z);
		}

		double max_scale_of_covered(int tile_id, int increased = 0) {
			const Tile &t = *idToIter[tile_id];

			double max_res = query_max_require_resolution(t.pos, t.resolution + increased);

			double real_res = 1. / (1 << (t.resolution + increased));
			return 1. + max(0., (max_res - real_res)) / real_res;
		}

		double calc_tile_score(int tile_id, int increased = 0) {
			const Tile &t = *idToIter[tile_id];

			double real_inte_res = calc_tile_covered_resolution(t.id, increased);	//Integration resolution after merging
			double requier_inte_res = query_avg_require_resolution(t.pos, t.resolution + increased);

			return requier_inte_res / real_inte_res;
		}

		double calc_area_score(const short p[3], int resolution) {
			double radius = 1 << resolution;
			Coord leftDown = coord_from_sampleMap_to_heatMap({ (double)p[0], (double)p[1], (double)p[2] });
			Coord rightUp = coord_from_sampleMap_to_heatMap({ p[0] + radius, p[1] + radius, p[2] + radius });

			double real_inte_res = 1. / radius * (rightUp.x - leftDown.x) * (rightUp.y - leftDown.y) * (rightUp.z - leftDown.z);
			double requier_inte_res = integrate_z(leftDown.x, rightUp.x, leftDown.y, rightUp.y, leftDown.z, rightUp.z);

			return requier_inte_res / real_inte_res;
		}

		bool is_mergeable_tile(int tile_id, bool regular) {
			const Tile &t = *idToIter[tile_id];

			if (regular &&
				(t.pos[0] + 0) % min((1 << (t.resolution + 1)), 128) ||
				(t.pos[1] + 0) % min((1 << (t.resolution + 1)), 128) ||
				(t.pos[2] + 0) % min((1 << (t.resolution + 1)), 128)) {
				return  false;
			}
			//if (!regular && (t.pos[0] % 2 || t.pos[1] % 2 || t.pos[2] % 2)) return false;

			if (!connection_validity(t)) return false;

			double maxScare = max_scale_of_covered(t.id, 1);//resolution preference - Avoid overshoot
			if (maxScare > max_scale_limitation + DOUBLE_EPSILON) {
				return false;
			}

			vector<int>ids = get_mergable_tiles(t);

			vector<int> merged_adj[6];
			if ((merge_adj_info({ ids[0], ids[2], ids[4], ids[6] }, 0, merged_adj[0])) == inf ||
				(merge_adj_info({ ids[1], ids[3], ids[5], ids[7] }, 1, merged_adj[1])) == inf ||
				(merge_adj_info({ ids[0], ids[1], ids[4], ids[5] }, 2, merged_adj[2])) == inf ||
				(merge_adj_info({ ids[2], ids[3], ids[6], ids[7] }, 3, merged_adj[3])) == inf ||
				(merge_adj_info({ ids[0], ids[1], ids[2], ids[3] }, 4, merged_adj[4])) == inf ||
				(merge_adj_info({ ids[4], ids[5], ids[6], ids[7] }, 5, merged_adj[5])) == inf) {
				return false;
			}

			return true;
		}

		bool is_mergeable_tile(int tile_id, int adapter_increase[6], vector<int> merged_adj[6]) {
			const Tile &t = *idToIter[tile_id];

			if (t.pos[0] % (1 << (t.resolution + 1)) ||
				t.pos[1] % (1 << (t.resolution + 1)) ||
				t.pos[2] % (1 << (t.resolution + 1))) {
				return  false;
			}

			if (!connection_validity(t)) return false;

			double maxScare = max_scale_of_covered(t.id, 1);//resolution preference - Avoid overshoot
			if (maxScare > max_scale_limitation + DOUBLE_EPSILON) {
				return false;
			}

			vector<int>ids = get_mergable_tiles(t);

			if ((adapter_increase[0] = merge_adj_info({ ids[0], ids[2], ids[4], ids[6] }, 0, merged_adj[0])) == inf ||
				(adapter_increase[1] = merge_adj_info({ ids[1], ids[3], ids[5], ids[7] }, 1, merged_adj[1])) == inf ||
				(adapter_increase[2] = merge_adj_info({ ids[0], ids[1], ids[4], ids[5] }, 2, merged_adj[2])) == inf ||
				(adapter_increase[3] = merge_adj_info({ ids[2], ids[3], ids[6], ids[7] }, 3, merged_adj[3])) == inf ||
				(adapter_increase[4] = merge_adj_info({ ids[0], ids[1], ids[2], ids[3] }, 4, merged_adj[4])) == inf ||
				(adapter_increase[5] = merge_adj_info({ ids[4], ids[5], ids[6], ids[7] }, 5, merged_adj[5])) == inf) {
				return false;
			}

			return true;
		}

		int is_splitable_tile(int tile_id) {
			const Tile &t = *idToIter[tile_id];

			int adapter_c = 0;
			if (t.resolution == 0)return inf;
			for (int i = 0; i < 6; ++i) {
				if (!t.adjTile[i].size())continue;
				int adjRes = idToIter[*t.adjTile[i].begin()]->resolution;
				if (t.resolution < adjRes) {
					return inf;
				}
				else if (t.resolution == adjRes) {
					++adapter_c;
				}
				else {
					--adapter_c;
				}
			}

			return adapter_c;
		}

		int delta_score_of_split(int tile_id) {
			const Tile &t = *idToIter[tile_id];

			double pre_score = calc_tile_score(tile_id);

			const short dir[8][3] = {
				{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
				{0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1},
			};

			short radius = 1 << (t.resolution - 1);
			double aft_score = 0.0;
			for (int i = 0; i < 8; ++i) {
				short p[3] = {
					t.pos[0] + dir[i][0] * radius,
					t.pos[1] + dir[i][1] * radius,
					t.pos[2] + dir[i][2] * radius,
				};
				aft_score += calc_area_score(p, t.resolution - 1);
			}

			return int((aft_score - pre_score) * 1000000);
		}

		int delta_score_of_merge(int tile_id) {
			const Tile &t = *idToIter[tile_id];

			double pre_score = 0;

			vector<int>ids = get_mergable_tiles(t);

			for (int i = 0; i < 8; ++i) {
				pre_score += calc_tile_score(ids[i]);
			}

			double aft_score = calc_tile_score(tile_id, 1);

			return int((aft_score - pre_score) * 1000000);
		}

		int delta_adapter_of_split(int tile_id) {
			const Tile &t = *idToIter[tile_id];
			int delta = 0;
			for (int i = 0; i < 6; ++i) {
				switch (t.adjTile[i].size()) {
				case 0: break;
				case 1: 
					delta += t.resolution == idToIter[*t.adjTile[i].begin()]->resolution; break;
				case 4: --delta; break;
				default: fatalif(1, "ERROR [ delta_adapter_of_split ].");
				}
			}
			return delta;
		}

		int delta_adapter_of_merge(int tile_id) {
			const Tile &t = *idToIter[tile_id];
			int delta = 0;
			for (int i = 0; i < 6; ++i) {
				switch (t.adjTile[i].size()) {
				case 0: break;
				case 1: 
					if (t.resolution == idToIter[*t.adjTile[i].begin()]->resolution) {
						delta += 1;
					}
					else {
						delta -= 1;
					}
				case 4: break;
				default: fatalif(1, "ERROR [ delta_adapter_of_merge ].");
				}
			}
			return delta;
		}

		int find_improve(int tile_id) {
			const Tile &t = *idToIter[tile_id];

			vector<int>merged_adj[6];	//Neighbor resolution constraint
			int adapter_increase[6];

			if (!is_mergeable_tile(tile_id, adapter_increase, merged_adj)) {
				return -1;
			}
			int decrease = delta_score_of_merge(tile_id);

			int adapter_1 = 0;
			for (int i = 0; i < 6; ++i) {
				adapter_1 += adapter_increase[i];
			}

			int increase = -inf, split_id = -1;
			for (int i = 0; i < 6; ++i) {
				if (merged_adj[i].size() != 1)continue;
				int adj_id = merged_adj[i][0], adapter_2;
				if ((adapter_2 = is_splitable_tile(adj_id)) == inf) {
					continue;
				}
				if (adapter_1 + adapter_2 + 2 + total_tile_used() > fabric_limitation)continue;
				
				int inc_score = delta_score_of_split(adj_id);
				if (increase < inc_score) {
					increase = inc_score;
					split_id = adj_id;
				}
			}

			if (decrease + increase <= 0)return -1;
			return split_id;
		}

		void local_search() {
			for (t_iter it = tiles.begin(); it != tiles.end();) {
				int cur_di = it->id;
				int adj_id = find_improve(cur_di);
				if (adj_id < 0) {
					++it; continue;
				}
				int split_id = split_one_tile(*idToIter[adj_id]);
				fatalif(split_id < 0, "Tile restricted and cannot be split.");

				IDSPOOR.push_back(adj_id);
				remove(*idToIter[adj_id]);

				int merged_id = merge_one_tile(*it, false);
				fatalif(merged_id < 0, "Should be mergeable.");

				IDSPOOR.push_back(cur_di);
				tiles.erase(it++);	//delete and get next

#ifdef LOCAL
				printf("Improve by merge tile %d and split tile %d, total pe used = %d after operation, score = %.6lf [%.6lf].\n", 
					cur_di, adj_id, total_tile_used(), score / max_scale(), score / max_scale() / fabric_limitation);
#endif // LOCAL

				fatalif(total_tile_used() > fabric_limitation, "Error total_tile_used [%d] > fabric_limitation [%d].", total_tile_used(), fabric_limitation);
			}
		}

		vector<int>get_mergable_tiles(const Tile &t) {
			int id_0 = t.id;
			int id_1 = *(t.adjTile[1].begin());
			int id_2 = *(idToIter[id_0]->adjTile[3].begin());
			int id_3 = *(idToIter[id_2]->adjTile[1].begin());
			int id_4 = *(idToIter[id_0]->adjTile[5].begin());
			int id_5 = *(idToIter[id_4]->adjTile[1].begin());
			int id_6 = *(idToIter[id_4]->adjTile[3].begin());
			int id_7 = *(idToIter[id_6]->adjTile[1].begin());
			return {
				id_0, id_1, id_2, id_3, id_4, id_5, id_6, id_7
			};
		}

		bool connection_validity(const Tile &t) {
			//2 3	//6 7
			//0 1	//4 5

			if (t.adjTile[1].size() != 1 || t.adjTile[3].size() != 1 || t.adjTile[5].size() != 1)return false;

			const Tile &t_1 = *idToIter[*t.adjTile[1].begin()];
			if (t_1.adjTile[0].size() != 1 || t_1.adjTile[3].size() != 1 || t_1.adjTile[5].size() != 1)return false;

			const Tile &t_2 = *idToIter[*t_1.adjTile[3].begin()];
			if (t_2.adjTile[0].size() != 1 || t_2.adjTile[2].size() != 1 || t_2.adjTile[5].size() != 1)return false;

			const Tile &t_3 = *idToIter[*t_2.adjTile[0].begin()];
			if (t_3.adjTile[1].size() != 1 || t_3.adjTile[2].size() != 1 || t_3.adjTile[5].size() != 1)return false;

			const Tile &t_4 = *idToIter[*t.adjTile[5].begin()];
			if (t_4.adjTile[1].size() != 1 || t_4.adjTile[3].size() != 1 || t_4.adjTile[4].size() != 1)return false;

			const Tile &t_5 = *idToIter[*t_4.adjTile[1].begin()];
			if (t_5.adjTile[0].size() != 1 || t_5.adjTile[3].size() != 1 || t_5.adjTile[4].size() != 1)return false;

			const Tile &t_6 = *idToIter[*t_5.adjTile[3].begin()];
			if (t_6.adjTile[0].size() != 1 || t_6.adjTile[2].size() != 1 || t_6.adjTile[4].size() != 1)return false;

			const Tile &t_7 = *idToIter[*t_6.adjTile[0].begin()];
			if (t_7.adjTile[1].size() != 1 || t_7.adjTile[2].size() != 1 || t_7.adjTile[4].size() != 1)return false;
			return true;
		}

		bool check_connection_vadility() {
			for (auto &t : tiles) {
				for (int i = 0; i < 6; ++i) {
					for (auto &id : t.adjTile[i]) {
						if (!v_find(idToIter[id]->adjTile[i ^ 1], t.id)) {
							printf("[A]Connection Error: ");
							printf("Tile %d have connection to tile %d on face %d, the opposite is not true.\n", t.id, id, i);
							exit(-1);
						}
					}
				}
			}
		}

		bool check_connection_vadility(const Tile &t) {
			for (int i = 0; i < 6; ++i) {
				for (auto &id : t.adjTile[i]) {
					if (!v_find(idToIter[id]->adjTile[i ^ 1], t.id)) {
						return false;
						printf("[B]Connection Error: ");
						printf("Tile %d have connection to tile %d on face %d, the opposite is not true.\n", t.id, id, i);
						exit(-1);
					}
				}
			}
			return true;
		}

		int merge_adj_info(vector<int> tile_ids, int face, vector<int> &res) {
			const Tile &t = *idToIter[tile_ids[0]];

			if (t.adjTile[face].size() == 0) {
				for (int i = 1; i < (int)tile_ids.size(); ++i) {
					if (idToIter[tile_ids[i]]->adjTile[face].size())return inf;
				}
				return 0;
			}
			else {
				if (t.adjTile[face].size() != 1)return inf;
				int adjID = *t.adjTile[face].begin();
				v_insert(res, adjID);
				int adjRes = idToIter[adjID]->resolution;

				for (int i = 1; i < (int)tile_ids.size(); ++i) {
					if (idToIter[tile_ids[i]]->adjTile[face].size() != 1)return inf;
					if (idToIter[*idToIter[tile_ids[i]]->adjTile[face].begin()]->resolution != adjRes) return inf;
					if (t.resolution < adjRes) {
						if (*idToIter[tile_ids[i]]->adjTile[face].begin() != adjID) {
							return inf;
						}
					}
					else if (t.resolution == adjRes) {
						v_insert(res, *idToIter[tile_ids[i]]->adjTile[face].begin());
					}
					else {
						return inf;
					}
				}
				return t.resolution < adjRes ? -1 : res.size() == 4 ? 1 : inf;
			}
		}

		int merge_one_tile(const Tile &t, bool regular) {
			//Successfully return the merged id, else -1
			//you should reomve 't' from 'tiles' after execute this function
			if (regular &&
				(t.pos[0] + 0) % min((1 << (t.resolution + 1)), 128) ||
				(t.pos[1] + 0) % min((1 << (t.resolution + 1)), 128) ||
				(t.pos[2] + 0) % min((1 << (t.resolution + 1)), 128)) {
				return  -1;
			}
			//if (!regular && (t.pos[0] % 2 || t.pos[1] % 2 || t.pos[2] % 2)) return -1;

			if (!connection_validity(t)) return -1;

			double maxScare = max_scale_of_covered(t.id, 1);//resolution preference - Avoid overshoot
			if (maxScare > max_scale_limitation + DOUBLE_EPSILON) {
				return -1;
			}

			vector<int>ids = get_mergable_tiles(t);

			vector<int>merged_adj[6];	//Neighbor resolution constraint
			int adapter_increase[6];

			if ((adapter_increase[0] = merge_adj_info({ ids[0], ids[2], ids[4], ids[6] }, 0, merged_adj[0])) == inf ||
				(adapter_increase[1] = merge_adj_info({ ids[1], ids[3], ids[5], ids[7] }, 1, merged_adj[1])) == inf ||
				(adapter_increase[2] = merge_adj_info({ ids[0], ids[1], ids[4], ids[5] }, 2, merged_adj[2])) == inf ||
				(adapter_increase[3] = merge_adj_info({ ids[2], ids[3], ids[6], ids[7] }, 3, merged_adj[3])) == inf ||
				(adapter_increase[4] = merge_adj_info({ ids[0], ids[1], ids[2], ids[3] }, 4, merged_adj[4])) == inf ||
				(adapter_increase[5] = merge_adj_info({ ids[4], ids[5], ids[6], ids[7] }, 5, merged_adj[5])) == inf) {
				return -1;
			}

			for (int i = 0; i < 8; ++i) {
				score -= calc_tile_score(ids[i]);
			}

			for (int i = 0; i < 6; ++i) {
				adapter_cnt += adapter_increase[i];
			}

			//merge!
			for (int j = 0; j < 8; ++j) {
				pair<double, int> p_erase = { max_scale_of_covered(idToIter[ids[j]]->id), idToIter[ids[j]]->id };
				max_scales.erase(max_scales.find(p_erase));
				if (j) {
					IDSPOOR.push_back(ids[j]);				//ID recycling
					remove(*idToIter[ids[j]]);
				}	//it Finally delete for j = 0

			}

			Tile addTile;
			addTile.id = IDSPOOR.back();
			IDSPOOR.pop_back();
			addTile.pos[0] = t.pos[0];
			addTile.pos[1] = t.pos[1];
			addTile.pos[2] = t.pos[2];
			addTile.resolution = t.resolution + 1;
			for (int j = 0; j < 6; ++j) {
				addTile.adjTile[j] = merged_adj[j];
			}
			insert(addTile);
			score += calc_tile_score(addTile.id);

			max_scales.insert({ max_scale_of_covered(addTile.id), addTile.id });

			//Update neighbor information
			for (int j = 0; j < 6; ++j) {
				int aface = j ^ 1;
				for (auto tile_id : merged_adj[j]) {
					Tile &adj_t = const_cast<Tile &>(*idToIter[tile_id]);
					for (int k = 0; k < 8; ++k) {
						if (v_find(adj_t.adjTile[aface], ids[k])) {
							v_erase(adj_t.adjTile[aface], ids[k]);
						}
					}
					fatalif(!v_insert(adj_t.adjTile[aface], addTile.id), "Error");
				}
			}

			return addTile.id;
		}

		int split_tile_helper(const Tile &t, int face, vector<Tile> &newTiles, vector<int> adjIdx) {
			int adapterIncreasCnt = 0;

			int aface = face ^ 1;
			if (t.adjTile[face].size() == 1) {
				adapterIncreasCnt += 1;

				Tile &adj_t = const_cast<Tile &>(*idToIter[*t.adjTile[face].begin()]);
				v_erase(adj_t.adjTile[aface], t.id);

				for (auto &idx : adjIdx) {
					v_insert(newTiles[idx].adjTile[face], adj_t.id);
					v_insert(adj_t.adjTile[aface], newTiles[idx].id);
				}
				fatalif(adj_t.adjTile[aface].size() != 4, "Four connections should be generated after increasing the resolution.[%d]",
					(int)adj_t.adjTile[aface].size());
			}
			else if (t.adjTile[face].size() == 4) {
				adapterIncreasCnt -= 1;

				for (auto &id : t.adjTile[face]) {
					Tile &adj_t = const_cast<Tile &>(*idToIter[id]);
					fatalif(adj_t.adjTile[aface].size() != 1, "Should be one.[%d]", (int)adj_t.adjTile[aface].size());

					v_erase(adj_t.adjTile[aface], t.id);

					fatalif(adj_t.adjTile[aface].size(), "Should be zero.[%d]", (int)adj_t.adjTile[aface].size());

					for (auto &idx : adjIdx) {
						bool leftRightAdj = (face == 0 || face == 1) && idToIter[id]->pos[1] == newTiles[idx].pos[1]
							&& idToIter[id]->pos[2] == newTiles[idx].pos[2];
						bool frontBackAdj = (face == 2 || face == 3) && idToIter[id]->pos[0] == newTiles[idx].pos[0]
							&& idToIter[id]->pos[2] == newTiles[idx].pos[2];
						bool  upDownAdj = (face == 4 || face == 5) && idToIter[id]->pos[1] == newTiles[idx].pos[1]
							&& idToIter[id]->pos[0] == newTiles[idx].pos[0];

						if (leftRightAdj || frontBackAdj || upDownAdj) {
							v_insert(newTiles[idx].adjTile[face], id);
							v_insert(adj_t.adjTile[aface], newTiles[idx].id);
						}
					}
					fatalif(adj_t.adjTile[aface].size() != 1, "A single connection should be produced after increasing the resolution.[%d]",
						(int)adj_t.adjTile[aface].size());
				}
			}
			else if (t.adjTile[face].size()) {
				fatalif(1, "[split_one_tile] Neighbor connection error.");
			}

			return adapterIncreasCnt;
		}

		int split_one_tile(const Tile &t) {	//Note: It is legal before splitting by default
			//Successfully return the splited tile id at leftdown, else -1

			if (t.resolution == 0)return -1;
			for (int i = 0; i < 6; ++i) {
				if (!t.adjTile[i].size())continue;
				int adjRes = idToIter[*t.adjTile[i].begin()]->resolution;
				if (t.resolution < adjRes) {
					return -1;
				}
			}

			//2 3	//6 7
			//0 1	//4 5
			const int dir[8][3] = {
				{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
				{0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1},
			};
			int radius = 1 << (t.resolution - 1);
			vector<Tile>newTiles(8);

			for (int i = 0; i < 8; ++i) {
				fatalif(IDSPOOR.size() == 0, "[split_one_tile] IDSPOOR do not have enough ID.");

				newTiles[i].id = IDSPOOR.back();
				IDSPOOR.pop_back();
				newTiles[i].resolution = t.resolution - 1;
				newTiles[i].pos[0] = t.pos[0] + dir[i][0] * radius;
				newTiles[i].pos[1] = t.pos[1] + dir[i][1] * radius;
				newTiles[i].pos[2] = t.pos[2] + dir[i][2] * radius;
			}

			//Internal connection
			{
				v_insert(newTiles[0].adjTile[1], newTiles[1].id);
				v_insert(newTiles[0].adjTile[3], newTiles[2].id);
				v_insert(newTiles[0].adjTile[5], newTiles[4].id);

				v_insert(newTiles[1].adjTile[0], newTiles[0].id);
				v_insert(newTiles[1].adjTile[3], newTiles[3].id);
				v_insert(newTiles[1].adjTile[5], newTiles[5].id);

				v_insert(newTiles[2].adjTile[1], newTiles[3].id);
				v_insert(newTiles[2].adjTile[2], newTiles[0].id);
				v_insert(newTiles[2].adjTile[5], newTiles[6].id);

				v_insert(newTiles[3].adjTile[0], newTiles[2].id);
				v_insert(newTiles[3].adjTile[2], newTiles[1].id);
				v_insert(newTiles[3].adjTile[5], newTiles[7].id);

				v_insert(newTiles[4].adjTile[1], newTiles[5].id);
				v_insert(newTiles[4].adjTile[3], newTiles[6].id);
				v_insert(newTiles[4].adjTile[4], newTiles[0].id);

				v_insert(newTiles[5].adjTile[0], newTiles[4].id);
				v_insert(newTiles[5].adjTile[3], newTiles[7].id);
				v_insert(newTiles[5].adjTile[4], newTiles[1].id);

				v_insert(newTiles[6].adjTile[1], newTiles[7].id);
				v_insert(newTiles[6].adjTile[2], newTiles[4].id);
				v_insert(newTiles[6].adjTile[4], newTiles[2].id);

				v_insert(newTiles[7].adjTile[0], newTiles[6].id);
				v_insert(newTiles[7].adjTile[2], newTiles[5].id);
				v_insert(newTiles[7].adjTile[4], newTiles[3].id);
			}

			int adapter_increase[6];

			//left neighbor
			adapter_increase[0] = split_tile_helper(t, 0, newTiles, { 0, 2, 4, 6 });

			//right neighbor
			adapter_increase[1] = split_tile_helper(t, 1, newTiles, { 1, 3, 5, 7 });

			//front neighbor
			adapter_increase[2] = split_tile_helper(t, 2, newTiles, { 0, 1, 4, 5 });

			//back neighbor
			adapter_increase[3] = split_tile_helper(t, 3, newTiles, { 2, 3, 6, 7 });

			//down neighbor
			adapter_increase[4] = split_tile_helper(t, 4, newTiles, { 0, 1, 2, 3 });

			// up neighbor
			adapter_increase[5] = split_tile_helper(t, 5, newTiles, { 4, 5, 6, 7 });

			for (int i = 0; i < 6; ++i) {
				adapter_cnt += adapter_increase[i];
			}

			max_scales.erase(max_scales.find({ max_scale_of_covered(t.id), t.id }));
			for (int i = 0; i < 8; ++i) {
				insert(newTiles[i]);
				score += calc_tile_score(newTiles[i].id);
				max_scales.insert({ max_scale_of_covered(newTiles[i].id), newTiles[i].id });
			}

			score -= calc_tile_score(t.id);
			return newTiles[0].id;
		}

		vector<int> split_require(const Tile &t) {
			fatalif(t.resolution <= 0, "Can not split.");

			vector<int>ret = { t.id };

			if (t.resolution == 0)return { -1 };
			for (int i = 0; i < 6; ++i) {
				if (!t.adjTile[i].size())continue;
				int adjRes = idToIter[*t.adjTile[i].begin()]->resolution;
				if (t.resolution < adjRes) {
					vector<int>nxt_rq = split_require(*idToIter[*t.adjTile[i].begin()]);
					ret.insert(ret.end(), nxt_rq.begin(), nxt_rq.end());
				}
			}

			return ret;
		}

		pair<vector<int>, vector<int>> split_compulsorily(const Tile t) {
			vector<int>ret;

			vector<int>s_req = split_require(t);
			unordered_map<int, bool>have_split;
			vector<int>real_req;
			for (int i = (int)s_req.size() - 1; i >= 0; --i) {
				if (have_split[s_req[i]])continue;
				have_split[s_req[i]] = true;
				real_req.push_back(s_req[i]);
			}

			for (int i = 0; i < (int)real_req.size(); ++i) {

				int split_id = split_one_tile(*idToIter[real_req[i]]);
				fatalif(split_id < 0, "Tile restricted and cannot be split.");

				IDSPOOR.push_back(real_req[i]);
				remove(*idToIter[real_req[i]]);

				ret.push_back(split_id);
			}

			return {real_req, ret};
		}

		void merge_tiles(int target_res, bool regular) {
			//check_connection_vadility();

			for (set<Tile>::iterator it = tiles.begin(); it != tiles.end();) {

				if (it->resolution > target_res)break;

				if (it->resolution != target_res) {
					++it; continue;
				}

				int merID = merge_one_tile(*it, regular);

				if (merID >= 0) {
					IDSPOOR.push_back(it->id);
					tiles.erase(it++);	//delete and get next
				}
				else {
					++it;
				}
			}

			//check_connection_vadility();
		}

		void merge_tiles_in_score_order(bool regular, double ratio, bool recycle) {
			bool surplus = false;

			if(recycle) empty_tiles_recycle();

			int limited = int(fabric_limitation * ratio);

			vector<pair<int, int>>psid(idToIter.size());
			auto cmp = [&](const int &a, const int &b) {
				//int avg1 = psid[a].first / (psid[a].second * 0.5 + 7);
				//int avg2 = psid[b].first / (psid[b].second * 0.5 + 7);
				//if (avg1 != avg2)return avg1 < avg2;
				//if (psid[a].second != psid[b].second) return psid[a].second < psid[b].second;
				if (psid[a].first != psid[b].first)return psid[a].first > psid[b].first;
				if (psid[a].second != psid[b].second) return psid[a].second < psid[b].second;

				return a > b;
			};

			set<int, decltype(cmp)>mergeQue(cmp);
			for (auto &t : tiles) {
				if (is_mergeable_tile(t.id, regular)) {
					psid[t.id] = {
						delta_score_of_merge(t.id),
						delta_adapter_of_merge(t.id)
					};
					mergeQue.insert(t.id);
				}
			}

			while (total_tile_used() > limited && mergeQue.size()) {
				int tile_id = *(mergeQue.begin());
				mergeQue.erase(tile_id);

				if (!is_mergeable_tile(tile_id, regular)) continue;

				//int merID = merge_one_tile(*idToIter[tile_id], regular);
				//fatalif(merID < 0, "Error: unmergeable tile %d.", tile_id);

				vector<int>ids = get_mergable_tiles(*idToIter[tile_id]);
				int merID = merge_one_tile(*idToIter[tile_id], regular);
				fatalif(merID < 0, "Error: unmergeable tile %d.", tile_id);
				for (auto id : ids) { if (mergeQue.count(id))mergeQue.erase(id); }

				remove(*idToIter[tile_id]);
				IDSPOOR.push_back(tile_id);

				if (is_mergeable_tile(merID, regular)) {
					psid[merID] = {
						delta_score_of_merge(merID),
						delta_adapter_of_merge(merID)
					};
					mergeQue.insert(merID);
				}

				for (int i = 0; i < 6; ++i) {
					const auto & adjs = idToIter[merID]->adjTile[i];
					for (const int adj_id : adjs) {
						if (mergeQue.find(adj_id) != mergeQue.end()) {
							mergeQue.erase(adj_id);
						}
						if (is_mergeable_tile(adj_id, regular)) {
							psid[adj_id] = {
								delta_score_of_merge(adj_id),
								delta_adapter_of_merge(adj_id)
							};
							mergeQue.insert(adj_id);
						}
					}
				}
			}
		}
		
		void merge_to_local_optima(bool regular) {
			int pre_tiles;
			int pre_adapters;
			int merge_resol = 0;

			do {
#ifdef  LOCAL
				db3(tiles.size(), adapter_cnt, score / max_scale());
#endif //  LOCAL

				pre_tiles = (int)tiles.size();
				pre_adapters = adapter_cnt;
				merge_tiles(merge_resol++, regular);

			} while (pre_tiles != (int)tiles.size() || pre_adapters != adapter_cnt);
		}

		void delete_a_empty_tile(int tile_id) {
			fatalif(query_max_require_resolution(tile_id) > DOUBLE_EPSILON, "Can't remove a real tile.");
			Tile &t = const_cast<Tile &>(*idToIter[tile_id]);
			for (int i = 0; i < 6; ++i) {
				for (auto &adj_id : t.adjTile[i]) {
					Tile &adj_t = const_cast<Tile &>(*idToIter[adj_id]);
					/*if (!adj_t.adjTile[i ^ 1].size()) {
						db("Error");
					}*/
					v_erase(adj_t.adjTile[i ^ 1], tile_id);
				}
				vector<int>().swap(t.adjTile[i]);
			}

			vector<int>epri;
			epri.push_back(-1);
			epri.push_back(t.pos[0] * 10);
			epri.push_back(t.pos[1] * 10);
			epri.push_back(t.pos[2] * 10);
			for (int i = 0; i < 3; ++i) {
				epri.push_back(1 << t.resolution);
			}
			empty_prisms.push_back(epri);
		}

		vector<bool> get_removable_tiles() {
			//Get the tiles that cannot be deleted according to the neighbor relationship
			vector<bool>removable(idToIter.size(), true);

			priority_queue<pair<int, int>>notEmpty;
			for (auto &t : tiles) {
				if (query_max_require_resolution(t.id) > DOUBLE_EPSILON) {
					notEmpty.push({ t.resolution, t.id });
					removable[t.id] = false;
				}
			}

			while (!notEmpty.empty()) {
				const Tile &t = *idToIter[notEmpty.top().second];
				notEmpty.pop();

				for (int i = 0; i < 6; ++i) {
					if (!t.adjTile[i].size())continue;
					fatalif(t.adjTile[i].size() != 1 && t.adjTile[i].size() != 4, "Neighbor relationship error [%d].", (int)t.adjTile[i].size());
					vector<int>::const_iterator it = t.adjTile[i].begin();

					if (t.adjTile[i].size() == 1) {
						const Tile &adj_t = *idToIter[*t.adjTile[i].begin()];
						if (adj_t.resolution > t.resolution && !removable[adj_t.id]) {
							for (auto &id : adj_t.adjTile[i ^ 1]) {
								if (removable[id]) {
									notEmpty.push({ idToIter[id]->resolution, id });
									removable[id] = false;
								}
							}
						}
					}
					else {
						int adj_0 = *it++; bool ne_0 = (query_max_require_resolution(adj_0) > DOUBLE_EPSILON);
						int adj_1 = *it++; bool ne_1 = (query_max_require_resolution(adj_1) > DOUBLE_EPSILON);
						int adj_2 = *it++; bool ne_2 = (query_max_require_resolution(adj_2) > DOUBLE_EPSILON);
						int adj_3 = *it++; bool ne_3 = (query_max_require_resolution(adj_3) > DOUBLE_EPSILON);

						if (ne_0 || ne_1 || ne_2 || ne_3) {
							if (!ne_0 && removable[adj_0]) {
								notEmpty.push({ idToIter[adj_0]->resolution, adj_0 });
								removable[adj_0] = false;
							}
							if (!ne_1 && removable[adj_1]) {
								notEmpty.push({ idToIter[adj_1]->resolution, adj_1 });
								removable[adj_1] = false;
							}
							if (!ne_2 && removable[adj_2]) {
								notEmpty.push({ idToIter[adj_2]->resolution, adj_2 });
								removable[adj_2] = false;
							}
							if (!ne_3 && removable[adj_3]) {
								notEmpty.push({ idToIter[adj_3]->resolution, adj_3 });
								removable[adj_3] = false;
							}
						}
					}
				}
			}

			return removable;
		}

		pair<int, int> get_tiles_adapters_reduce(vector<bool> &removable) {
			int tile_reduce = 0;
			int adaptersReduce = 0;
			for (auto &t : tiles) {
				if (removable[t.id])tile_reduce += 1;
				for (int i = 0; i < 6; ++i) {
					if (removable[t.id] && (int)t.adjTile[i].size() == 4) {
						adaptersReduce += 1;
					}
					else if (!removable[t.id] && (int)t.adjTile[i].size() == 4) {
						adaptersReduce += removable[*t.adjTile[i].begin()];
					}
				}
			}
			return make_pair(tile_reduce, adaptersReduce);
		}

		void empty_tiles_recycle() {

			vector<bool> removable = get_removable_tiles();
			pair<int, int>reduce = get_tiles_adapters_reduce(removable);
			adapter_cnt -= reduce.second;

			for (t_iter it = tiles.begin(); it != tiles.end();) {
				if (removable[it->id]) {
					max_scales.erase(max_scales.find({ max_scale_of_covered(it->id), it->id }));
					delete_a_empty_tile(it->id);
					IDSPOOR.emplace_back(it->id);
					tiles.erase(it++);
				}
				else {
					it++;
				}
			}
#ifdef LOCAL
			printf("\tRecycle %d empty tiles and %d adapters, total PEs used = %d now, score = %.5lf.\n",
				reduce.first, reduce.second, total_tile_used(), score / max_scale());
#endif // LOCAL
		}

		int count_max_scale() {
			int max_scale_cnt = 0;
			double mscale = max_scale();
			for (auto &p : max_scales) {
				max_scale_cnt += dequal(p.first, mscale);
			}
			return max_scale_cnt;
		}

		int print_tiles_used_info() {
			int useless_tiles = 0;

			unordered_map<int, int>mp;
			for (auto &t : tiles) {
				useless_tiles += query_max_require_resolution(t.id) < DOUBLE_EPSILON;
				if (query_max_require_resolution(t.id) > DOUBLE_EPSILON) {
					mp[t.resolution]++;
				}
			}

			int max_scale_cnt = count_max_scale();

#ifdef LOCAL
			printf("\n=========================Tile Information=========================\n\n");
			printf("Merge infromation: "); for (int i = 0; mp[i]; ++i)printf("[%d, %d] ", i, mp[i]); putchar('\n');
			printf("Empty prisms cnt = %d, Useless tiles cnt = %d\n", (int)empty_prisms.size(), useless_tiles);
			printf("Compute tiles cnt = %d, Adapters cnt = %d\n", (int)tiles.size(), adapter_cnt);
			printf("Total PEs used: %d, there are %d empty tiles but can not remove.\n",
				(int)tiles.size() + adapter_cnt, useless_tiles);
			printf("Total Score = %.3lf, Scale = %.3lf [%d / %d], Pure Accuracy = %.3lf, score = %.5lf.\n",
				score, max_scale(), max_scale_cnt, (int)max_scales.size(), score / max_scale(), score / max_scale() / fabric_limitation);
			printf("\n==================================================================\n\n");
#endif // LOCAL

			return (int)tiles.size() + adapter_cnt;
		}

		bool split_max_scale_tiles(double aimSocre = 1.0) {
			bool reduced = false;

			empty_tiles_recycle();

			int limited = fabric_limitation;
			if (total_tile_used() <= limited && score / max_scale() / limited >= aimSocre) {
				return false;
			}

			while (true) {
				pair<double, int>p = *(--max_scales.end());

				/*if (p.first <= 2 + DOUBLE_EPSILON) break;*/

				if (dequal(p.first, 1.))break;

				if (total_tile_used() <= limited && score / max_scale() / limited >= aimSocre) {
					return false;
				}

				vector<int>split_ids = split_compulsorily(*idToIter[p.second]).second;

				if (total_tile_used() > limited) {//split success
					for (int j = (int)split_ids.size() - 1; j >= 0; --j) {
						int merged_id = merge_one_tile(*idToIter[split_ids[j]], false);
						fatalif(merged_id < 0, "Should be mergeable.");

						remove(*idToIter[split_ids[j]]);
						IDSPOOR.push_back(split_ids[j]);
					}
					break;
				}

				reduced = true;
			}

			return reduced;
		}

		bool use_surplus_tiles(double aimSocre = 1.0) {

			bool surplus = false;

			//check_connection_vadility();

			empty_tiles_recycle();

			int limited = fabric_limitation;

			if (total_tile_used() <= limited && score / max_scale() / limited >= aimSocre) {
				return false;
			}

			//auto 

			struct Comparator {
				int dScore, dAdapter, seed;
				Comparator(){}
				Comparator(int _ds, int _da, int _sd) {
					dScore = _ds; dAdapter = _da; seed = _sd;
				}
			};

			vector<Comparator>psid(idToIter.size());
			auto cmp = [&](const int &a, const int &b) {
				int avg1 = psid[a].dScore / (psid[a].dAdapter / 2. + 7);
				int avg2 = psid[b].dScore / (psid[b].dAdapter / 2. + 7);
				//if (fabs(avg1 - avg2) > DOUBLE_EPSILON)return avg1 < avg2;

				if (fabs(psid[a].dScore - psid[b].dScore) > DOUBLE_EPSILON) {
					return psid[a].dScore < psid[b].dScore;
				}
				
				//return psid[a].dAdapter > psid[b].dAdapter;

				if (psid[a].dAdapter != psid[b].dAdapter) return psid[a].dAdapter > psid[b].dAdapter;
				//if (avg1 != avg2)return avg1 < avg2;
				return psid[a].seed < psid[b].seed;
				return a > b;
			};

			set<int, decltype(cmp)>splitQue(cmp);
			for (auto &t : tiles) {
				if (t.resolution > 0 && is_splitable_tile(t.id) != inf) {
					psid[t.id] = {
						//int(calc_tile_score(t.id) * 1000000),
						delta_score_of_split(t.id),
						delta_adapter_of_split(t.id),
						rand() * rand()
					};
					splitQue.insert(t.id);
				}
			}

			while(total_tile_used() < limited && splitQue.size()) {
				if (score / max_scale() / limited >= aimSocre) {
					empty_tiles_recycle(); break;
				}
				
				int tile_id = *(--splitQue.end());

				//db3(psid[tile_id].dScore, psid[tile_id].dAdapter, idToIter[tile_id]->resolution);

				vector<int>upd;
				for (int i = 0; i < 6; ++i) {
					for (const int fadj : idToIter[tile_id]->adjTile[i]) {
						const Tile &t = *idToIter[fadj];
						if (t.resolution > 0) {
							upd.push_back(fadj);
						}
					}
				}

				auto res = split_compulsorily(*idToIter[tile_id]);
				vector<int>origin_ids = res.first;
				vector<int>split_ids = res.second;

				for (auto &orx : origin_ids) {
					if (splitQue.find(orx) != splitQue.end()) {
						splitQue.erase(orx);
					}
				}

				if (origin_ids.size() == 1) {
					for (auto &adj_id : upd) {
						if (splitQue.find(adj_id) != splitQue.end()) {
							splitQue.erase(adj_id);
						}
						if (is_splitable_tile(adj_id) != inf) {
							psid[adj_id] = {
								//int(calc_tile_score(adj_id) * 1000000),
								delta_score_of_split(adj_id),
								delta_adapter_of_split(adj_id),
								rand() * rand()
							};
							splitQue.insert(adj_id);
						}
					}
				}

				if (total_tile_used() > limited) {//split success

					for (int j = (int)split_ids.size() - 1; j >= 0; --j) {
						int merged_id = merge_one_tile(*idToIter[split_ids[j]], false);
						fatalif(merged_id < 0, "Should be mergeable.");

						remove(*idToIter[split_ids[j]]);
						IDSPOOR.push_back(split_ids[j]);

						if (j == (int)split_ids.size() - 1)continue;
						psid[merged_id] = { 
							//int(calc_tile_score(merged_id) * 1000000),
							delta_score_of_split(merged_id), 
							delta_adapter_of_split(merged_id) ,
							rand() * rand()
						};
						splitQue.insert(merged_id);

					}
					//continue;
					break;
				}
				else {
					for (int j = (int)split_ids.size() - 1; j >= 0; --j) {
						vector<int>ids = get_mergable_tiles(*idToIter[split_ids[j]]);
						for (auto id: ids) {
							if (idToIter[id]->resolution > 0) {
								psid[id] = { 
									//int(calc_tile_score(id) * 1000000),
									delta_score_of_split(id), 
									delta_adapter_of_split(id) ,
									rand() * rand()
								};
								splitQue.insert(id);
							}
						}
					}
					surplus = true;
				}

			}

			return  surplus;
		}

		void adjustment(double _max_scale_limitation = 1., double aimSocre = 1.0) {
			max_scale_limitation = _max_scale_limitation;

			merge_to_local_optima(true); 

			//merge_tiles_in_score_order(true, 1, false);

			//merge_to_local_optima();

			//print_tiles_used_info();

//#ifdef LOCAL
//			printf("*Start reduce max_scale, max_scale = %.10f now.\n", max_scale());
//#endif // LOCAL
//
//			while (true) {
//				//print_tiles_used_info();
//				if (!split_max_scale_tiles(aimSocre))break;
//			}
//
//#ifdef LOCAL
//			printf("*Reduce max_scale OK, max_scale = %.10f now.\n", max_scale());
//#endif // LOCAL

#ifdef LOCAL
			printf("*Start use surplus tiles, total tiles = %d now.\n", total_tile_used());
#endif // LOCAL

			//empty_tiles_recycle();

			while (use_surplus_tiles(aimSocre)) {
				//print_tiles_used_info();
			}

			//check_connection_vadility();

			int iteration = 1;

			while (iteration--) {
#ifdef LOCAL
				printf("*Start merge low-score tiles, total tiles = %d now.\n", total_tile_used());
#endif // LOCAL

				merge_tiles_in_score_order(false, 0.5, false);

#ifdef LOCAL
				printf("*Merge low-score tiles OK, total tiles = %d now.\n", total_tile_used());
#endif // LOCAL

				//check_connection_vadility();

#ifdef LOCAL
				printf("*Start use surplus tiles, total tiles = %d now.\n", total_tile_used());
#endif // LOCAL

				while (use_surplus_tiles(aimSocre)) {
					//print_tiles_used_info();
				}

#ifdef LOCAL
				printf("*Use surplus tiles OK, total tiles = %d now.\n", total_tile_used());
#endif // LOCAL
			}

//#ifdef LOCAL
//			printf("*Start local search, total tiles = %d now.\n", total_tile_used());
//#endif // LOCAL
//			
//			local_search();
//
//#ifdef LOCAL
//			printf("*Local search OK, total tiles = %d now.\n", total_tile_used());
//#endif // LOCAL
		}

		inline int total_tile_used() {
			return (int)tiles.size() + adapter_cnt;
		}

		void format_to_prisms(vector<vector<int>> &pris) {
			for (int i = 0; i < (int)empty_prisms.size(); ++i) {
				pris.push_back(empty_prisms[i]);
			}

			for (auto &t : tiles) {
				vector<int>pri;

				pri.push_back(t.resolution);

				pri.push_back(t.pos[0] * 10);
				pri.push_back(t.pos[1] * 10);
				pri.push_back(t.pos[2] * 10);

				for (int i = 0; i < dimension; ++i) {
					pri.push_back(1);
				}

				pris.push_back(pri);
			}
		}

		void print() {
			for (auto &t : tiles) {
				printf("ID: %d, resolution: %d, avg_res: %.3lf, max_res: %.3lf, Coord: [%d, %d, %d]\n",
					t.id, t.resolution, query_avg_require_resolution(t.id), query_max_require_resolution(t.id),
					t.pos[0], t.pos[1], t.pos[2]);
				for (int i = 0; i < 6; ++i) {
					printf("Adj face: %d: ", i);
					for (auto &val : t.adjTile[i]) {
						printf("%d ", val);
					}
					printf("\n");
				}
			}
		}
	};

	vector<int> get_possible_sample_steps() {
		fatalif(!heatMap.size(), "No heatMap!\n");

		vector<int> ret;

		int xSize = cubeShape[0];
		int ySize = cubeShape[1];
		int zSize = cubeShape[2];

		long long estimateTiles = xSize * ySize * zSize;

		for (int m_iss = 1, step = 1;; ++step) {
			if (xSize * step % 10 || ySize * step % 10 || zSize * step % 10)continue;
			if (estimateTiles * step * step * step / 1000 <= fabricShape[0] * fabricShape[1]) {
				m_iss = step;
			}
			else if (ret.size() < 1000) {
				if (!ret.size() && (xSize * m_iss % 10 == 0 || ySize * m_iss % 10 == 0 || zSize * m_iss % 10 == 0)) {
					ret.push_back(m_iss);
				}
				ret.push_back(step);
			}
			else {
				break;
			}
		}

		return ret;
	}

	vector<int> get_empty_space_range(Solution &sol) {
		int xSize = sol.sampleShape[0];
		int ySize = sol.sampleShape[1];
		int zSize = sol.sampleShape[2];

		vector<int>ret;

		//0-left
		bool not_find_real = true;
		for (int x = 0; x < xSize && not_find_real; ++x) {
			for (int z = 0; z < zSize && not_find_real; ++z) {
				for (int y = 0; y < ySize && not_find_real; ++y) {
					const short pos[3] = { (short)x, (short)y, (short)z };
					if (sol.query_max_require_resolution(pos, 0) > DOUBLE_EPSILON) {
						ret.push_back(max(-1, x - 1));
						not_find_real = false;
						break;
					}
				}
			}
		}

		//1-right
		not_find_real = true;
		for (int x = xSize - 1; x >= 0 && not_find_real; --x) {
			for (int z = 0; z < zSize && not_find_real; ++z) {
				for (int y = 0; y < ySize && not_find_real; ++y) {
					const short pos[3] = { (short)x, (short)y, (short)z };
					if (sol.query_max_require_resolution(pos, 0) > DOUBLE_EPSILON) {
						ret.push_back(min(xSize, x + 1));
						not_find_real = false;
						break;
					}
				}
			}
		}

		//2-front
		not_find_real = true;
		for (int y = 0; y < ySize && not_find_real; ++y) {
			for (int z = 0; z < zSize && not_find_real; ++z) {
				for (int x = 0; x < xSize && not_find_real; ++x) {
					const short pos[3] = { (short)x, (short)y, (short)z };
					if (sol.query_max_require_resolution(pos, 0) > DOUBLE_EPSILON) {
						ret.push_back(max(-1, y - 1));
						not_find_real = false;
						break;
					}
				}
			}
		}

		//3-back
		not_find_real = true;
		for (int y = ySize - 1; y >= 0 && not_find_real; --y) {
			for (int z = 0; z < zSize && not_find_real; ++z) {
				for (int x = 0; x < xSize && not_find_real; ++x) {
					const short pos[3] = { (short)x, (short)y, (short)z };
					if (sol.query_max_require_resolution(pos, 0) > DOUBLE_EPSILON) {
						ret.push_back(min(ySize, y + 1));
						not_find_real = false;
						break;
					}
				}
			}
		}

		//4-down
		not_find_real = true;
		for (int z = 0; z < zSize && not_find_real; ++z) {
			for (int y = 0; y < ySize && not_find_real; ++y) {
				for (int x = 0; x < xSize && not_find_real; ++x) {
					const short pos[3] = { (short)x, (short)y, (short)z };
					if (sol.query_max_require_resolution(pos, 0) > DOUBLE_EPSILON) {
						ret.push_back(max(-1, z - 1));
						not_find_real = false;
						break;
					}
				}
			}
		}

		//5-up
		not_find_real = true;
		for (int z = zSize - 1; z >= 0 && not_find_real; --z) {
			for (int y = 0; y < ySize && not_find_real; ++y) {
				for (int x = 0; x < xSize && not_find_real; ++x) {
					const short pos[3] = { (short)x, (short)y, (short)z };
					if (sol.query_max_require_resolution(pos, 0) > DOUBLE_EPSILON) {
						ret.push_back(min(zSize, z + 1));
						not_find_real = false;
						break;
					}
				}
			}
		}

		return ret;
	}

	vector<int> delete_surrounding_space(Solution &sol) {
		int xSize = sol.sampleShape[0];
		int ySize = sol.sampleShape[1];
		int zSize = sol.sampleShape[2];

		vector<int>empty_space = get_empty_space_range(sol);
		if (empty_space[0] >= 0) {	//left-sapce
			vector<int>epri;
			epri.push_back(-1);
			epri.push_back(0);
			epri.push_back(0);
			epri.push_back(0);
			epri.push_back(empty_space[0] + 1);
			epri.push_back(ySize);
			epri.push_back(zSize);
			sol.empty_prisms.push_back(epri);
		}

		if (empty_space[1] < xSize) {	//right-sapce
			vector<int>epri;
			epri.push_back(-1);
			epri.push_back(empty_space[1] * 10);
			epri.push_back(0);
			epri.push_back(0);
			epri.push_back(xSize - empty_space[1]);
			epri.push_back(ySize);
			epri.push_back(zSize);
			sol.empty_prisms.push_back(epri);
		}

		if (empty_space[2] >= 0) {	//front-sapce
			vector<int>epri;
			epri.push_back(-1);
			epri.push_back(empty_space[0] * 10 + 10);
			epri.push_back(0);
			epri.push_back(0);
			epri.push_back(empty_space[1] - empty_space[0] - 1);
			epri.push_back(empty_space[2] + 1);
			epri.push_back(zSize);
			sol.empty_prisms.push_back(epri);
		}

		if (empty_space[3] < ySize) {	//back-sapce
			vector<int>epri;
			epri.push_back(-1);
			epri.push_back(empty_space[0] * 10 + 10);
			epri.push_back(empty_space[3] * 10);
			epri.push_back(0);
			epri.push_back(empty_space[1] - empty_space[0] - 1);
			epri.push_back(ySize - empty_space[3]);
			epri.push_back(zSize);
			sol.empty_prisms.push_back(epri);
		}

		if (empty_space[4] >= 0) {	//down-sapce
			vector<int>epri;
			epri.push_back(-1);
			epri.push_back(empty_space[0] * 10 + 10);
			epri.push_back(empty_space[2] * 10 + 10);
			epri.push_back(0);
			epri.push_back(empty_space[1] - empty_space[0] - 1);
			epri.push_back(empty_space[3] - empty_space[2] - 1);
			epri.push_back(empty_space[4] + 1);
			sol.empty_prisms.push_back(epri);
		}

		if (empty_space[5] < zSize) {	//up-sapce
			vector<int>epri;
			epri.push_back(-1);
			epri.push_back(empty_space[0] * 10 + 10);
			epri.push_back(empty_space[2] * 10 + 10);
			epri.push_back(empty_space[5] * 10);
			epri.push_back(empty_space[1] - empty_space[0] - 1);
			epri.push_back(empty_space[3] - empty_space[2] - 1);
			epri.push_back(zSize - empty_space[5]);
			sol.empty_prisms.push_back(epri);
		}

		return empty_space;
	}

	void get_jump_over_list(vector<bool>&jumpOverList, int xyz_min[3], int xyz_max[3], int jump_res, Solution &ret) {
		int xSize = ret.sampleShape[0];
		int ySize = ret.sampleShape[1];
		int zSize = ret.sampleShape[2];

		int radius = 1 << jump_res;

		for (int z = xyz_min[2]; z < xyz_max[2]; ++z) {
			for (int y = xyz_min[1]; y < xyz_max[1]; ++y) {
				for (int x = xyz_min[0]; x < xyz_max[0]; ++x) {
					int id = ret.sampleMapId(x, y, z);
					short pos[3] = { (short)x, (short)y, (short)z };

					if (jumpOverList[id] || x % radius || y % radius || z % radius)continue;

					if (x + radius > xyz_max[0] || y + radius > xyz_max[1] || z + radius > xyz_max[2])continue;

					Coord bias_min = {
						//1. * pos[0], 1. * pos[1], 1. * pos[2]
						max(1. * xyz_min[0], 1. * pos[0]),
						max(1. * xyz_min[1], 1. * pos[1]),
						max(1. * xyz_min[2], 1. * pos[2])
					};
					Coord bias_max = {
						//1. * pos[0] + radius, 1. * pos[1] + radius, 1. * pos[2] + radius,
						min(1. * xyz_max[0], 1. * pos[0] + radius),
						min(1. * xyz_max[1], 1. * pos[1] + radius),
						min(1. * xyz_max[2], 1. * pos[2] + radius)
					};

					bias_min = ret.coord_from_sampleMap_to_heatMap(bias_min);
					bias_max = ret.coord_from_sampleMap_to_heatMap(bias_max);

					if (find_max_heat_value_z(bias_min.x, bias_max.x, bias_min.y, bias_max.y, bias_min.z, bias_max.z) < DOUBLE_EPSILON) {
						for (int z_d = 0; z_d < radius; ++z_d) {
							for (int y_d = 0; y_d < radius; ++y_d) {
								for (int x_d = 0; x_d < radius; ++x_d) {
									int cur_id = ret.sampleMapId(x + x_d, y + y_d, z + z_d);
									jumpOverList[cur_id] = true;
								}
							}
						}

						vector<int>epri;
						epri.push_back(-1);
						for (int i = 0; i < 3; ++i)epri.push_back(pos[i] * 10);
						for (int i = 0; i < 3; ++i)epri.push_back(radius);

						ret.empty_prisms.push_back(epri);
					}


				}
			}
		}
	}

	Solution get_init_solution(int iss) {
		Solution ret;

		ret.init_sampleShape(iss);
		//ret.init_resolution_info();

		int xSize = ret.sampleShape[0];
		int ySize = ret.sampleShape[1];
		int zSize = ret.sampleShape[2];

		vector<int>empty_space = { -1, xSize, -1, ySize, -1, zSize };
		//empty_space = delete_surrounding_space(ret);

		int dir[6][3] = { {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1} };

		int xyz_min[3] = { empty_space[0] + 1, empty_space[2] + 1, empty_space[4] + 1 };
		int xyz_max[3] = { empty_space[1], empty_space[3], empty_space[5] };

		vector<bool>jumpOverList(xSize * ySize * zSize);
		get_jump_over_list(jumpOverList, xyz_min, xyz_max, 5, ret);
		get_jump_over_list(jumpOverList, xyz_min, xyz_max, 4, ret);
		get_jump_over_list(jumpOverList, xyz_min, xyz_max, 3, ret);
		//get_jump_over_list(jumpOverList, xyz_min, xyz_max, 2, ret);

#ifdef LOCAL
		db(ret.empty_prisms.size());
#endif // LOCAL
		for (int z = xyz_min[2]; z < xyz_max[2]; ++z) {
			for (int y = xyz_min[1]; y < xyz_max[1]; ++y) {
				for (int x = xyz_min[0]; x < xyz_max[0]; ++x) {
					Tile t;
					t.id = ret.sampleMapId(x, y, z);
					t.resolution = 0;
					t.pos[0] = x;
					t.pos[1] = y;
					t.pos[2] = z;

					if (jumpOverList[t.id])continue;

					for (int d = 0; d < 6; ++d) {
						int nx = x + dir[d][0];
						int ny = y + dir[d][1];
						int nz = z + dir[d][2];
						if (nx < xyz_min[0] || nx >= xyz_max[0] ||
							ny < xyz_min[1] || ny >= xyz_max[1] ||
							nz < xyz_min[2] || nz >= xyz_max[2]) {
							continue;
						}

						int adj_id = ret.sampleMapId(nx, ny, nz);
						if (jumpOverList[adj_id])continue;

						v_insert(t.adjTile[d], adj_id);
					}

					ret.insert(t);
					ret.score += ret.calc_tile_score(t.id);
					ret.max_scales.insert({ ret.max_scale_of_covered(t.id), t.id });
				}
			}
		}

		ret.adapter_cnt = 0;
		ret.max_scale_limitation = 1.;
		ret.fabric_limitation = fabricShape[0] * fabricShape[1];
		ret.IDSPOOR.push_back(xSize * ySize * zSize);	//A sentry
		ret.IDSPOOR.push_back(xSize * ySize * zSize + 1);	//A sentry

#ifdef LOCAL
		printf("get_init_solution OK.\n");
#endif // LOCAL
		return ret;
	}

	Solution find_feasible_solution(int iss, double scale_limitation = 2.01, double aimSocre = 1.0) {

		Solution ret = get_init_solution(iss);

		double l = 1., r = scale_limitation;
		for (double l = 1.; l <= r; l += 0.1) {
			ret.adjustment(l, aimSocre);
			if (ret.total_tile_used() <= ret.fabric_limitation) {
				return ret;
			}
		}

		return Solution();
	}

	Solution find_suitable_solution(int iss, double scale, double aimSocre = 1.0) {
		int total_pes = fabricShape[0] * fabricShape[1];

		Solution ret = get_init_solution(iss);

		ret.adjustment(scale, aimSocre);

		return ret;
	}

	int estimated_suitable_sample_step() {
		int cnt[4] = { 0, 0, 0, 0 };
		for (int i = 0; i < (int)heatMap.size(); ++i) {
			if (heatMap[i] < DOUBLE_EPSILON)continue;
			if (heatMap[i] > 0.5)cnt[0]++;
			else if (heatMap[i] > 0.25)cnt[1]++;
			else if (heatMap[i] > 0.125)cnt[2]++;
			else cnt[3]++;
		}

		//db2(cnt[0], cnt[1]);
		//db2(cnt[2], cnt[3]);
		

		int costPEs = cnt[0] / 8 + cnt[1] / 64 + cnt[2] / 512 + cnt[3] / 4096;

		//db(costPEs);
		//db(pow(1. * fabricShape[0] * fabricShape[1] / costPEs, 1. / 3) * 10);
		//db(int(pow(1. * fabricShape[0] * fabricShape[1] / costPEs, 1. / 3) * 10));

		//exit(-2);

		return int(pow(1. * fabricShape[0] * fabricShape[1] / costPEs, 1. / 3) * 10);
	}

	void solve(const int &_dim,
		const vector<int> &_cubeShape,
		const vector<int> &_fabricShape,
		const double &_alpha, const double &_beta,
		const vector<double> &_heatMap, //inputs

		int &inv_sampling_step,
		vector<vector<int>>&pris) { //returns

		dimension = _dim;
		alpha = _alpha;
		beta = _beta;
		cubeShape = _cubeShape;
		heatMap = _heatMap;
		fabricShape = _fabricShape;

		int e_iss = estimated_suitable_sample_step();
		vector<int> p_iss = get_possible_sample_steps();

		int start_id = lower_bound(p_iss.begin(), p_iss.end(), e_iss) - p_iss.begin() - 1;

		//db3(e_iss, start_id, p_iss[start_id]);

		Solution bestSol = find_suitable_solution(p_iss[start_id], 2.0 + DOUBLE_EPSILON);

#ifdef LOCAL
		bestSol.print_tiles_used_info();
#endif // LOCAL

		inv_sampling_step = bestSol.inverse_sample_step;
		bestSol.format_to_prisms(pris);

		return;
	}
}
