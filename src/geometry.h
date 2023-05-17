#pragma once
#include "header.h"

const int maxn = (int)1e6 + 10;
const double eps = 1e-10;
struct Point {
	double x, y;
	Point() {}
	Point(double _x, double _y) {
		x = _x; y = _y;
	}
};
typedef Point Vector;
Vector operator + (Vector A, Vector B) {
	return Vector(A.x + B.x, A.y + B.y);
}
Vector operator - (Vector A, Vector B) {
	return Vector(A.x - B.x, A.y - B.y);
}
Vector operator * (Vector A, double d) {
	return Vector(A.x*d, A.y*d);
}
Vector operator / (Vector A, double d) {
	return Vector(A.x / d, A.y / d);
}
double Dot(Vector A, Vector B) {
	return A.x*B.x + A.y*B.y;
}
double Cross(Vector A, Vector B) {
	return A.x*B.y - A.y*B.x;
}
double Length(Vector A) {
	return sqrt(Dot(A, A));
}
int dcmp(double x) {
	return fabs(x) < eps ? 0 : x < 0 ? -1 : 1;
}
Point P[maxn], Convexhull[maxn];
bool operator <(Point p1, Point p2) {
	return dcmp(p1.x - p2.x) < 0 || (dcmp(p1.x - p2.x) == 0 && dcmp(p1.y - p2.y) < 0);
}
int Andrew(int n) {      //Assume n>2
	sort(P, P + n);
	int top = 0;
	for (int i = 0; i < n; i++) {
		while (top > 1 && dcmp(Cross(P[i] - Convexhull[top - 2], Convexhull[top - 1] - Convexhull[top - 2])) >= 0)top--;
		Convexhull[top++] = P[i];
	}
	int k = top;
	for (int i = n - 2; i >= 0; i--) {
		while (top > k&&dcmp(Cross(P[i] - Convexhull[top - 2], Convexhull[top - 1] - Convexhull[top - 2])) >= 0)top--;
		Convexhull[top++] = P[i];
	}
	return top - 1;
}

double areaOfConvexhull(const vector<pair<int, int>>& pts) {
	int n = (int)pts.size();
	for (int i = 0; i < n; ++i) {
		P[i] = { 1. * pts[i].first, 1. * pts[i].second };
	}

	int m = Andrew(n);

	double area = 0.0;
	for (int i = 1; i < m - 1; ++i) {
		area += Cross(Convexhull[0] - Convexhull[i], Convexhull[0] - Convexhull[i + 1]);
		//area += 0.5 * Cross(Convexhull[0] - Convexhull[i], Convexhull[0] - Convexhull[i + 1]);
	}

	return fabs(area);
}

pair<pair<int, int>, pair<int, int>> diameterOfConvexhull(const vector<pair<int, int>>& pts) {
	int n = (int)pts.size();
	for (int i = 0; i < n; ++i) {
		P[i] = { 1. * pts[i].first, 1. * pts[i].second };
	}

	int m = Andrew(n);

	auto diaDist = [&](int l, int r) {
		return (Convexhull[l].x - Convexhull[r].x) * (Convexhull[l].x - Convexhull[r].x)
			 + (Convexhull[l].y - Convexhull[r].y) * (Convexhull[l].y - Convexhull[r].y);
	};

	pair<int, int>A, B;
	A = { (int)Convexhull[0].x, (int)Convexhull[0].y };
	B = { (int)Convexhull[1].x, (int)Convexhull[1].y };
	double maxDiameter = diaDist(0, 1);
	double currDiameter = maxDiameter;
	int p1 = 0, p2 = 1;

	while(p1 < m) {
		for (; diaDist(p1, (p2 + 1) % m) >= currDiameter;) {
			p2 = (p2 + 1) % m;
			currDiameter = diaDist(p1, p2);
		}
		if (currDiameter > maxDiameter) {
			maxDiameter = currDiameter;
			A = { (int)Convexhull[p1].x, (int)Convexhull[p1].y }; 
			B = { (int)Convexhull[p2].x, (int)Convexhull[p2].y };
		}
		currDiameter = diaDist(++p1, p2);
	}

	return { A, B };
}