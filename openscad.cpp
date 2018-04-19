
#include <iostream>
#include <fstream>
#include <algorithm>


#include "triangularize.h"
#include "openscad.h"

#define SCAL(x) (int((x) * 1e6)/1e6)


void render_triangle(std::ofstream & ofscad,
		int iv1, int iv2, int iv3, double shell_thickness)
{

	char buf[1024];
	const double h = 0.5 * shell_thickness;

	TriSurface::Vertex & v1 = trisurf.vertex_vec[iv1];
	TriSurface::Vertex & v2 = trisurf.vertex_vec[iv2];
	TriSurface::Vertex & v3 = trisurf.vertex_vec[iv3];


	Point_3 p1 = v1.p;
	Point_3 p2 = v2.p;
	Point_3 p3 = v3.p;

	if (CGAL::collinear(p1, p2, p3)) {
		return;
	}

	Point_3 u1 = p1 + (p2 - p1) * 0.5;
	Point_3 u2 = p2 + (p3 - p2) * 0.5;
	Point_3 u3 = p3 + (p1 - p3) * 0.5;

	Vector_3 n1 = v1.n;
	Vector_3 n2 = v2.n;
	Vector_3 n3 = v3.n;

	if (fabs(n1 * n1 - 1) > 1e-6) {
		std::cout << "n1 * n1 = " << n1 * n1 << std::endl;
	}
	if (fabs(n2 * n2 - 1) > 1e-6) {
		std::cout << "n2 * n2 = " << n2 * n2 << std::endl;
	}
	if (fabs(n3 * n3 - 1) > 1e-6) {
		std::cout << "n3 * n3 = " << n3 * n3 << std::endl;
	}


	Point_3 q1 = p1 - h * n1;
	Point_3 q2 = p2 - h * n2;
	Point_3 q3 = p3 - h * n3;

	Point_3 r1 = p1 + h * n1;
	Point_3 r2 = p2 + h * n2;
	Point_3 r3 = p3 + h * n3;

#if 1


	sprintf(buf, "%s(%.12f, %.12f, %.12f,  %.12f, %.12f, %.12f,  %.12f, %.12f, %.12f,"
			"  %.12f, %.12f, %.12f,  %.12f, %.12f, %.12f,  %.12f, %.12f, %.12f,"
			"  %.12f, %.12f, %.12f,  %.12f, %.12f, %.12f,  %.12f, %.12f, %.12f);",
			"rendertri",
			SCAL(q1.x()), SCAL(q1.y()), SCAL(q1.z()), SCAL(q2.x()), SCAL(q2.y()), SCAL(q2.z()),
			SCAL(q3.x()), SCAL(q3.y()), SCAL(q3.z()),
			SCAL(r1.x()), SCAL(r1.y()), SCAL(r1.z()), SCAL(r2.x()), SCAL(r2.y()), SCAL(r2.z()),
			SCAL(r3.x()), SCAL(r3.y()), SCAL(r3.z()),
			SCAL(u1.x()), SCAL(u1.y()), SCAL(u1.z()), SCAL(u2.x()), SCAL(u2.y()), SCAL(u2.z()),
			SCAL(u3.x()), SCAL(u3.y()), SCAL(u3.z()));
#endif

#if 0
	sprintf(buf, "%s(%d, %d, %d,  %d, %d, %d,  %d, %d, %d,"
			"  %d, %d, %d,  %d, %d, %d,  %d, %d, %d,"
			"  %d, %d, %d,  %d, %d, %d,  %d, %d, %d);",
			"rendertri",
			SCAL(q1.x()), SCAL(q1.y()), SCAL(q1.z()), SCAL(q2.x()), SCAL(q2.y()), SCAL(q2.z()),
			SCAL(q3.x()), SCAL(q3.y()), SCAL(q3.z()),
			SCAL(r1.x()), SCAL(r1.y()), SCAL(r1.z()), SCAL(r2.x()), SCAL(r2.y()), SCAL(r2.z()),
			SCAL(r3.x()), SCAL(r3.y()), SCAL(r3.z()),
			SCAL(u1.x()), SCAL(u1.y()), SCAL(u1.z()), SCAL(u2.x()), SCAL(u2.y()), SCAL(u2.z()),
			SCAL(u3.x()), SCAL(u3.y()), SCAL(u3.z()));
#endif

	ofscad << buf << std::endl;

}

bool fun_less(const TriSurface::Face & f1, const TriSurface::Face & f2)
{
	return f1.aux < f2.aux;
}

double lin_fun(Point_3 p, Vector_3 l) {
	Vector_3 vp(p.x(), p.y(), p.z());
	return vp * l;
}

void compute_face_aux(Vector_3 ll)
{

	for(int i = 0; i < trisurf.face_vec.size(); ++i) {
		TriSurface::Face & ff = trisurf.face_vec[i];
		Point_3 p1 = trisurf.vertex_vec[ff.vh[0]].p;
		Point_3 p2 = trisurf.vertex_vec[ff.vh[1]].p;
		Point_3 p3 = trisurf.vertex_vec[ff.vh[2]].p;


		double aux1 = lin_fun(p1, ll);
		double aux2 = lin_fun(p2, ll);
		double aux3 = lin_fun(p3, ll);

		double aux_res = std::max(std::max(aux1, aux2), aux3);

		ff.aux = aux_res;
	}
}

void prepare_trisurf()
{
	Vector_3 ll(double(rand())/RAND_MAX, double(rand())/RAND_MAX, double(rand())/RAND_MAX);
	ll = ll/sqrt(ll * ll);

	compute_face_aux(ll);

	std::sort(trisurf.face_vec.begin(), trisurf.face_vec.end(), fun_less);

}



void build_openscad(std::ofstream & ofscad, double shell_thickness)
{

	char buf[1024];

	setlocale(LC_NUMERIC, "en_US.UTF-8");

	prepare_trisurf();

	int cnt = 0;

#if 0
	ofscad << "module rendertri(xx1, yy1, zz1, xx2, yy2, zz2, xx3, yy3, zz3, x1, y1, z1, x2, y2, z2, x3, y3, z3,"
			"ux, uy, uz, vx, vy, vz, wx, wy, wz){" << std::endl;

	ofscad << "polyhedron([[xx1, yy1, zz1], [xx2, yy2, zz2], [xx3, yy3, zz3]," << std::endl;
	ofscad << "				[x1, y1, z1], [x2, y2, z2], [x3, y3, z3],"
							"[ux, uy, uz], [vx, vy, vz], [wx, wy, wz]]," << std::endl;
	ofscad << "[[0,1,2],[3,5,4], [6,0,3], [6,3,4], [6,4,1], [6,1,0],"
			"					 [7,1,4], [7,4,5], [7,5,2], [7,2,1],"
			"					 [8,2,5], [8,5,3], [8,3,0], [8,0,2]]);" << std::endl;
	ofscad << "}" << std::endl;
#endif

	ofscad << std::endl << std::endl;
	ofscad << "union(){" << std::endl;

	double aux_last = -DBL_MAX;
	double aux_akt;

	for(int iface = 0; iface < trisurf.face_vec.size(); ++iface) {

		TriSurface::Face & face = trisurf.face_vec[iface];
		int iv1 = face.vh[0];
		int iv2 = face.vh[1];
		int iv3 = face.vh[2];

		aux_akt = face.aux;

		render_triangle(ofscad, iv1, iv2, iv3, shell_thickness);

		if (cnt > 0 && cnt >= 100 && aux_akt > aux_last) {

			std::cout << "processing iface: " << iface << std::endl;

			sprintf(buf, "};\necho(\"processed %d\");\n union(){\n", cnt);
			ofscad << buf << std::endl;
			cnt = 0;
		} else {
			++cnt;
		}

		aux_last = aux_akt;
	}

	ofscad << "};" << std::endl;

}
