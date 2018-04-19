
#include <QString>

#include <iostream>

#include <fstream>

#include <stdlib.h>
#include <stdio.h>

#include "asxp.h"
#include "global_headers.h"

#include "openscad.h"
#include "gtssurface.h"

#include "triangularize.h"

namespace SMS = CGAL::Surface_mesh_simplification;

int taubin_steps;
int manifold_sel;

QProgressBar *mainStatusProgressBar;


#define CUT_RADIUS 7.0



static double normal_thresh_local;



TriSurface trisurf;

void TriSurface::clear()
{
	vertex_vec.clear();
	face_vec.clear();
	newvert_index = 0;
	newface_index = 0;
}

int TriSurface::addVertex(TriSurface::Vertex newv)
{
	vertex_vec.push_back(newv);
	int newvert_index_old = newvert_index;
	++newvert_index;
	return newvert_index_old;
}

int TriSurface::addFace(int v1, int v2, int v3)
{
	TriSurface::Face f;
	f.vh[0] = v1;
	f.vh[1] = v2;
	f.vh[2] = v3;

	face_vec.push_back(f);

	int newface_index_old = newface_index;
	++newface_index;
	return newface_index_old;
}

template< class C >
void EquivalenceClasses<C>::insert(C & v1, C & v2)
{
	if (coloring_map.count(v1) == 0) {
		if (coloring_map.count(v2) == 0) {
			coloring_map[v1] = global_cnt;
			coloring_map[v2] = global_cnt;
			++global_cnt;
		} else {
			coloring_map[v1] = coloring_map[v2];
		}
	} else if (coloring_map.count(v2) == 0) {
		coloring_map[v2] = coloring_map[v1];
	} else {
		int i1 = coloring_map[v1];
		int i2 = coloring_map[v2];
		for(Colmap_iterator it = coloring_map.begin(); it != coloring_map.end(); ++it){
			if (it->second == i1 || it->second == i2) {
				coloring_map[it->first] = global_cnt;
			}
		}
		++global_cnt;
	}
};

template< class C >
int EquivalenceClasses<C>::num_classes()
{
	int cnt = 0;
	int ntst = -1;
	Colmap_iterator it = coloring_map.begin();
	while (it != coloring_map.end()) {
		if (it->second != ntst) {
			ntst = it->second;
			++cnt;
		}
		++it;
	}
	return cnt;
};


static std::vector<int> outp3;
static std::vector<Vector_3> outn3;

std::map< int, TriVertexHandle > convert_vertex_map;

Tr* trx;            // 3D-Delaunay triangulation
C2t3* c2t3;   // 2D-complex in 3D-Delaunay triangulation


std::map< Vertex_handle, int > vtoi_map;
std::vector< Vector_3 > normals_verti_vec;

std::vector< Vector_3 > vnew_vec;
std::vector< int > incid_vec;

std::vector< int > status_code_vec;
std::vector< C2t3::Vertex_iterator> itov_vec;


//std::ofstream ofgenpoly("genpoly.log");


static int tri_vertex_cnt = 0;
static int orig_vertex_cnt = 0;

static bool is_in_sphere(double r2, Point_3 & p)
{
	return (p - CGAL::ORIGIN)*(p - CGAL::ORIGIN) <= r2;
}

static double f_sphere(double r2, Point_3 p)
{
	return p.x() * p.x() + p.y() * p.y() + p.z() * p.z() - r2;
}

// computed with maple
static void get_abc_quad(double f0, double f1, double f2, double & a, double & b, double & c)
{
	c = f0;
	b = -(3./2.)*f0+2*f1-(1./2.)*f2;
	a = (1./2.)*f0-f1+(1./2.)*f2;
}

//
// t1 >= t2 holds in result
static void solve_abc_quad(double a, double b, double c, double & t1, double & t2)
{
	if (a == 0) {
		t1 = -DBL_MAX;
		if (b != 0) {
			t2 = -c/b;
		} else {
			t2 = DBL_MAX;
		}
	} else {
		double p = b/a;
		double q = c/a;
		double root2 = p * p / 4.0 - q;
		if (root2 < 0) {
			t1 = -DBL_MAX;
			t2 = DBL_MAX;
			return;
		}
		t1 = -p/2 - sqrt(root2);
		t2 = -p/2 + sqrt(root2);
	}
}

std::map< std::pair<int, int >, int > new_vertex_map;

typedef std::pair<int, int> int_index_pair;

// p1 is in sphere, p2 out
static void find_intersect_sphere_1(double r2, int iv1, int iv2,
					Vector_3 & n1, Vector_3 & n2, int & iv3, Vector_3 & n3, TriSurface & trisurf)
{

	if (new_vertex_map.count(int_index_pair(iv1, iv2)) > 0) {
		iv3 = new_vertex_map[int_index_pair(iv1, iv2)];

		return;
	}

	Point_3 p1 = itov_vec[iv1]->point();
	Point_3 p2 = itov_vec[iv2]->point();
	Point_3 p3;

	double f0 = f_sphere(r2, p1);
	double f1 = f_sphere(r2, p2);
	double f2 = f_sphere(r2, p1 + 2*(p2 - p1));

	double a, b, c;

	get_abc_quad(f0, f1, f2, a, b, c);

	double t1, t2;

	solve_abc_quad(a, b, c, t1, t2);

	p3 = p1 + t2 * (p2 - p1);
	n3 = n1 + t2 * (n2 - n1);

	n3 = n3/sqrt(n3 * n3);

	// ofgenpoly << "extra vertex building" << std::endl;

	TriSurface::Vertex newv3;

	newv3.p = p3;
	newv3.n = n3;

	int newv3_index = trisurf.addVertex(newv3);

	iv3 = newv3_index;

	new_vertex_map[int_index_pair(iv1, iv2)] = newv3_index;
	new_vertex_map[int_index_pair(iv2, iv1)] = newv3_index;



}

std::map<int, int> transfer_index_map;


int transfer_vertex(int iv, TriSurface & trisurf) {
	if (transfer_index_map.count(iv) > 0) {
		return transfer_index_map[iv];
	}

	Point_3 p1 = itov_vec[iv]->point();

	TriSurface::Vertex newvert;

	newvert.p = p1;

	newvert.n = getNormal(iv);

	int newvert_index = trisurf.addVertex(newvert);

	transfer_index_map[iv] = newvert_index;

	return newvert_index;
}




// we output the points in the *half-open* interval [p(iv1), p(iv2))
static void cut_against_sphere(double r2, int iv1, int iv2, Vector_3 & n1, Vector_3 & n2,
			int & res_code, std::vector<int> & outp3a, std::vector<Vector_3> & outn3a, TriSurface & trisurf)
{

	Point_3 p1 = itov_vec[iv1]->point();
	Point_3 p2 = itov_vec[iv2]->point();

	res_code = -1;

	bool isisp1 = is_in_sphere(r2, p1);
	bool isisp2 = is_in_sphere(r2, p2);

	if (isisp1 && isisp2) {
		int newv1 = transfer_vertex(iv1, trisurf);

		outp3a.push_back(newv1);
		outn3a.push_back(n1);
		res_code = 0;
	}
	if (isisp1 && !isisp2) {

		int iv3;
		Vector_3 n3;

		int newv1 = transfer_vertex(iv1, trisurf);

		outp3a.push_back(newv1);
		outn3a.push_back(n1);

		find_intersect_sphere_1(r2, iv1, iv2, n1, n2, iv3, n3, trisurf);
		outp3a.push_back(iv3);
		outn3a.push_back(n3);
		res_code = 1;
	}
	if (!isisp1 && isisp2) {

		// as p(iv1) is not in sphere and we compute a *half-open* interval
		// only the middle point p(v3) is given to output

		int iv3;
		Vector_3 n3;
		find_intersect_sphere_1(r2, iv2, iv1, n2, n1, iv3, n3, trisurf);

		outp3a.push_back(iv3);
		outn3a.push_back(n3);
		res_code = 2;
	}

}

static int face_counter = 0;

// is_special is true if one of the vertices belongs to the singular chain
static void displayFace(int iv1, int iv2, int iv3,
		Vector_3 & n1, Vector_3 & n2, Vector_3 & n3, TriSurface & trisurf)
{

	//std::cout << "displayFace ( " << iv1 << ", " << iv2 << ", " << iv3 << ")" << std::endl;

	const double r2cut = CUT_RADIUS * CUT_RADIUS;

	outp3.clear();
	outn3.clear();

	int res_code[3];

	cut_against_sphere( r2cut, iv1, iv2, n1, n2, res_code[0], outp3, outn3, trisurf);
	cut_against_sphere( r2cut, iv2, iv3, n2, n3, res_code[1], outp3, outn3, trisurf);
	cut_against_sphere( r2cut, iv3, iv1, n3, n1, res_code[2], outp3, outn3, trisurf);

	if (outp3.size() > 2) {
		int ivv1 = outp3[0];

		for(size_t i = 1; i < outp3.size() - 1; ++i) {

			int ivv2 = outp3[i];
			int ivv3 = outp3[i+1];

			trisurf.addFace(ivv1, ivv2, ivv3);
			++face_counter;


		}
	}
}


void build_trisurface(C2t3 * c2t3, Tr * trx1, TriSurface & trisurf)
{

	int cnt = 0;
	new_vertex_map.clear();
	transfer_index_map.clear();
	trisurf.clear();

	//ofgenpoly << "start build_trisurface" << std::endl;

	for (C2t3::Facet_iterator fit = c2t3->facets_begin(); fit != c2t3->facets_end(); ++fit) {

	    const Tr::Cell_handle& cell = fit->first;
	    const int index = fit->second;

		// points on the facet

		const Vertex_handle & v1 = cell->vertex(trx1->vertex_triple_index(index, 0));
		const Vertex_handle & v2 = cell->vertex(trx1->vertex_triple_index(index, 1));
		const Vertex_handle & v3 = cell->vertex(trx1->vertex_triple_index(index, 2));

		int iv1 = getVerttoi(v1);
		int iv2 = getVerttoi(v2);
		int iv3 = getVerttoi(v3);

		Point_3 & p1 = v1->point();
		Point_3 & p2 = v2->point();
		Point_3 & p3 = v3->point();

		Vector_3 & n1 = getNormal(iv1);
		Vector_3 & n2 = getNormal(iv2);
		Vector_3 & n3 = getNormal(iv3);

		if (!CGAL::collinear(p1, p2, p3)) {
			Vector_3 nn = CGAL::unit_normal(p1, p2, p3);

			if (nn * n1 < 0) {
				displayFace(iv1, iv3, iv2, n1, n3, n2, trisurf);
			} else {
				displayFace(iv1, iv2, iv3, n1, n2, n3, trisurf);
			}
		}
		++cnt;
	}

	std::cout << "faces processed = " << cnt << std::endl;
	std::cout << "vertices = " << trisurf.newvert_index << std::endl;
	std::cout << "faces = " << trisurf.newface_index << std::endl;
}

static double sign_center = 1;

FT sphere_function (Point_3 p) {
  const FT x = p.x();
  const FT y = p.y();
  const FT z = p.z();

  //const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
  const FT a = 2;
  const FT b = 1;
  const FT c = 1.5;

  double resval;
  eval_poly_f3(x, y, z, resval);

  return sign_center * resval;
  //return x*x/a/a + y*y/b/b + z*z/c/c -1;
}

Vector_3 sphere_normal(Point_3 p) {
	const FT a = 2;
	const FT b = 1;
	const FT c = 1.5;
	Vector_3 nn = Vector_3(p.x()/a/a, p.y()/b/b, p.z()/c/c);
	return nn/sqrt(nn * nn);
}

int getVerttoi(const Vertex_handle & v)
{
	return vtoi_map[v];
}


Vector_3 & getNormal(int i)
{
	return normals_verti_vec[i];
}


void laplacian_smooth(C2t3* c2t3x, double lambda)
{
	//std::cout << "start laplacian_smooth..." << std::endl;

	for(C2t3::Edge_iterator eit = c2t3x->edges_begin(); eit != c2t3x->edges_end(); ++eit) {

		const Vertex_handle & v1 = eit->first->vertex(eit->second);
		const Vertex_handle & v2 = eit->first->vertex(eit->third);

		Point_3 p1 = v1->point();
		Point_3 p2 = v2->point();

		int iv1 = vtoi_map[v1];
		int iv2 = vtoi_map[v2];

		Vector_3 vec1 = Vector_3(p1.x(), p1.y(), p1.z());
		Vector_3 vec2 = Vector_3(p2.x(), p2.y(), p2.z());

		vnew_vec[iv1] = vnew_vec[iv1] + vec2;
		vnew_vec[iv2] = vnew_vec[iv2] + vec1;

		incid_vec[iv1] = incid_vec[iv1] + 1;
		incid_vec[iv2] = incid_vec[iv2] + 1;

	}

	std::ofstream oflapp("laplacian.dat");

	for(C2t3::Vertex_iterator vit = c2t3x->vertices_begin(); vit != c2t3x->vertices_end(); ++vit) {

		int iv = vtoi_map[vit];

		if (incid_vec[iv] != 0) {

			Point_3 pmid = vit->point();

			Vector_3 vDelta = 1/double(incid_vec[iv]) * vnew_vec[iv];
			Vector_3 vMid( pmid.x(), pmid.y(), pmid.z() );
			Vector_3 vv = (1-lambda) * vMid + lambda * vDelta;

			oflapp << "function(Point old) = " << sphere_function(pmid) << std::endl;
			oflapp << "function(Point new) = " << sphere_function(Point_3(vv.x(), vv.y(), vv.z())) << std::endl;

			vit->set_point(Point_3(vv.x(), vv.y(), vv.z()));
		}

	}

	for(size_t i = 0; i < vnew_vec.size(); ++i) {
		vnew_vec[i] = Vector_3(0,0,0);
		incid_vec[i] = 0;
	}
}

void total_smooth(C2t3* c2t3x)
{

  for(int i = 0; i < taubin_steps; ++i) {

	  if (i % 1 == 0) {
		  std::cout << "i = " << i << " start laplacian_smooth" << std::endl;
	  }

	  laplacian_smooth(c2t3x, 0.330);
	  laplacian_smooth(c2t3x, -0.331);
  }


}

void init_complexes()
{
	trx = 0;
	c2t3 = 0;
}

void remove_complexes()
{

	if (trx != 0) {
		delete trx;
		trx = 0;
	}

	if (c2t3 != 0) {
		delete c2t3;
		c2t3 = 0;
	}
}

void check_surface_distance()
{
	std::ofstream ofroots("distance_to_surface.dat");

	double abs_max = -DBL_MAX;

	int cnt = 0;

	double root_average = 0;

	for (C2t3::Vertex_iterator it = c2t3->vertices_begin();
			it != c2t3->vertices_end(); ++it) {

		int iv = vtoi_map[it];
		Vector_3 nn = getNormal(iv);

		Point_3 & p = it->point();

		int root_list_len;
		double root_list[21];

		root_list_poly_point_normal(p.x(), p.y(), p.z(), nn.x(), nn.y(), nn.z(),
				root_list, root_list_len);

		double min_abs = DBL_MAX;
		double fval = DBL_MAX;
		int jmin = -1;

		for (int j = 0; j < root_list_len; ++j) {

			double abs_root = fabs(root_list[j]);

			if (abs_root < min_abs) {
				min_abs = abs_root;
				jmin = j;
			}
		}

		if (jmin != -1) {
			Point_3 pp = p + root_list[jmin] * nn;
			fval = sphere_function(pp);
		}

		if (min_abs > abs_max) {
			abs_max = min_abs;
		}

		root_average += min_abs;

		ofroots << "root[" << cnt << "]  = " << min_abs << " val = " << fval
				<< std::endl;

		++cnt;
	}

	root_average /= cnt;

	ofroots << "max_all_root = " << abs_max << std::endl;
	ofroots << "root_average = " << root_average << std::endl;

}

class LogHisto {

	static const int histo_len = 64;
	static const int histo_offset = 32;

public:

	LogHisto();
	void enter_val(double x);
	void display_histo_list(std::ostream & os);

	int histo_list[histo_len];

};

LogHisto::LogHisto() {
	memset(histo_list, 0, sizeof(histo_list));

}

void LogHisto::enter_val(double x) {
	double log_x = log(x)/log(10);
	int log_x_int = floor(log_x);
	log_x_int = log_x_int + histo_offset;
	if (log_x_int < 0) {
		log_x_int = 0;
	} else if (log_x_int >= histo_len) {
		log_x_int = histo_len - 1;
	}
	histo_list[log_x_int] += 1;
}

void LogHisto::display_histo_list(std::ostream & os) {
	for(int i = 0; i < histo_len; ++i) {
		os << "histo_list[" << i - histo_offset << "] = " << histo_list[i] << std::endl;
	}
}


void check_surface_distance_trisurf(TriSurface & trisurf)
{
	std::ofstream ofroots("distance_to_surface.dat");

	double abs_max = -DBL_MAX;

	int cnt = 0;

	double root_average = 0;

	LogHisto normal_dist_histo;
	LogHisto fval_histo;

	for (int j = 0; j < trisurf.vertex_vec.size(); ++j) {

		Vector_3 nn = trisurf.vertex_vec[j].n;

		Point_3 & p = trisurf.vertex_vec[j].p;

		double fval = sphere_function(p);

		fval_histo.enter_val(fabs(fval));

		int root_list_len;
		double root_list[21];

		root_list_poly_point_normal(p.x(), p.y(), p.z(), nn.x(), nn.y(), nn.z(),
				root_list, root_list_len);

		double min_abs = DBL_MAX;
		fval = DBL_MAX;
		int jmin = -1;

		for (int j = 0; j < root_list_len; ++j) {

			double abs_root = fabs(root_list[j]);

			if (abs_root < min_abs) {
				min_abs = abs_root;
				jmin = j;
			}
		}

		if (jmin != -1) {
			Point_3 pp = p + root_list[jmin] * nn;
			fval = sphere_function(pp);
		}

		if (jmin != -1 && min_abs > abs_max) {
			abs_max = min_abs;
		}

		if (jmin != -1) {
			root_average += min_abs;

			normal_dist_histo.enter_val(min_abs);

			ofroots << "root[" << cnt << "]  = " << min_abs << " val = " << fval
					<< std::endl;
		} else {
			ofroots << "root*[" << cnt << "]  = " << min_abs << " val = " << fval
					<< std::endl;
		}

		++cnt;
	}

	root_average /= cnt;

	ofroots << "max_all_root = " << abs_max << std::endl;
	ofroots << "root_average = " << root_average << std::endl;

	ofroots << "normal_dist_histo" << std::endl;
	normal_dist_histo.display_histo_list(ofroots);

	ofroots << "fval_histo" << std::endl;
	fval_histo.display_histo_list(ofroots);

}


void compute_normals(C2t3* c2t3x)
{

	for (C2t3::Facet_iterator fit = c2t3x->facets_begin();
			fit != c2t3x->facets_end(); ++fit) {

		const Tr::Cell_handle& cell = fit->first;
		const int index = fit->second;

#if 0
		const Tr::Facet & facet_mirror = c2t3x.opposite_facet(*fit);
		const int index_mirror = facet_mirror.second;
		const Tr::Cell_handle & cell_mirror = facet_mirror.first;
#endif

		Vertex_handle v1 = cell->vertex(trx->vertex_triple_index(index, 0));
		Vertex_handle v2 = cell->vertex(trx->vertex_triple_index(index, 1));
		Vertex_handle v3 = cell->vertex(trx->vertex_triple_index(index, 2));

		int iv1 = vtoi_map[v1];
		int iv2 = vtoi_map[v2];
		int iv3 = vtoi_map[v3];

		if (!CGAL::collinear(v1->point(), v2->point(), v3->point())) {

			Vector_3 n = CGAL::unit_normal(v1->point(), v2->point(),
					v3->point());

#if 0
			Point_3 & pv = cell->vertex(index)->point();
			Point_3 & pv_mirror = cell_mirror->vertex(index_mirror)->point();

			double val_pv = sphere_function(pv);
			double val_pv_mirror = sphere_function(pv_mirror);
#endif

			//assert(val_pv * val_pv_mirror <= 0);

			Point_3 p = v1->point();

			double root_list[21];
			int root_list_len;

			root_list_poly_point_normal(p.x(), p.y(), p.z(), n.x(), n.y(),
					n.z(), root_list, root_list_len);

			double min_abs = DBL_MAX;
			int jmin = -1;

			for (int j = 0; j < root_list_len; ++j) {

				double abs_root = fabs(root_list[j]);

				if (abs_root < min_abs) {
					min_abs = abs_root;
					jmin = j;
				}
			}

			// == 0 should never occur, but one doesn't know for sure
			if (sphere_function(p) == 0) {
				if (root_list_len > 1) {
					int jnext = jmin < root_list_len - 1 ? jmin + 1: jmin - 1;
					double h = 0.5 * root_list[jnext];
					double fh = sphere_function(p + h * n);
					if ((h < 0 && fh > 0) || (h > 0 && fh < 0)) {
						n = -n;
					}
				} else {
					if (sphere_function(p + n) < 0) {
						n = -n;
					}
				}
			} else if (root_list[jmin] * sphere_function(p) > 0) {
				n = -n;
			}

			normals_verti_vec[iv1] = normals_verti_vec[iv1] + n;
			normals_verti_vec[iv2] = normals_verti_vec[iv2] + n;
			normals_verti_vec[iv3] = normals_verti_vec[iv3] + n;
		}
	}

	for (std::vector<Vector_3>::iterator it = normals_verti_vec.begin();
			it != normals_verti_vec.end(); ++it) {
		Vector_3 & nn = *it;
		if (nn * nn > 0) {
			nn = nn / sqrt(nn * nn);
		}
	}

}

void compute_normals_trisurf(TriSurface & trisurf)
{

	for (int j = 0; j < trisurf.face_vec.size(); ++j) {

		TriSurface::Face & ff = trisurf.face_vec[j];

		int iv1 = ff.vh[0];
		int iv2 = ff.vh[1];
		int iv3 = ff.vh[2];

		Point_3 & p1 = trisurf.vertex_vec[iv1].p;
		Point_3 & p2 = trisurf.vertex_vec[iv2].p;
		Point_3 & p3 = trisurf.vertex_vec[iv3].p;

		if (!CGAL::collinear(p1, p2, p3)) {

			Vector_3 n = CGAL::unit_normal(p1, p2, p3);

			//assert(val_pv * val_pv_mirror <= 0);

			Point_3 & p = p1;

			double root_list[21];
			int root_list_len;

			root_list_poly_point_normal(p.x(), p.y(), p.z(), n.x(), n.y(),
					n.z(), root_list, root_list_len);

			double min_abs = DBL_MAX;
			int jmin = -1;

			for (int j = 0; j < root_list_len; ++j) {

				double abs_root = fabs(root_list[j]);

				if (abs_root < min_abs) {
					min_abs = abs_root;
					jmin = j;
				}
			}

			// == 0 should never occur, but one doesn't know for sure
			if (sphere_function(p) == 0) {
				if (root_list_len > 1) {
					int jnext = jmin < root_list_len - 1 ? jmin + 1: jmin - 1;
					double h = 0.5 * root_list[jnext];
					double fh = sphere_function(p + h * n);
					if ((h < 0 && fh > 0) || (h > 0 && fh < 0)) {
						n = -n;
					}
				} else {
					if (sphere_function(p + n) < 0) {
						n = -n;
					}
				}
			} else if (root_list[jmin] * sphere_function(p) > 0) {
				n = -n;
			}

			trisurf.vertex_vec[iv1].n = trisurf.vertex_vec[iv1].n + n;
			trisurf.vertex_vec[iv2].n = trisurf.vertex_vec[iv2].n + n;
			trisurf.vertex_vec[iv3].n = trisurf.vertex_vec[iv3].n + n;
		}
	}

	for (int j = 0; j < trisurf.vertex_vec.size(); ++j) {
		Vector_3 & nn = trisurf.vertex_vec[j].n;
		if (nn * nn > 0) {
			nn = nn / sqrt(nn * nn);
		}
		//std::cout << "nn = " << nn << " |nn|^2 = " << nn * nn << std::endl;
	}

}


void init_c2t3_aux_data(C2t3* c2t3x)
{
	//std::ofstream ofvv("vertex_values.dat");

	int cnt = 0;
	for (C2t3::Vertex_iterator vit = c2t3x->vertices_begin();
			vit != c2t3x->vertices_end(); ++vit) {

		Point_3 & pit = vit->point();
		double val = sphere_function(pit);
		//ofvv << "cnt = " << cnt << " val = " << val << std::endl;

		vtoi_map[vit] = cnt;
		++cnt;
		normals_verti_vec.push_back(Vector_3(0, 0, 0));
		vnew_vec.push_back(Vector_3(0, 0, 0));
		incid_vec.push_back(0);
		itov_vec.push_back(vit);
	}

	orig_vertex_cnt = cnt;

}

int cnt_pre_add_points;

void try_add_point(Tr* trx, double x1, double y1, double z1)
{
	double fval;
	double fnormal[3];
	double fhessian[9];

	eval_poly_poly_f3(x1, y1, z1, fval, fnormal, fhessian);

	double nnorm = fnormal[0] * fnormal[0] + fnormal[1] * fnormal[1] + fnormal[2] * fnormal[2];

	nnorm = sqrt(nnorm);

	if (nnorm < normal_thresh_local) {
		trx->insert(Point_3(x1, y1, z1));
		++cnt_pre_add_points;
	}

}
void prepare_initial_points(Tr* trx)
{
	const int root_list_len_max = 21;
	const int step = 1;

	cnt_pre_add_points = 0;

	for(int x = 0; x < 800; x += step ) {
		for(int y = 0; y < 800; y += step ) {

			double x1 = xrast_to_x(x);
			double y1 = yrast_to_y(y);

			for(int i = 0; i < (*pnsel_buf)[x][y]; ++i) {

				double z1 = (*pzfull_buf)[x][y][i];

				try_add_point(trx, x1, y1, z1);
			}
		}
	}

	for(int z = 0; z < 800; z += step) {
		for(int x = 0; x < 800; x += step) {

			double z1 = xrast_to_x(z);
			double x1 = xrast_to_x(x);

			double root_list[root_list_len_max];
			int root_list_len;

			root_list_poly_point_normal(x1, 0, z1, 0, 1, 0, root_list, root_list_len);

			for(int i = 0; i < root_list_len; ++i) {
				double y1 = root_list[i];

				try_add_point(trx, x1, y1, z1);

			}

		}


	}

	for(int z = 0; z < 800; z += step) {
		for(int y = 0; y < 800; y += step) {

			double z1 = xrast_to_x(z);
			double y1 = xrast_to_x(y);

			double root_list[root_list_len_max];
			int root_list_len;

			root_list_poly_point_normal(0, y1, z1, 1, 0, 0, root_list, root_list_len);

			for(int i = 0; i < root_list_len; ++i) {
				double x1 = root_list[i];

				try_add_point(trx, x1, y1, z1);


			}

		}

	}

	std::cout << "points added = " << cnt_pre_add_points << std::endl;
}


void triangularize_surface(C2t3* c2t3, Tr* trx, TriSurface & trisurfx,
		Surface_3 & surface,
		CGAL::Surface_mesh_default_criteria_3<Tr> &  criteria, int nstart)
{

	// meshing surface

	prepare_initial_points(trx);

	switch (manifold_sel) {
	case 0:
		std::cout << "Manifold tag()" << std::endl;
		CGAL::make_surface_mesh(*c2t3, surface, criteria,
				CGAL::Manifold_tag(), nstart);
		break;
	case 1:
		std::cout << "Manifold_with_boundary tag()" << std::endl;
		CGAL::make_surface_mesh(*c2t3, surface, criteria,
				CGAL::Manifold_with_boundary_tag(), nstart);
		break;
	case 2:
		std::cout << "Non_manifold tag()" << std::endl;
		CGAL::make_surface_mesh(*c2t3, surface, criteria,
				CGAL::Non_manifold_tag(), nstart);
		break;
	default:
		break;
	}

	mainStatusProgressBar->setValue(30);


	std::cout << "Final number of points in trx: " << trx->number_of_vertices() << "\n";

	//access vectors are initialized

	init_c2t3_aux_data(c2t3);

	total_smooth(c2t3);

	mainStatusProgressBar->setValue(40);


	std::cout << "compute normals" << std::endl;

	compute_normals(c2t3);

	//check_surface_distance();

	mainStatusProgressBar->setValue(50);


	std::cout << "build polyhedron" << std::endl;

	//CGAL::output_surface_facets_to_polyhedron(c2t3, c2t3_polyhedron);

	build_trisurface(c2t3, trx, trisurfx);

	//std::cout << "check surface distance trisurf..." << std::endl;

	//check_surface_distance_trisurf(trisurf);

	std::cout << "finished triangularize." << std::endl;

	mainStatusProgressBar->setValue(100);

}





int triangularize(double normal_thresh) {

	std::cout << "start triangularize..." << std::endl;

	normal_thresh_local = normal_thresh;

	std::cout << "normal_thresh = " << normal_thresh_local << std::endl;

	mainStatusProgressBar->setValue(0);

	Point_3 orig0(0, 0, 0);
	Point_3 orig = orig0;

	sign_center = 1;

	while (fabs(sphere_function(orig)) < 1e-9) {
		FT x = 0.1 * double(rand()) / RAND_MAX - 0.05;
		FT y = 0.1 * double(rand()) / RAND_MAX - 0.05;
		FT z = 0.1 * double(rand()) / RAND_MAX - 0.05;
		orig = orig0 + Vector_3(x, y, z);
	}

	double val_center = sphere_function(orig);

	if (val_center > 0) {
		sign_center = -1;
	}

	std::cout << "orig found..." << std::endl;
	mainStatusProgressBar->setValue(10);


	// defining the surface
	Surface_3 surface(sphere_function,             // pointer to function
			Sphere_3(orig, 64.0)); // bounding sphere
	// Note that "64." above is the *squared* radius of the bounding sphere!

#if TRY_HIRES

	const double resol = 0.02;

#else

	const double resol = 0.08;
	const int nstart = 800;

#endif


	// defining meshing criteria
	CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound, gold standard: 30.
			resol,  // radius bound
			resol); // distance bound

	remove_complexes();

	trx = new Tr();
	c2t3 = new C2t3(*trx);

	vtoi_map.clear();
	normals_verti_vec.clear();

	vnew_vec.clear();
	incid_vec.clear();

	itov_vec.clear();

	triangularize_surface(c2t3, trx, trisurf, surface, criteria, nstart);

	return 0;

}

