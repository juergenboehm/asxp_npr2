#ifndef __triangularize_h
#define __triangularize_h

#include <QProgressBar>

// C2t3 Mesh
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>


#include <CGAL/IO/output_surface_facets_to_polyhedron.h>

//Polyhedron stuff

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>


//Mesh simplification stuff

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>


// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::Vector_3 Vector_3;
typedef GT::FT FT;

typedef FT (*Function)(Point_3);

typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

typedef typename C2t3::Triangulation Tr1;

typedef typename C2t3::Vertex_handle Vertex_handle;
typedef typename C2t3::Facet_iterator Facet_iterator;

typedef typename Tr1::Point Point;

//typedef CGAL::Polyhedron_3<GT> Polyhedron;
//typedef Polyhedron::HalfedgeDS HalfedgeDS;


extern Tr* trx;            // 3D-Delaunay triangulation
extern C2t3* c2t3;   // 2D-complex in 3D-Delaunay triangulation

//extern Polyhedron* c2t3_polyhedron;

extern std::map< Vertex_handle, int > verttoi_map;
extern std::vector< Vector_3 > normals_verti_vec;

extern int taubin_steps;
extern int manifold_sel;

FT sphere_function (Point_3 p);
Vector_3 sphere_normal(Point_3 p);


void remove_complexes();
void init_complexes();

int triangularize(double normal_thresh);

Vector_3 & getNormal(int i);
int getVerttoi(const Vertex_handle & v);


template <class C>
class EquivalenceClasses {
public:

	EquivalenceClasses() : global_cnt(0) {};
	void insert(C & v1, C & v2);
	int num_classes();

	std::map<C, int> coloring_map;
	int global_cnt;
	typedef typename std::map<C, int>::iterator Colmap_iterator;
};

extern QProgressBar *mainStatusProgressBar;


class TriSurface {

public:

	TriSurface(): newvert_index(0), newface_index(0) {};

	class Vertex {
	public:
		Point_3 p;
		Vector_3 n;
		int status;
	};

	typedef std::vector<Vertex >::iterator VertexIt;

	class Face {
	public:
		int vh[3];
		double aux;
	};

	typedef std::vector<Face >::iterator FaceIt;

	void clear();
	int addVertex(TriSurface::Vertex newv);
	int addFace(int v1, int v2, int v3);

	int newvert_index;
	int newface_index;

	std::vector<Vertex > vertex_vec;
	std::vector<Face > face_vec;


};

typedef int TriVertexHandle;
typedef TriSurface::Vertex TriVertex;
typedef TriSurface::Face TriFace;

extern TriSurface trisurf;

void compute_normals_trisurf(TriSurface & trisurf);


#endif
