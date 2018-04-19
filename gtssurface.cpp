
#include <iostream>
#include <math.h>

#include "glib.h"
#include "gts.h"


#include "triangularize.h"
#include "gtssurface.h"

typedef enum { NUMBER, COST } StopOptions;
typedef enum { COST_LENGTH, COST_OPTIMIZED, COST_ANGLE } CostOptions;
typedef enum { MIDVERTEX, OPTIMIZED } MidvertexOptions;

static gboolean stop_number_verbose (gdouble cost, guint number, guint * min)
{
  static guint nmax = 0, nold = 0;
  static GTimer * timer = NULL, * total_timer = NULL;

  g_return_val_if_fail (min != NULL, TRUE);

  if (timer == NULL) {
    nmax = nold = number;
    timer = g_timer_new ();
    total_timer = g_timer_new ();
    g_timer_start (total_timer);
  }

  if (number != nold && number % 121 == 0 &&
      number < nmax && nmax > *min) {
    gdouble total_elapsed = g_timer_elapsed (total_timer, NULL);
    gdouble remaining;
    gdouble hours, mins, secs;
    gdouble hours1, mins1, secs1;

    g_timer_stop (timer);

    hours = floor (total_elapsed/3600.);
    mins = floor ((total_elapsed - 3600.*hours)/60.);
    secs = floor (total_elapsed - 3600.*hours - 60.*mins);

    remaining = total_elapsed*((nmax - *min)/(gdouble) (nmax - number) - 1.);
    hours1 = floor (remaining/3600.);
    mins1 = floor ((remaining - 3600.*hours1)/60.);
    secs1 = floor (remaining - 3600.*hours1 - 60.*mins1);

    fprintf (stderr,
	     "\rEdges: %10u %3.0f%% %6.0f edges/s "
	     "Elapsed: %02.0f:%02.0f:%02.0f "
	     "Remaining: %02.0f:%02.0f:%02.0f ",
	     number,
	     100.*(nmax - number)/(nmax - *min),
	     (nold - number)/g_timer_elapsed (timer, NULL),
	     hours, mins, secs,
	     hours1, mins1, secs1);
    fflush (stderr);

    nold = number;
    g_timer_start (timer);
  }
  if (number < *min) {
    g_timer_destroy (timer);
    g_timer_destroy (total_timer);
    timer = NULL;
    total_timer = NULL;
    return TRUE;
  }
  return FALSE;
}

static gboolean stop_cost_verbose (gdouble cost, guint number, gdouble * max)
{
  g_return_val_if_fail (max != NULL, TRUE);

  if (number % 121 == 0) {
    fprintf (stderr, "\rEdges: %10u Cost: %10g ", number, cost);
    fflush (stderr);
  }
  if (cost > *max)
    return TRUE;
  return FALSE;
}

static gboolean stop_log_cost (gdouble cost, guint number)
{
  fprintf (stderr, "%d %g\n", number, cost);
  return FALSE;
}


void to_gtssurface(TriSurface & trisurf, GtsSurface* & s)
{

	s = gts_surface_new (gts_surface_class (),
			       gts_face_class (),
			       gts_edge_class (),
			       gts_vertex_class ());

	std::vector< GtsVertex* > gtsVertVec;

	for(int j = 0; j < trisurf.vertex_vec.size(); ++j) {
		TriSurface::Vertex & vv = trisurf.vertex_vec[j];
		gdouble x = vv.p.x();
		gdouble y = vv.p.y();
		gdouble z = vv.p.z();

		GtsVertex * v = gts_vertex_new(s->vertex_class, x, y, z);
		gtsVertVec.push_back(v);
	}

	for (int j = 0; j < trisurf.face_vec.size(); ++j) {

		TriSurface::Face & ff = trisurf.face_vec[j];
		int iv1 = ff.vh[0];
		int iv2 = ff.vh[1];
		int iv3 = ff.vh[2];

		GtsVertex* vv1 = gtsVertVec[iv1];
		GtsVertex* vv2 = gtsVertVec[iv2];
		GtsVertex* vv3 = gtsVertVec[iv3];


		GtsEdge * e1 = GTS_EDGE (gts_vertices_are_connected (vv1, vv2));
		GtsEdge * e2 = GTS_EDGE (gts_vertices_are_connected (vv2, vv3));
		GtsEdge * e3 = GTS_EDGE (gts_vertices_are_connected (vv3, vv1));

		if (e1 == NULL) {
		  e1 = gts_edge_new (s->edge_class, vv1, vv2);
		}

		if (e2 == NULL) {
		  e2 = gts_edge_new (s->edge_class, vv2, vv3);
		}

		if (e3 == NULL) {
		  e3 = gts_edge_new (s->edge_class, vv3, vv1);
		}
		gts_surface_add_face (s, gts_face_new (s->face_class, e1, e2, e3));
	}
	gts_surface_print_stats (s, stdout);

}

static void build_list (gpointer data, GSList ** list)
{
  /* always use O(1) g_slist_prepend instead of O(n) g_slist_append */
  *list = g_slist_prepend (*list, data);
}


void to_trisurface(GtsSurface* s, TriSurface & trisurf)
{
	   GSList *vertices=NULL;
	   GSList *triangles=NULL;
	   GSList *i=NULL;

	   GtsVertex *v1,*v2,*v3;
	   gint i1,i2,i3;

	   GtsPoint *p;

	   trisurf.clear();

	   gts_surface_foreach_vertex (s, (GtsFunc) build_list, &vertices);

	   guint vert_cnt = g_slist_length(vertices);

	   guint cnt = 0;

	   i=vertices;

	   std::map<GtsVertex*, int> vertex_map;

	   while(i){
	      p=GTS_POINT(i->data);
	      //fprintf(stdout,"%f %f %f\n",p->x,p->y,p->z);

	      vertex_map[GTS_VERTEX(i->data)] = cnt;

	      TriSurface::Vertex vv;
	      vv.p = Point_3(p->x, p->y, p->z);
	      vv.n = Vector_3(0., 0., 0.);
	      trisurf.vertex_vec.push_back(vv);

	      if (cnt % 500 == 0) {
	    	  fprintf(stdout, "\rvertcnt = %d               ", cnt);
	      }

	      //mainStatusProgressBar->setValue(int(50 + 20 * double(cnt)/vert_cnt));

	      i=i->next;
	      ++cnt;
	   }

	   gts_surface_foreach_face (s, (GtsFunc) build_list, &triangles);

	   guint tri_cnt = g_slist_length(triangles);
	   cnt = 0;

	   i=triangles;
	   while(i){
	      gts_triangle_vertices(GTS_TRIANGLE(i->data),&v1,&v2,&v3);
	      //i1=g_slist_index(vertices,v1); //+1;
	      //i2=g_slist_index(vertices,v2); //+1;
	      //i3=g_slist_index(vertices,v3); //+1;

	      i1 = vertex_map[v1];
	      i2 = vertex_map[v2];
	      i3 = vertex_map[v3];

	      //fprintf(stdout,"%d %d %d\n",i1,i2,i3);
	      TriSurface::Face ff;
	      ff.vh[0] = i1;
	      ff.vh[1] = i2;
	      ff.vh[2] = i3;
	      trisurf.face_vec.push_back(ff);

	      if (cnt % 1000 == 0) {
	    	  fprintf(stdout, "\rtricnt = %d              ", cnt);
	      }

	      //mainStatusProgressBar->setValue(int(70 + 20 * double(cnt)/tri_cnt));


	      i=i->next;
	      ++cnt;
	   }

	   g_slist_free(vertices);
	   g_slist_free(triangles);


}

void reduce_surface(GtsSurface* & s, int numedge_goal)
{
	const gdouble PI = 3.14159265;

	GtsKeyFunc cost_func = NULL;
	GtsCoarsenFunc coarsen_func = NULL;
	GtsStopFunc stop_func = NULL;

	gpointer coarsen_data = NULL, cost_data = NULL;
	gpointer stop_data = NULL;

	gdouble fold = PI/180.;
	guint number = 0;

	number = numedge_goal;


	GtsVolumeOptimizedParams params = { 0.5, 0.5, 0. };


	cost_func = (GtsKeyFunc) gts_volume_optimized_cost;
    cost_data = &params;

    coarsen_func = (GtsCoarsenFunc) gts_volume_optimized_vertex;
    coarsen_data = &params;

	stop_func = (GtsStopFunc) stop_number_verbose;
    stop_data = &number;

	gts_surface_coarsen (s,
			 cost_func, cost_data,
			 coarsen_func, coarsen_data,
			 stop_func, stop_data, fold);

	std::cout << "reduced surface:" << std::endl;

	gts_surface_print_stats (s, stdout);


}

void get_surface_quality(GtsSurface *s1, GtsSurface* s, double & max_distance)
{
	GSList *vertices=NULL;
	GSList *i = NULL;
	GtsPoint *p;
	GtsBBox* bbox = NULL;



	GNode *tree_s1 = gts_bb_tree_surface (s1);

	gts_surface_foreach_vertex (s, (GtsFunc) build_list, &vertices);

	guint vert_cnt = g_slist_length(vertices);

	guint cnt = 0;

	i=vertices;

	gdouble max_dist = -DBL_MAX;

	while(i){

	  p=GTS_POINT(i->data);

	  gdouble new_dist = sqrt(gts_bb_tree_point_distance(tree_s1, p, (GtsBBoxDistFunc)gts_point_triangle_distance2, &bbox));

	  if (new_dist > max_dist) {
		  max_dist = new_dist;
	  }

	  if (cnt % 500 == 0) {
		  fprintf(stdout, "\rchecking vertex distance = %9d maximum = %f	              ", cnt, max_dist);
	  }

	  //mainStatusProgressBar->setValue(int(50 + 20 * double(cnt)/vert_cnt));

	  i=i->next;
	  ++cnt;
	}

	g_slist_free(vertices);
	gts_bb_tree_destroy(tree_s1, TRUE);


	fprintf(stdout, "\n\nmaximal distance s to s1 = %f\n\n", max_dist);

	max_distance = max_dist;




}

gint count_vertices(GtsVertex *v, gint* data) {
	*data = *data + 1;
	return *data;
}

void reduce_gts_surface(TriSurface & trisurf, double red_perc, double & max_dist)
{
	GtsSurface* s;
	GtsSurface* s1;

	std::cout << "reduce_gts start..." << std::endl;

	to_gtssurface(trisurf, s);

	int num_edges = 0;


	s1 = gts_surface_new (gts_surface_class (),
		       gts_face_class (),
		       gts_edge_class (),
		       gts_vertex_class ());

	gts_surface_copy(s1, s);


	gts_surface_foreach_edge(s, (GtsFunc)count_vertices, &num_edges);

	std::cout << "edges = " << num_edges << std::endl;

	std::cout << "start reduce_surface..." << std::endl;

	reduce_surface(s, (int)(red_perc * num_edges));

	get_surface_quality(s1, s, max_dist);

	to_trisurface(s, trisurf);

	gts_object_destroy(GTS_OBJECT(s1));
	gts_object_destroy(GTS_OBJECT(s));

	std::cout << "compute normals trisurf..." << std::endl;

	compute_normals_trisurf(trisurf);

	std::cout << "reduce_gts end." << std::endl;
}

