
#include <iostream>

#include <QOpenGLFunctions>

#include <GL/glu.h>
#include "triangularize.h"


#include "glwidget.h"



ScreenGLWidget::ScreenGLWidget(QWidget *parent) : QGLWidget(parent)
{
	euler_phi = 0;
	euler_theta = 0;
	euler_psi = 0;
}


void ScreenGLWidget::rotate()
{
	repaint();
}


GLfloat n[6][3] = {  /* Normals for the 6 faces of a cube. */
  {-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0},
  {0.0, -1.0, 0.0}, /* changed z sign  */ {0.0, 0.0, -1.0}, {0.0, 0.0, 1.0} };
GLint faces[6][4] = {  /* Vertex indices for the 6 faces of a cube. */
  {0, 1, 2, 3}, {3, 2, 6, 7}, {7, 6, 5, 4},
  {4, 5, 1, 0}, {5, 6, 2, 1}, {7, 4, 0, 3} };
GLfloat v[8][3];  /* Will be filled in with X,Y,Z vertexes. */


void ScreenGLWidget::initializeGL()
{
	// Set up the rendering context, load shaders and other resources, etc.:
	QOpenGLFunctions *f = QOpenGLContext::currentContext()->functions();
	//f->glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

	/* Setup cube vertex data. */
	v[0][0] = v[1][0] = v[2][0] = v[3][0] = -1;
	v[4][0] = v[5][0] = v[6][0] = v[7][0] = 1;
	v[0][1] = v[1][1] = v[4][1] = v[5][1] = -1;
	v[2][1] = v[3][1] = v[6][1] = v[7][1] = 1;
	v[0][2] = v[3][2] = v[4][2] = v[7][2] = 1;
	v[1][2] = v[2][2] = v[5][2] = v[6][2] = -1;


	GLfloat light_diffuse[] = {1.0, 1.0, 1.0, 1.0};  /* White diffuse light. */
	GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};  /* Infinite light location. */


	/* Enable a single OpenGL light. */
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glEnable(GL_LIGHT0);

	GLfloat light_diffuse1[] = {1.0, 1.0, 1.0, 1.0};  /* White diffuse light. */
	GLfloat light_position1[] = {-1.0, -1.0, -1.0, 0.0};  /* Infinite light location. */

	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_POSITION, light_position1);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);


	GLfloat ambientfront[] = { 0.2, 0.0, 0.0, 1.0 };
	GLfloat diffusefront[] = { 0.8, 0.0, 0.0, 1.0 };

	GLfloat ambientback[] = { 0.0, 0.2, 0.0, 1.0 };
	GLfloat diffuseback[] = { 0.0, 0.8, 0.0, 1.0 };

	GLfloat specular [] = { 1.0, 1.0, 1.0, 1.0 };


	glMaterialfv ( GL_FRONT, GL_AMBIENT, ambientfront );
	glMaterialfv ( GL_FRONT, GL_DIFFUSE, diffusefront );
	glMaterialfv ( GL_FRONT, GL_SPECULAR, specular );


	glMaterialfv ( GL_BACK, GL_AMBIENT, ambientback );
	glMaterialfv ( GL_BACK, GL_DIFFUSE, diffuseback );
	glMaterialfv ( GL_BACK, GL_SPECULAR, specular );



	glMaterialf ( GL_FRONT_AND_BACK, GL_SHININESS, 100.0 );

	glLightModeli ( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE );


	//glShadeModel(GL_FLAT);
	glShadeModel(GL_SMOOTH);

	/* Use depth buffering for hidden surface elimination. */
	glEnable(GL_DEPTH_TEST);

	/* Setup the view of the cube. */
	glMatrixMode(GL_PROJECTION);

	gluPerspective( /* field of view in degree */ 40.0,
	/* aspect ratio */ 1.0,
	/* Z near */ 1.0, /* Z far */ 100.0);

#if 0
	glMatrixMode(GL_MODELVIEW);

	gluLookAt(0.0, 0.0, 5.0,  /* eye is at (0,0,5) */
	0.0, 0.0, 0.0,      /* center is at (0,0,0) */
	0.0, 1.0, 0.);      /* up is in positive Y direction */

	/* Adjust cube position to be asthetic angle. */
	glTranslatef(0.0, 0.0, -1.0);
	glRotatef(60, 1.0, 0.0, 0.0);
	glRotatef(-20, 0.0, 0.0, 1.0);
#endif


}

void ScreenGLWidget::resizeGL(int w, int h)
{

	glViewport(0, 0, (GLint)w, (GLint)h);

	// Update projection matrix and other size related settings:
	//m_projection.setToIdentity();
	//m_projection.perspective(45.0f, w / float(h), 0.01f, 100.0f);
}

#if 0

static void crossProduct(double* a, double* b, double* c)
{

	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];

}

static void normalize(double *a, double mul)
{
	double s = 0;
	for(int i = 0; i < 3; ++i) {
		s += a[i] * a[i];
	}
	s = sqrt(s);
	if (s != 0) {
		for(int i = 0; i < 3; ++i) {
			a[i] /= s;
			a[i] *= mul;
		}
	}
}
#endif


void ScreenGLWidget::paintGL()
{

	//std::cout << "do paintGL()" << std::endl;

	// Draw the scene:
	// QOpenGLFunctions *f = QOpenGLContext::currentContext()->functions();

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();

	gluLookAt(0.0, 0.0, 25.0,  /* eye is at (0,0,25) */
	0.0, 0.0, 0.0,      /* center is at (0,0,0) */
	0.0, 1.0, 0.);      /* up is in positive Y direction */

	/* Adjust cube position to be aesthetic angle. */
	glTranslatef(0.0, 0.0, -1.0);

#define TO_DEG(x) ((x)/3.14159265 * 180)

	//std::cout << "GL: euler_phi = " << TO_DEG(euler_phi) << " euler_theta = " << TO_DEG(euler_theta) << std::endl;

	glRotatef(TO_DEG(euler_phi), 1.0, 0.0, 0.0);
	glRotatef(TO_DEG(euler_theta), 0.0, 0.0, 1.0);

	glBegin(GL_TRIANGLES);

	// now do C2t3 facets

	if (c2t3 == 0 || trx == 0) {
		return;
	}

	int cnt1 = c2t3->number_of_facets();
	int cnt2 = 0;

	std::cout << "paintGL: nfaces = " << trisurf.face_vec.size() << std::endl;

	for(int iface = 0; iface < trisurf.face_vec.size(); ++iface) {
		TriSurface::Face & face = trisurf.face_vec[iface];
		TriSurface::Vertex & v1 = trisurf.vertex_vec[face.vh[0]];
		TriSurface::Vertex & v2 = trisurf.vertex_vec[face.vh[1]];
		TriSurface::Vertex & v3 = trisurf.vertex_vec[face.vh[2]];


		Point_3 & p1 = v1.p;
		Point_3 & p2 = v2.p;
		Point_3 & p3 = v3.p;

		Vector_3 & n1 = v1.n;
		Vector_3 & n2 = v2.n;
		Vector_3 & n3 = v3.n;

		glNormal3d(n1.x(), n1.y(), n1.z());
		glVertex3d(p1.x(), p1.y(), p1.z());

		glNormal3d(n2.x(), n2.y(), n2.z());
		glVertex3d(p2.x(), p2.y(), p2.z());

		glNormal3d(n3.x(), n3.y(), n3.z());
		glVertex3d(p3.x(), p3.y(), p3.z());

	}

#if TRY_OLD_C2T3
	for (C2t3::Facet_iterator fit = c2t3->facets_begin(); fit != c2t3->facets_end(); ++fit) {

	    const Tr::Cell_handle& cell = fit->first;
	    const int index = fit->second;

		// points on the facet

		const Vertex_handle & v1 = cell->vertex(trx->vertex_triple_index(index, 0));
		const Vertex_handle & v2 = cell->vertex(trx->vertex_triple_index(index, 1));
		const Vertex_handle & v3 = cell->vertex(trx->vertex_triple_index(index, 2));

		int iv1 = getVerttoi(v1);
		int iv2 = getVerttoi(v2);
		int iv3 = getVerttoi(v3);

		Point_3 & p1 = v1->point();
		Point_3 & p2 = v2->point();
		Point_3 & p3 = v3->point();

		Vector_3 & n1 = getNormal(iv1);
		Vector_3 & n2 = getNormal(iv2);
		Vector_3 & n3 = getNormal(iv3);

		glNormal3d(n1.x(), n1.y(), n1.z());
		glVertex3d(p1.x(), p1.y(), p1.z());

		glNormal3d(n2.x(), n2.y(), n2.z());
		glVertex3d(p2.x(), p2.y(), p2.z());

		glNormal3d(n3.x(), n3.y(), n3.z());
		glVertex3d(p3.x(), p3.y(), p3.z());
	}
#endif

	glEnd();

//#define SHOW_NORMALS
#ifdef SHOW_NORMALS

	glBegin(GL_LINES);

	for (C2t3::Vertex_iterator vit = c2t3.vertices_begin(); vit != c2t3.vertices_end(); ++vit) {

		Point_3 & p1 = vit->point();

		Vector_3 n1 = getNormal(getVerttoi(vit));
		Point_3 q1 = p1 + n1;

		glVertex3d(p1.x(), p1.y(), p1.z());
		glVertex3d(q1.x(), q1.y(), q1.z());

	}


	glEnd();

#endif

}
