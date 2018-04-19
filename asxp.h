
#ifndef __asxp_h
#define __asxp_h
 
#include <QBrush>
#include <QFont>
#include <QPen>
#include <QString>

#include "asxp_arrays.h"
#include "streamline.h"

 
class QPainter;
class QPaintEvent;

class QLabel;

class Streamplot;
 
class PaintHelper
 {
public:
    PaintHelper();
    ~PaintHelper();
 
public:

    void getImageInfo(int x, int y, double & a, double & b, double & c, double & d );

    void paint(QPainter *painter, QPaintEvent *event, bool mouse_akt,
    			int mousex, int mousey, double scale_im);
    void compute_colmat(double a, double b, int xmax, int ymax);

    void compute_colmat_gpu(double a, double b, int xmax, int ymax);

    void compute_streamfield(int xmax, int ymax);
    void compute_streamfield_CGAL(int xmax, int ymax);

    void init_colmat();


public:
    QBrush background;
    QBrush circleBrush;
    QFont textFont;
    QPen circlePen;
    QPen textPen;

    double euler_phi;
    double euler_theta;
    double euler_psi;

    double zsobel_thresh;
    double nsobel_thresh;

    int shade_type;

    bool displayFlowLines;
    bool displayCrossField;

    int streamLineColor;

    double dsep;

    int *colmat_r;
    int *colmat_g;
    int *colmat_b;

    Streamplot* full_plot1;
    Streamplot* full_plot2;

    bool colmat_valid;

    int streamgen_type; // 0 for own 1 for CGAL

    QLabel* displayLabel;
};

class FindSilhouette : public PointClassifier {

public:

	virtual ~FindSilhouette() {};

	virtual int operator()(double x, double y);
	void set_xy_max(int xmaxa, int ymaxa);

	int xmax;
	int ymax;
};

class FindBackground : public PointClassifier {

public:

	virtual ~FindBackground() {};

	virtual int operator()(double x, double y);
	void set_xy_max(int xmaxa, int ymaxa);

	int xmax;
	int ymax;
};


int init_f5(QString f5str);

void eval_poly_f3(double x, double y, double z, double & resval);

int eval_poly_poly_f3(double x, double y, double z, double & f, double* fnormal, double* fhessian);


int root_list_poly_point_normal(double x, double y, double z,
		double nx, double ny, double nz, double* root_list, int & root_list_len);


 
extern int phong_exponent;

extern Array2d_double* pz_buf;
extern Array2d_double* pn_buf;


extern Array3d_double* pzfull_buf;
extern Array2d_int* pjsel_buf;
extern Array2d_int* pnsel_buf;

extern Array2d_bool* pis_silhouette_mat;

extern Array3d_double* pshape_mat_buf;

extern Array3d_double* pvbase_1_buf;
extern Array3d_double* pvbase_2_buf;

extern Array3d_double* pv1_buf;
extern Array3d_double* pv2_buf;

extern Array2d_double* pl1_buf;
extern Array2d_double* pl2_buf;




#endif // __asxp_h
