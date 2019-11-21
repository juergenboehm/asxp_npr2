
#ifndef __asxp_h
#define __asxp_h
 
#include <QBrush>
#include <QFont>
#include <QPen>
#include <QString>

#include "screenwidget.h"

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



    // the paint routines (with part routines factored out)

	void paint_compute_colmats(double a, double b, int xmax, int ymax);
	void paint_reset_full_plots();

	void paint_fullplot12(QPainter* painter);

    void paint_display_crossfield(int xmax, int ymax, QPainter* painter);

    void paint_from_colmat(QPainter* painter, int x, int y, int xmax, int ymax, double divider);

	void paint_silhouette_line(QPainter* painter, int xmax, int ymax);

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


private:

    int* pdata;

    //[gl_win_size * gl_win_size];

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




 
extern int phong_exponent;




#endif // __asxp_h
