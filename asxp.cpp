
#include <QtGui>
#include <QString>
#include <QLabel>

#include <list>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <timing.h>

#include "asxp_arrays.h"



#include "poly.h"
#include "parser.h"
#include "roots.h"

#include "screenwidget.h"

#include "procpoly.h"

#include "asxp.h"

#include "eigen.h"

#include "streamline.h"
#include "pointlist.h"
#include "streamplot.h"

#define NOCUDA

#include "pimage.h"


int phong_exponent;

Scale xrast_to_x, yrast_to_y, x_to_xrast, y_to_yrast;

 

#define Q_SQ(x) ((x) * (x))
#define Q_NORM(x,y) (Q_SQ(x) + Q_SQ(y))
#define Q_SKALPROD(x1, y1, x2, y2) ((x1)*(x2) + (y1)*(y2))

#define TO_GRAD(x) ((x)/(2 * 3.14159265) * 360)
#define SIGN(x) ((x) > 0 ? 1 : (((x) == 0) ? 0 : -1))

#define UNUSED(expr) do { (void)(expr); } while (0)



inline int & mref(int* mat, int i, int j)
{
	return mat[i * 2 * gl_win_size + j];
}


PaintHelper::PaintHelper()
 {
    QLinearGradient gradient(QPointF(50, -20), QPointF(80, 20));
    gradient.setColorAt(0.0, Qt::white);
    gradient.setColorAt(1.0, QColor(0xa6, 0xce, 0x39));
 
    background = QBrush(QColor(200, 200, 200));
    circleBrush = QBrush(gradient);
    circlePen = QPen(Qt::black);
    circlePen.setWidth(1);
    textPen = QPen(Qt::white);
    textFont.setPixelSize(50);

    const int colmat_size = 4 * gl_win_size * gl_win_size * sizeof(int);

    colmat_r = (int*)malloc(colmat_size);
    colmat_g = (int*)malloc(colmat_size);
    colmat_b = (int*)malloc(colmat_size);

    init_colmat();

    full_plot1 = 0;
    full_plot2 = 0;

    dsep = 10;

    displayFlowLines = false;
    displayCrossField = true;

    streamLineColor = 0;

    streamgen_type = 1;

    colmat_valid = true;

    pdata = new int[gl_win_size * gl_win_size];
}

void PaintHelper::init_colmat()
{
	for(int x = 0; x < gl_win_size; ++x) {
		for(int y = 0; y < gl_win_size; ++y) {
			mref(colmat_r, x, y) = 200; //64;
			mref(colmat_g, x, y) = 200;// 32;
			mref(colmat_b, x, y) = 200;  // 64;
		}
	}
}

PaintHelper::~PaintHelper()
{
	free(colmat_r);
	free(colmat_g);
	free(colmat_b);

	free(pdata);
}




void PaintHelper::getImageInfo(int x, int y, double & a, double & b, double & c, double & d )
{

	a = pp.l1_buf[x][y];
	b = pp.l2_buf[x][y];

	c = a * b;
	d = a + b;

	a = pp.nsel_buf[x][y];
	b = pp.jsel_buf[x][y];

}



//#define SCALE 10.0
#define SCALE 50.0

void calc_brightness( double nz, bool & is_backside, double & brightness, double & highlight) {

	double nabs = fabs(nz);

	is_backside = (nz < 0);

	brightness = nabs/2;

	double phong_kernel = 2 * nz * nz - 1;
	double spec_coef = 0.3 * pow(phong_kernel, phong_exponent);

	highlight = spec_coef;

}


void shade_color(double nz, double & color_red, double & color_green, double & color_blue ) {

	double brightness;
	double highlight;
	bool is_backside;

	calc_brightness(nz, is_backside, brightness, highlight);

	if (is_backside) {
		color_red = brightness + highlight;
		color_green = 0;
		color_blue = 0;

	} else {
		color_red = 0;
		color_green = brightness + highlight;
		color_blue = 0;
	}

}

const int win_size =  gl_win_size;

Array2d_double floyd_stein_buf(boost::extents[win_size + 1][win_size + 1]);

void shade_dots(double nz, double & color_red, double & color_green, double & color_blue ) {

	double brightness;
	double highlight;
	bool is_backside;

	calc_brightness(nz, is_backside, brightness, highlight);

	color_red = brightness + highlight;
	color_green = brightness + highlight;
	color_blue = brightness + highlight;

	double rand_thr = ((double (rand()))/(RAND_MAX));


	if (color_red > rand_thr) {
		color_red = 0.0;
		color_green = 0.0;
		color_blue = 0.0;
	} else {
		color_red = 1.0;
		color_green = 1.0;
		color_blue = 1.0;
	}


}

void shade_grey(double nz, double & c_grey)
{
	double brightness;
	double highlight;
	bool is_backside;

	calc_brightness(nz, is_backside, brightness, highlight);

	c_grey = brightness + highlight;

}

void shade_floyd_steinberg(int x, int y, double & color_red, double & color_green, double & color_blue)
{
	bool is_set = floyd_stein_buf[x][y] > 0.5;

	double error;

	if (is_set) {
		error = floyd_stein_buf[x][y] - 1;

		color_red = color_green = color_blue = 1;

	} else {

		error = floyd_stein_buf[x][y];

		color_red = color_green = color_blue = 0;
	}

	floyd_stein_buf[x + 1][y] = floyd_stein_buf[x + 1][y] + (3 * error) / 8;
	floyd_stein_buf[x][y + 1] = floyd_stein_buf[x][y + 1] + (3 * error) / 8;
	floyd_stein_buf[x + 1][y + 1] = floyd_stein_buf[x + 1][y + 1] + error / 4;

}



//  y
//  ^
//  |   a        b
//  |
//  |    ..sx *
//  |        sy
//  |         .
//  |   c     .  d
//  |
//  -----------------> x



static double bilin_interpolate(double a, double b, double c, double d, double sx, double sy) {

	double val_a = a * sy * (1 - sx);
	double val_b = b * sx * sy;
	double val_c = c * (1 - sx) * (1 - sy);
	double val_d = d * sx * (1 - sy);

	double val = val_a + val_b + val_c + val_d;

	return val;
}




void PaintHelper::compute_colmat(double a, double b, int xmax, int ymax)
{

	double local_scale = gl_win_size/xmax;

	double parm[2];
	int win_offset = xmax/2;
	bool disc_zero;

	UNUSED(win_offset);
	UNUSED(disc_zero);

	UNUSED(local_scale);

	int x;
	int y;

	cout << "xmax = " << xmax << "ymax = " << ymax << endl;


	parm[0] = a;
	parm[1] = b;

    pp.eval_parms_poly(pp.f5, pp.f3, parm);
    pp.rotate_mat(euler_phi, euler_theta, euler_psi);

	pp.init_f3_diff(pp.f3);

#ifdef RANDOM_COLMAT
	for (x = 0; x < ScreenWidget::gl_win_size; ++x) {
		for( y = 0; y < ScreenWidget::gl_win_size; ++y) {
			mref(colmat_r, x, y) = RAND_COL;
			mref(colmat_g, x, y) = RAND_COL;
			mref(colmat_b, x, y) = RAND_COL;
		}
	}

	return;
#endif

	pp.fill_arrays(xmax, ymax, shade_type);

	double k_gauss_min, k_gauss_max;

	if (!(shade_type == 1 || shade_type == 2)) {
		pp.calc_silhouette(xmax, ymax, k_gauss_min, k_gauss_max);
	}

	//cout << "k_gauss_min = " << k_gauss_min << " k_gauss_max = " << k_gauss_max << endl;

	if (shade_type == 2) {
		for(x = 0; x < xmax; ++x) {
			for(y = 0; y < ymax; ++y) {

				double n = pp.n_buf[x][y];
				double z = pp.z_buf[x][y];

				if (z > M_INF) {

					double c_grey; // between 0.0 and 1.0

					shade_grey(n, c_grey);

					floyd_stein_buf[x][y] = c_grey;
				} else
					floyd_stein_buf[x][y] = 1;
			}
		}
	}

    for (x = 0; x < xmax; ++x) {
    	for (y = 0; y < ymax; ++y) {

    		double z;
    		double n;

    		z = pp.z_buf[x][y];
    		n = pp.n_buf[x][y];

    		if (z > M_INF) {

    			double color_red = 1.0;
    			double color_green = 1.0;
    			double color_blue = 1.0;

    			bool is_silhouette = false;
    			bool is_true_silhouette = false;

    			if (!(shade_type == 1 || shade_type == 2)) {

    				is_silhouette = pp.is_silhouette_mat[x][y];

					is_true_silhouette = false;

					if (fabs(n) <= nsobel_thresh * 0.2) {
						is_silhouette = true;
						is_true_silhouette = true;
					}

    			}

    			if (is_silhouette == true) {

    				pp.v1_buf[x][y][0] = 0;
    				pp.v1_buf[x][y][1] = 0;

    				pp.v2_buf[x][y][0] = 0;
    				pp.v2_buf[x][y][1] = 0;

    				color_red = 0.0;
    				color_green = 0.0;
    				color_blue = 0.0;

    				if (is_true_silhouette) {
    					color_red = 1.0;
    				}

    			} else {

    				if (shade_type == 1) {
    					shade_color(n, color_red, color_green, color_blue);
    				}

    				if (shade_type == 2) {
    					shade_floyd_steinberg(x, y, color_red, color_green, color_blue);
    				}

    			}

    			int col_red = (int) (250 * color_red);
    			int col_green = (int) (250 * color_green);
    			int col_blue = (int) (250 * color_blue);

    			mref(colmat_r, x, y) = col_red;
    			mref(colmat_g, x, y) = col_green;
    			mref(colmat_b, x, y) = col_blue;

    		} else {

    			// background color
    			mref(colmat_r, x, y) = 250; // 64
    			mref(colmat_g, x, y) = 250; // 32
    			mref(colmat_b, x, y) = 250; // 64

    		}
    	}
    }


}



int translate_poly(cudaPoly3 & f3_h, Poly3 & f3)
{

	Poly3::vecIt it;
	int i;

	memset(f3_h.xexpo, 0, 40 * sizeof(int));
	memset(f3_h.yexpo, 0, 40 * sizeof(int));
	memset(f3_h.zexpo, 0, 40 * sizeof(int));


	for(i = 0, it = f3.terms.begin(); it < f3.terms.end(); ++it, ++i) {
		f3_h.coefs[i] = it->coef;
		f3_h.xexpo[i] = it->expos[0];
		f3_h.yexpo[i] = it->expos[1];
		f3_h.zexpo[i] = it->expos[2];
	}

	f3_h.len = i;
}


void PaintHelper::compute_colmat_gpu(double a, double b, int xmax, int ymax)
{

	double local_scale = gl_win_size/xmax;

	double parm[2];
	int win_offset = xmax/2;
	bool disc_zero;

	int x;
	int y;

	cout << "xmax = " << xmax << "ymax = " << ymax << endl;


	parm[0] = a;
	parm[1] = b;

    pp.eval_parms_poly(pp.f5, pp.f3, parm);

    cudaPoly3 f3_cuda_h;

    translate_poly(f3_cuda_h, pp.f3);

	gpu_compute_colmat(a, b, xmax, ymax, f3_cuda_h, euler_phi, euler_theta, euler_psi,
			colmat_r, colmat_g, colmat_b);

}



// (x, y) raster coordinates

void FindSilhouette::set_xy_max(int xmaxa, int ymaxa) {
	xmax = xmaxa;
	ymax = ymaxa;
}

int FindSilhouette::operator()(double x, double y) {

	int cx = int(x);
	int cy = int(y);

	if ((cx < 0 || cx >= xmax) || (cy < 0 || cy >= ymax)) {
		return true;
	} else {
		return pp.is_silhouette_mat[cx][cy];
	}
}

FindSilhouette is_silhouette_point;

void FindBackground::set_xy_max(int xmaxa, int ymaxa) {
	xmax = xmaxa;
	ymax = ymaxa;
}

int FindBackground::operator()(double x, double y) {

	int cx = int(x);
	int cy = int(y);

	if ((cx < 0 || cx >= xmax) || (cy < 0 || cy >= ymax)) {
		return true;
	} else {
		return (pp.z_buf[cx][cy] == M_INF);
	}
}

FindBackground is_background_point;


void PaintHelper::compute_streamfield_CGAL(int xmax, int ymax)
{

	const double dsep_multiplier = 4 * 0.6;

	if (!full_plot1) {
		full_plot1 = new Streamplot(&pp.v1_buf, &pp.v2_buf, xmax, ymax, dsep * dsep_multiplier,
								xrast_to_x, yrast_to_y, x_to_xrast, y_to_yrast,
								&is_background_point);

		full_plot1->set_mode(1);

		full_plot1->compute_stream_field_CGAL();
	}

	if (!full_plot2) {

		full_plot2 = new Streamplot(&pp.v1_buf, &pp.v2_buf, xmax, ymax, dsep * dsep_multiplier,
								xrast_to_x, yrast_to_y, x_to_xrast, y_to_yrast,
								&is_background_point);

		full_plot2->set_mode(2);

		full_plot2->compute_stream_field_CGAL();
	}

}

void PaintHelper::compute_streamfield(int xmax, int ymax)
{

	bool init_streamfield1 = true;
	bool init_streamfield2 = true;

	if (streamgen_type == 1) {
		compute_streamfield_CGAL(xmax, ymax);
		return;
	}

	for(int tstx = 0; tstx < xmax; tstx += 10) {

		for(int tsty = 0; tsty < ymax; tsty += 10) {

			char buf[128];

			sprintf(buf, "tstx = %d, tsty = %d", tstx, tsty);

			displayLabel->setText(QString(buf));
			displayLabel->repaint();

			cout << "tstx = " << tstx << " tsty = " << tsty << endl;

			while (1) {

				if (pp.z_buf[tstx][tsty] == M_INF) {
					break;
				}

				double xstart, ystart;

				xstart = xrast_to_x(tstx);
				ystart = yrast_to_y(tsty);

#if 0
				Streamline flow_line(&v1_buf, &v2_buf, xmax, ymax, x_to_xrast, y_to_yrast, &is_silhouette_point);

				flow_line.integrate_from(xstart, ystart);
#endif


				if (!full_plot1) {

					full_plot1 = new Streamplot(&pp.v1_buf, &pp.v2_buf, xmax, ymax, dsep,  xrast_to_x, yrast_to_y, x_to_xrast, y_to_yrast,
											&is_silhouette_point);

					full_plot1->set_mode(1);

					init_streamfield1 = true;

				}

				cout << "compute streamfield 1" << endl;

				full_plot1->compute_stream_field(xstart, ystart, init_streamfield1);

				init_streamfield1 = false;

				if (!full_plot2) {

					full_plot2 = new Streamplot(&pp.v1_buf, &pp.v2_buf, xmax, ymax, dsep,  xrast_to_x, yrast_to_y, x_to_xrast, y_to_yrast,
											&is_silhouette_point);

					full_plot2->set_mode(2);


					init_streamfield2 = true;
				}

				cout << "compute streamfield 2" << endl;

				full_plot2->compute_stream_field(xstart, ystart, init_streamfield2);

				init_streamfield2 = false;

				break;

			}
		}
	}



}



void PaintHelper::paint_fullplot12(QPainter* painter) {

	switch (streamLineColor) {
	case 0:
		painter->setPen(QPen(QColor(0, 0, 0)));
		break;
	case 1:
		painter->setPen(QPen(Qt::red));
		break;
	default:
		painter->setPen(QPen(Qt::blue));
		break;
	}
	// dummy while, only for break
	Streamplot* full_plot_arr[2];

	full_plot_arr[0] = full_plot1;
	full_plot_arr[1] = full_plot2;

	for (int ifpa = 0; ifpa < 2; ++ifpa) {

		Streamplot* full_plot = full_plot_arr[ifpa];

		if (full_plot == 0) {
			continue;
		}

		cout << "no of streamlines = " << full_plot->streaml_res.size() << endl;
		// cin.get();

		Point2DList fringe_d = full_plot->debug_points;

		for (Point2DListIterator pit = fringe_d.begin(); pit != fringe_d.end();
				++pit) {
			painter->drawPoint(pit->x, pit->y);
		}

		for (list<Streamline>::iterator psl = full_plot->streaml_res.begin();
				psl != full_plot->streaml_res.end(); ++psl) {

			Point2DList& line_samples(psl->line_samples);
			Point2DListIterator pit = line_samples.begin();

			cout << "Anfang size = " << line_samples.size() << endl;
			cout << endl << "Ende" << endl;

			// line_samples contains screen coordinates

			pit = line_samples.begin();

			double x0 = pit->x;
			double y0 = pit->y;

			++pit;

			while (pit != line_samples.end()) {

				double x1 = pit->x;
				double y1 = pit->y;

				++pit;

				painter->drawLine(x0, y0, x1, y1);

				x0 = x1;
				y0 = y1;
			}
		}

		cout << "no of streamlines = " << full_plot->streaml_res.size() << endl;
	}
}

void PaintHelper::paint_display_crossfield(int xmax, int ymax, QPainter* painter) {

	for (int x = 0; x < xmax; ++x) {
		for (int y = 0; y < ymax; ++y) {

			if (pp.z_buf[x][y] == M_INF) {
				continue;
			}

			if (x % 16 != 0 || y % 16 != 0) {
				continue;
			}

			const double c = 7;

			double v1x = c * pp.v1_buf[x][y][0];
			double v1y = c * pp.v1_buf[x][y][1];

			double v2x = c * pp.v2_buf[x][y][0];
			double v2y = c * pp.v2_buf[x][y][1];

			painter->setPen(QPen(Qt::blue));

			painter->drawLine(x, y, x + v1x, y + v1y);
			painter->drawLine(x, y, x - v1x, y - v1y);

			painter->setPen(QPen(Qt::green));

			painter->drawLine(x, y, x + v2x, y + v2y);
			painter->drawLine(x, y, x - v2x, y - v2y);

		}
	}
}

void PaintHelper::paint_from_colmat(QPainter* painter, int x, int y, int xmax, int ymax, double divider) {

	Timing time0;
	time0.start();

	int stride = gl_win_size / divider;

	for (x = 0; x < xmax; ++x) {

		for (y = 0; y < ymax; ++y) {

			int col_red = mref(colmat_r, x, y);
			int col_green = mref(colmat_g, x, y);
			int col_blue = mref(colmat_b, x, y);
			pdata[y * stride + x] = (255 << 24) + (col_red << 16)
					+ (col_green << 8) + col_blue;
		}
	}


	time0.stop();

	cout << "time to copy = " << time0.elapsed() << endl;

	int im_size = gl_win_size / divider;

	QImage qimage((uchar*) (pdata), im_size, im_size, QImage::Format_ARGB32);

	painter->drawImage(0, 0, qimage);
}

void PaintHelper::paint_silhouette_line(QPainter* painter, int xmax, int ymax) {

	painter->fillRect(0, 0, xmax, ymax, QBrush(QColor(255, 255, 255)));

	painter->setPen(QPen(QColor(0, 0, 0)));

	painter->setBrush(QBrush(QColor(0, 0, 0)));

	for (Point2DListIterator it = pp.silhouette_pointl.begin();
			it != pp.silhouette_pointl.end(); ++it) {

		painter->drawPoint(it->x, it->y);

	}
}

void PaintHelper::paint_reset_full_plots() {

	if (full_plot1 != 0) {
		delete full_plot1;
		full_plot1 = 0;
	}

	if (full_plot2 != 0) {
		delete full_plot2;
		full_plot2 = 0;
	}
}

void PaintHelper::paint_compute_colmats(double a, double b, int xmax, int ymax) {

	bool isPrint = (shade_type == -1);

	if (isPrint) {
		shade_type = 0;
	}

	Timing time0;
	time0.start();

	//compute_colmat(a, b, gl_win_size/divider, gl_win_size/divider);
	if (shade_type != 1) {
		compute_colmat(a, b, xmax, ymax);
	} else {
		compute_colmat_gpu(a, b, xmax, ymax);
	}

	time0.stop();

	double elapsed = time0.elapsed();
	cout << "time to render = " << elapsed << endl;

	colmat_valid = true;

	paint_reset_full_plots();

	if (isPrint) {
		compute_streamfield(xmax, ymax);
	}
}

void PaintHelper::paint(QPainter *painter, QPaintEvent *event, bool mouse_moved,
		int mousex, int mousey, double scale_im)
 {

	UNUSED(event);

	const int win_offset = gl_win_size/2;

	double a = (mousex - win_offset)/10.0;
	double b = (mousey - win_offset)/10.0;

	cout << "colmat_valid = " << colmat_valid << endl;


	double divider = scale_im;

	int xmax = gl_win_size/divider;
	int ymax = gl_win_size/divider;

	cout << "xmax = " << xmax << " ymax = " << ymax << endl;

	double local_scale_x = gl_win_size/((double)xmax);
	double local_scale_y = gl_win_size/((double)ymax);
	int win_offset_x = xmax/2;
	int win_offset_y = ymax/2;

	is_silhouette_point.set_xy_max(xmax, ymax);

	// do
	// x1 = (xrast - win_offset_x)/SCALE * local_scale;
	// and analogously to y1 =

	xrast_to_x.set(1/SCALE * local_scale_x, -win_offset_x/SCALE * local_scale_x);
	yrast_to_y.set(1/SCALE * local_scale_y, -win_offset_y/SCALE * local_scale_y);

	// xrast = win_offset + x * SCALE/local_scale;
	// and analogously to yrast =

	x_to_xrast.set(SCALE/local_scale_x, win_offset_x);
	y_to_yrast.set(SCALE/local_scale_y, win_offset_y);

	cout << "xrast_to_x(x_to_xrast(1)) = " << xrast_to_x(x_to_xrast(1)) << endl;

	//cin.get();


	assert(xrast_to_x(x_to_xrast(1)) == 1);
	assert(yrast_to_y(y_to_yrast(1)) == 1);
	assert(xrast_to_x(x_to_xrast(0)) == 0);
	assert(yrast_to_y(y_to_yrast(0)) == 0);

	cout << "divider = " << divider << endl;

	if (!colmat_valid) {

		paint_compute_colmats(a, b, xmax, ymax);
	}

    //if (!mouse_moved)
    //	painter->fillRect(event->rect(), background);

	// painter->translate(win_size/2, win_size/2);

	painter->scale(scale_im, scale_im);


    //painter->save();
    //painter->setBrush(circleBrush);
    //painter->setPen(textPen);

    if (mouse_moved) {
		char buf[128];
		sprintf(buf, "%.3f %.3f", a, b);
		painter->setPen(QPen(Qt::blue));
		painter->fillRect(99,100 - 12, 150, 13, background);
		painter->drawText(100,100, QString(buf));
    }

    int x;
    int y;

    bool is_pen_image = displayFlowLines && (shade_type == 0);

    if (! is_pen_image) {

		paint_from_colmat(painter, x, y, xmax, ymax, divider);
    }

    if (colmat_valid && displayFlowLines) {

    	// draw silhouette points

    	if (shade_type == 0) {

			paint_silhouette_line(painter, xmax, ymax);
    	}

    	// draw cross-field

    	if (displayCrossField) {

			paint_display_crossfield(xmax, ymax, painter);
    	}
    	// get mouse (x,y) as start-point for integrate flow-line

    	if (full_plot1 && full_plot2) {

			paint_fullplot12(painter);
    	}


	}

    char buf[128];
    sprintf(buf, "%.3f %.3f", a, b);
    painter->setPen(QPen(Qt::blue));
    painter->drawText(100,100, QString(buf));

}
