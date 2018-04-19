
#include <QTimer>
#include <QMouseEvent>
#include <Qt>
#include <QLabel>


#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "screenwidget.h"
#include "asxp.h"

using namespace std;

const double redscale_factor = 8.0;

ScreenWidget::ScreenWidget(PaintHelper *helper, QWidget *parent)
//    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent), paint_helper(helper)
	: QWidget( parent), paint_helper(helper)
 {
    elapsed = 0;
    mouse_moved = false;
    colmat_valid = true;
    setFixedSize(gl_win_size, gl_win_size);
    //setAutoFillBackground(false);

    //setBackgroundRole(QPalette::Base);
    setAutoFillBackground(true);

    setMouseTracking(true);

    scale = 1.0;

    mousex = 0;
    mousey = 0;

    euler_phi = 0;
    euler_theta = 0;
    euler_psi = 0;

    zsobel_thresh = 0;
    nsobel_thresh = 0;

    displayFlowLines = false;

    shade_type = 0; // silhouette

    is_pressed_non_left = false;

    setEulerPaintHelper();
}

double ScreenWidget::get_redscale_factor()
{
	if (paint_helper->shade_type == 2) {
		return 8.0;
	} else if (paint_helper->shade_type == 1) {
		return 1.0;
	} else {
		return 8.0;
	}
}

void ScreenWidget::setEulerPaintHelper()
{
	paint_helper->euler_phi = euler_phi;
	paint_helper->euler_theta = euler_theta;
	paint_helper->euler_psi = euler_psi;

	paint_helper->zsobel_thresh = zsobel_thresh;
	paint_helper->nsobel_thresh = nsobel_thresh;

	paint_helper->shade_type = shade_type;

	paint_helper->displayFlowLines = displayFlowLines;

	std::cout << "zsobel_thresh = " << zsobel_thresh << " nsobel_thresh = " << nsobel_thresh << std::endl;

}

void ScreenWidget::setMousexy(QMouseEvent* event)
{
	mousex = event->x();
	mousey = event->y();
}

void ScreenWidget::animate()
 {
    elapsed = (elapsed + qobject_cast<QTimer*>(sender())->interval()) % 10000;
    repaint();
}

void ScreenWidget::paintEvent(QPaintEvent *event)
 {
    QPainter painter;
    painter.begin(this);

    //painter.setRenderHint(QPainter::Antialiasing);

    paint_helper->colmat_valid = colmat_valid;
    paint_helper->paint(&painter, event, mouse_moved, mousex, mousey, scale);
    painter.end();
}

void ScreenWidget::mousePressEvent(QMouseEvent *event)
{
	if (!(event->buttons() & Qt::LeftButton)) {
		is_pressed_non_left = true;
		return;
	}

	is_pressed_non_left = false;

	setMousexy(event);
	setEulerPaintHelper();

	printf("mouse press x = %d, y= %d\n", mousex, mousey);
	mouse_moved = false;
	repaint();
}

void ScreenWidget::mouseMoveEvent(QMouseEvent *event)
{
	if (!(event->buttons() & Qt::LeftButton)) {

		int mousex = event->x();
		int mousey = event->y();

		//cout << "mouse x raw = " << mousex << " mouse y raw = " << mousey << endl;

		mousex = ::max(0, ::min(mousex, 799));
		mousey = ::max(0, ::min(mousey, 799));

		double l1, l2, kg, kh;

		paint_helper->getImageInfo(mousex, mousey, l1, l2, kg, kh);

		char buf[256];

		sprintf(buf, "nsel = %f jsel = %f", l1, l2);

		//cout << buf << endl;

		main_info_label->setText(buf);
		main_info_label->repaint();

		return;

	}

	if (displayFlowLines){
		return;
	}

	is_pressed_non_left = false;

	setMousexy(event);
	setEulerPaintHelper();

	mouse_moved = true;
	printf("mouse move x = %d, y= %d\n", mousex, mousey);

	colmat_valid = false;

	scale = get_redscale_factor();
	repaint();
	scale = 1.0;
	colmat_valid = true;
	mouse_moved = false;
}

void ScreenWidget::rotate()
{
	setEulerPaintHelper();

	mouse_moved = false;

	colmat_valid = false;
	scale = get_redscale_factor();
	repaint();
	scale = 1.0;
	colmat_valid = true;

}

void ScreenWidget::refresh()
{
	setEulerPaintHelper();

	mouse_moved = false;

	colmat_valid = false;

	scale = 1.0;

	std::cout << "refresh: DisplayFlowlines = " << displayFlowLines << std::endl;

	repaint();
	colmat_valid = true;

}


void ScreenWidget::mouseReleaseEvent(QMouseEvent *event)
{
	if (is_pressed_non_left) {
		return;
	}

	if (displayFlowLines) {
		return;
	}

	setMousexy(event);
	setEulerPaintHelper();

	printf("mouse release x = %d, y= %d\n", mousex, mousey);

	mouse_moved = false;
	colmat_valid = false;
	repaint();
	colmat_valid = true;
}

