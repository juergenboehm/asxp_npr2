


#include <QtWidgets>
#include <QtWidgets/QGroupBox>


#include "screenwidget.h"
#include "glwidget.h"
#include "sliders.h"

 SlidersGroup::SlidersGroup(Qt::Orientation orientation, ScreenWidget * screen_maina, ScreenGLWidget* screen_gla,
                            QWidget *parent)
     : QGroupBox(parent), screen_main(screen_maina), screen_gl(screen_gla)
 {
     slider_phi = new QSlider(orientation);
     slider_phi->setFocusPolicy(Qt::StrongFocus);
     slider_phi->setTickPosition(QSlider::TicksBothSides);
     slider_phi->setTickInterval(10);
     slider_phi->setSingleStep(1);

     slider_theta = new QSlider(orientation);
     slider_theta->setFocusPolicy(Qt::StrongFocus);
     slider_theta->setTickPosition(QSlider::TicksBothSides);
     slider_theta->setTickInterval(10);
     slider_theta->setSingleStep(1);

     slider_psi = new QSlider(orientation);
     slider_psi->setFocusPolicy(Qt::StrongFocus);
     slider_psi->setTickPosition(QSlider::TicksBothSides);
     slider_psi->setTickInterval(10);
     slider_psi->setSingleStep(1);

     slider_nsobel = new QSlider(orientation);
     slider_nsobel->setFocusPolicy(Qt::StrongFocus);
     slider_nsobel->setTickPosition(QSlider::TicksBothSides);
     slider_nsobel->setTickInterval(10);
     slider_nsobel->setSingleStep(1);

     slider_zsobel = new QSlider(orientation);
     slider_zsobel->setFocusPolicy(Qt::StrongFocus);
     slider_zsobel->setTickPosition(QSlider::TicksBothSides);
     slider_zsobel->setTickInterval(10);
     slider_zsobel->setSingleStep(1);


     slider_phi->setMinimum(0);
     slider_theta->setMinimum(0);
     slider_psi->setMinimum(0);
     slider_nsobel->setMinimum(0);
     slider_zsobel->setMinimum(0);

     slider_phi->setMaximum(360);
     slider_theta->setMaximum(360);
     slider_psi->setMaximum(360);
     slider_nsobel->setMaximum(1000);
     slider_zsobel->setMaximum(1000);

     connect(slider_phi, SIGNAL(valueChanged(int)), slider_phi, SLOT(setValue(int)));
     connect(slider_theta, SIGNAL(valueChanged(int)), slider_theta, SLOT(setValue(int)));
     connect(slider_psi, SIGNAL(valueChanged(int)), slider_psi, SLOT(setValue(int)));
     connect(slider_zsobel, SIGNAL(valueChanged(int)), slider_zsobel, SLOT(setValue(int)));
     connect(slider_nsobel, SIGNAL(valueChanged(int)), slider_nsobel, SLOT(setValue(int)));

     connect(slider_phi, SIGNAL(valueChanged(int)), this, SLOT(setValue(int)));
     connect(slider_theta, SIGNAL(valueChanged(int)), this, SLOT(setValue(int)));
     connect(slider_psi, SIGNAL(valueChanged(int)), this, SLOT(setValue(int)));
     connect(slider_zsobel, SIGNAL(valueChanged(int)), this, SLOT(setValue(int)));
     connect(slider_nsobel, SIGNAL(valueChanged(int)), this, SLOT(setValue(int)));

     connect(slider_phi, SIGNAL(sliderReleased()), this, SLOT(releaseSlider()));
     connect(slider_theta, SIGNAL(sliderReleased()), this, SLOT(releaseSlider()));
     connect(slider_psi, SIGNAL(sliderReleased()), this, SLOT(releaseSlider()));
     connect(slider_zsobel, SIGNAL(sliderReleased()), this, SLOT(releaseSlider()));
     connect(slider_nsobel, SIGNAL(sliderReleased()), this, SLOT(releaseSlider()));

     QBoxLayout::Direction direction;

     if (orientation == Qt::Horizontal)
         direction = QBoxLayout::TopToBottom;
     else
         direction = QBoxLayout::LeftToRight;

     QBoxLayout *slidersLayout = new QBoxLayout(direction);
     slidersLayout->addWidget(slider_phi);
     slidersLayout->addWidget(slider_theta);
     slidersLayout->addWidget(slider_psi);
     slidersLayout->addWidget(slider_zsobel);
     slidersLayout->addWidget(slider_nsobel);
     setLayout(slidersLayout);

     block_setValue = false;

     setMinimumSize(200, 800);
     setMaximumSize(200, 800);

 }

void SlidersGroup::setSliders(double euler_phi, double euler_theta, double euler_psi)
{

#define TO_DEG(x) ((x)/3.14159265 * 180.0)

	block_setValue = true;

	slider_phi->setValue(TO_DEG(euler_phi));
	slider_theta->setValue(TO_DEG(euler_theta));
	slider_psi->setValue(TO_DEG(euler_psi));

	repaint();

	block_setValue = false;

}


 void SlidersGroup::releaseSlider()
 {
	 if (tabval == 0) {
		 screen_main->refresh();
	 }
 }

 void SlidersGroup::setValue(int value)
 {
	 const double Pi = 3.14159265;

	 double euler_phi;
	 double euler_theta;
	 double euler_psi;

	 double zsobel;
	 double nsobel;

	 if (block_setValue) {
		 return;
	 }

     int phi_val = slider_phi->value();
     int theta_val = slider_theta->value();
     int psi_val = slider_psi->value();

     int zsobel_val = slider_zsobel->value();
     int nsobel_val = slider_nsobel->value();

#define TO_RAD(x) ((double (x))/360 * 2 * Pi)

     euler_phi = TO_RAD(phi_val);
     euler_theta = TO_RAD(theta_val);
     euler_psi = TO_RAD(psi_val);

#define TO_UNIT(x) ((double (x))/1000)

     zsobel = TO_UNIT(zsobel_val);
     nsobel = TO_UNIT(nsobel_val);

     if (tabval == 0) {
		 screen_main->euler_phi = euler_phi;
		 screen_main->euler_theta = euler_theta;
		 screen_main->euler_psi = euler_psi;

	     screen_main->zsobel_thresh = zsobel;
	     screen_main->nsobel_thresh = nsobel;
	     screen_main->rotate();
     }

     if (tabval == 1) {
		 screen_gl->euler_phi = euler_phi;
		 screen_gl->euler_theta = euler_theta;
		 screen_gl->euler_psi = euler_psi;
		 screen_gl->rotate();
     }
 }

 void SlidersGroup::setMinimum(int value)
 {
     slider_phi->setMinimum(value);
     slider_theta->setMinimum(value);
     slider_psi->setMinimum(value);

     slider_zsobel->setMinimum(value);
     slider_nsobel->setMinimum(value);

 }

 void SlidersGroup::setMaximum(int value)
 {
     slider_phi->setMaximum(value);
     slider_theta->setMaximum(value);
     slider_psi->setMaximum(value);

     slider_zsobel->setMaximum(value);
     slider_nsobel->setMaximum(value);

 }

 void SlidersGroup::invertAppearance(bool invert)
 {
     slider_phi->setInvertedAppearance(invert);
     slider_theta->setInvertedAppearance(invert);
     slider_psi->setInvertedAppearance(invert);

     slider_zsobel->setInvertedAppearance(invert);
     slider_nsobel->setInvertedAppearance(invert);

 }

 void SlidersGroup::invertKeyBindings(bool invert)
 {
     slider_phi->setInvertedControls(invert);
     slider_theta->setInvertedControls(invert);
     slider_psi->setInvertedControls(invert);

     slider_zsobel->setInvertedControls(invert);
     slider_nsobel->setInvertedControls(invert);
 }
