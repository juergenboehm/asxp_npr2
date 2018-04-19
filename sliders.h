
#ifndef __sliders_h
#define __sliders_h

#include <QGroupBox>
#include <QGroupBox>

 class QDial;
 class QScrollBar;
 class QSlider;
 class ScreenWidget;
 class ScreenGLWidget;

 class SlidersGroup : public QGroupBox
 {
     Q_OBJECT

 public:
     SlidersGroup(Qt::Orientation orientation, ScreenWidget * screen_maina, ScreenGLWidget * screen_gla,
                  QWidget *parent = 0);

 signals:
     void valueChanged(int value);

 public slots:
     void setValue(int value);
     void releaseSlider();


     void setMinimum(int value);
     void setMaximum(int value);
     void invertAppearance(bool invert);
     void invertKeyBindings(bool invert);

public:

     void setSliders(double euler_phi, double euler_theta, double euler_psi);

     int tabval;

     bool block_setValue;

     QSlider *slider_phi;
     QSlider *slider_theta;
     QSlider *slider_psi;

     QSlider *slider_zsobel;
     QSlider *slider_nsobel;

     ScreenWidget* screen_main;
     ScreenGLWidget* screen_gl;
 };


#endif
