
#ifndef __glwidget_h
#define __glwidget_h

#include <QtOpenGL/QGLWidget>


class ScreenGLWidget : public QGLWidget
{
	Q_OBJECT

public:
    ScreenGLWidget(QWidget *parent);

    void rotate();

    double euler_phi;
    double euler_theta;
    double euler_psi;

protected:

    void initializeGL();

    void resizeGL(int w, int h);

    void paintGL();

};



















#endif
