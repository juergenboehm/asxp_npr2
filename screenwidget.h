
#ifndef __screenwidget_h
#define __screenwidget_h

//#include <QGLWidget>
#include <QtWidgets>

class PaintHelper;

class QLabel;

class QPaintEvent;
class QWidget;

class ScreenWidget : public QWidget
 {
    Q_OBJECT

public:
    ScreenWidget(PaintHelper *helper, QWidget *parent);

    double get_redscale_factor();


    void setMousexy(QMouseEvent* event);
    void setEulerPaintHelper();

    void rotate();
    void refresh();


public slots:
    void animate();

protected:
    void paintEvent(QPaintEvent *event);

    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);


public:
    PaintHelper *paint_helper;
    QPoint pmouse_akt;

    int elapsed;

    bool mouse_moved;
    bool colmat_valid;

    double euler_phi;
    double euler_theta;
    double euler_psi;

    double zsobel_thresh;
    double nsobel_thresh;

    bool displayFlowLines;

    int shade_type;

    int mousex, mousey;

    bool is_pressed_non_left;

    QLabel *main_info_label;


    double scale;


 };

const int gl_win_size = 800;


#endif
