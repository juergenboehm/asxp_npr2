
#ifndef __mainwindow_h
#define __mainwindow_h

#include <QMainWindow>

#include "configfile.h"
#include "asxp.h"

class QLabel;
class QWidget;

class QCheckBox;
class QComboBox;

class QAction;
class QActionGroup;
class QLabel;
class QMenu;
class QGroupBox;

class QLineEdit;
class QSpinBox;
class QDoubleSpinBox;

class QTabWidget;

class SlidersGroup;

class ScreenWidget;
class ScreenGLWidget;

class QProgressBar;


class Window : public QMainWindow
 {
    Q_OBJECT

public:
    Window();
    virtual ~Window();

protected:
    void contextMenuEvent(QContextMenuEvent *event);


private slots:

    void newFile();
    void open();
    void save();
    void print();
    void undo();
    void redo();
    void cut();
    void copy();
    void paste();
    void bold();
    void italic();

    void setInputData(QString fname);

    void taubinStepProc(int sel_ind);
    void manifoldProc(int sel_ind);

    void displaySilhouette();
    void displayColors();
    void displayDots();

    void displayFlowLines();
    void displayFlowLines2(int checked);

    void resetFlowfield();

    void calcTriangulation();
    void exportOpenscad();
    void reduceTriangulation();

    void selectSurf(int sel_ind);
    void flowstreamMeth(int sel_ind);
    void displayCrossfield(int checked);

    void doLineInput();

    void center();
    void setLineSpacing();
    void setParagraphSpacing();
    void about();
    void aboutQt();

    void processTab(int val);



private:
    void createMenus();
    void createActions();

    PaintHelper helper;

    SlidersGroup* euler_sliders;
    QGroupBox* euler_sliders_box;
    QTabWidget* screen_tab;
    QTabWidget* right_tab;

    ScreenWidget *screen_main;
    ScreenGLWidget *screen_main_gl;


    QMenu *fileMenu;
    QMenu *editMenu;
    QMenu *viewMenu;
    QMenu *helpMenu;


    QActionGroup *alignmentGroup;
    QAction *newAct;
    QAction *openAct;
    QAction *saveAct;
    QAction *printAct;
    QAction *exitAct;
    QAction *undoAct;
    QAction *redoAct;
    QAction *cutAct;
    QAction *copyAct;
    QAction *pasteAct;
    QAction *boldAct;
    QAction *italicAct;

    QAction *displaySilhouetteAct;
    QAction *displayColorsAct;
    QAction *displayDotsAct;

    QAction *triangularizeAct;
    QAction *openscadAct;
    QAction *reducetriAct;

    QAction *displayFlowLinesAct;


    QAction *centerAct;
    QAction *setLineSpacingAct;
    QAction *setParagraphSpacingAct;
    QAction *aboutAct;
    QAction *aboutQtAct;

    QProgressBar *statusProgressBar;

    QLabel *infoLabel;

    QLabel *displayLabel;

    QLineEdit * line_input;

    QComboBox *select_surf_box;

    QComboBox *taubin_step_box;

    QComboBox *manifold_sel_box;

    QSpinBox *red_tri_perc_box;
    QDoubleSpinBox *thickness_box;
    QDoubleSpinBox *normal_thresh_box;


    QCheckBox *display_flowlines_cb;

    InputData config_obj;


};



#endif
