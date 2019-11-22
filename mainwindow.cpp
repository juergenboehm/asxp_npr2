
#include <iostream>

#include <stdlib.h>
#include <stdio.h>


#include <QToolButton>
#include <QLineEdit>
#include <QString>

#include <QtWidgets>


#include <QPrinter>
#include <QPrintDialog>

#include <QtSvg/QSvgGenerator>


#include <QTimer>

#include <QToolBar>


#include <QTabWidget>

#include "configfile.h"

#include "streamplot.h"

#include "sliders.h"
#include "screenwidget.h"
#include "glwidget.h"
#include "openscad.h"
#include "asxp.h"
#include "procpoly.h"

#include "triangularize.h"
#include "gtssurface.h"


#include "mainwindow.h"

using namespace std;

QTabWidget* screen_tab_global;



Window::Window()
    : QMainWindow()
 {

	QWidget *widget = new QWidget();
	setCentralWidget(widget);


	screen_main = new ScreenWidget(&helper, this);
	screen_main_gl = new ScreenGLWidget(this);

    euler_sliders = new SlidersGroup(Qt::Vertical, screen_main, screen_main_gl);
    euler_sliders_box = new QGroupBox(tr("Eulerwinkel"));

    QGridLayout *euler_sliders_layout = new QGridLayout();

    euler_sliders_layout->addWidget(euler_sliders, 0, 0);

    euler_sliders_box->setLayout(euler_sliders_layout);
    euler_sliders_box->setFixedSize(200, 845);


    screen_tab = new QTabWidget();
    screen_tab_global = screen_tab;

    screen_tab->setFixedSize(810, 850);

    screen_tab->addTab(screen_main, tr("Main"));
    screen_tab->addTab(screen_main_gl, tr("AltGL"));

    connect(screen_tab, SIGNAL(currentChanged(int)), this, SLOT(processTab(int)));

    euler_sliders->tabval = 0;


    // right controls is in its own tab on the right

    QVBoxLayout *right_controls = new QVBoxLayout();

    // taubin label

	QLabel* taubin_label = new QLabel(tr("Taubin steps"));

	// taubin step box

    taubin_step_box = new QComboBox();

    int taubin_vals[]  = { 0, 10, 20, 30, 50, 70, 100 };
    int taubin_len = sizeof(taubin_vals)/sizeof(int);

    for(int i = 0; i < taubin_len; ++ i) {

    	taubin_step_box->addItem(QString::number(taubin_vals[i]));
    }

	connect(taubin_step_box, SIGNAL(currentIndexChanged(int)), this, SLOT(taubinStepProc(int)));

	// manifold select box

	manifold_sel_box = new QComboBox();

	manifold_sel_box->addItem(tr("Manifold"));
	manifold_sel_box->addItem(tr("Manifold with boundary"));
	manifold_sel_box->addItem(tr("Non Manifold"));

	connect(manifold_sel_box, SIGNAL(currentIndexChanged(int)), this, SLOT(manifoldProc(int)));


	// flowfield reset button


    QToolButton *ffr_button = new QToolButton(this);
    ffr_button->setText(tr("Reset Flowstreams"));

    connect(ffr_button, SIGNAL(clicked()), this, SLOT(resetFlowfield()));


    // flowstream method combo box


    QComboBox *flowstream_meth_box = new QComboBox();

    flowstream_meth_box->addItem(tr("Mebarki"));
    flowstream_meth_box->addItem(tr("Jobard"));


	connect(flowstream_meth_box, SIGNAL(currentIndexChanged(int)), this, SLOT(flowstreamMeth(int)));

	// display crossfield checkbox

	QCheckBox *display_crossfield_cb = new QCheckBox(tr("display Crossfield"));
	display_crossfield_cb->setCheckState(Qt::Checked);
	connect(display_crossfield_cb, SIGNAL(stateChanged(int)), this, SLOT(displayCrossfield(int)));


	// display flowlines checkbox

	display_flowlines_cb = new QCheckBox(tr("display Flowlines"));
	display_flowlines_cb->setCheckState(Qt::Unchecked);
	connect(display_flowlines_cb, SIGNAL(stateChanged(int)), this, SLOT(displayFlowLines2(int)));

	display_flowlines_cb->setEnabled(false);


    // display label

    displayLabel = new QLabel(tr("<i>Choose a menu option, or right-click to "
                              "invoke a context menu</i>"));
    displayLabel->setFrameStyle(QFrame::StyledPanel | QFrame::Sunken);
    displayLabel->setAlignment(Qt::AlignCenter);

    helper.displayLabel = displayLabel;

    // line input

    line_input = new QLineEdit(this);

    connect(line_input, SIGNAL(returnPressed()), this, SLOT(doLineInput()));

    // red_tri_perc_box label

	QLabel* red_tri_perc_box_label = new QLabel(tr("Tri-reduction percentage"));


    // reduceTriangle % input

    red_tri_perc_box = new QSpinBox();
    red_tri_perc_box->setMinimum(1);
    red_tri_perc_box->setMaximum(100);
    red_tri_perc_box->setSingleStep(1);
    // Will increment the current value with 1 (if you use up arrow key) (if you use down arrow key => -1)
    red_tri_perc_box->setValue(40);// Default/begining value

    // Thickness of shell label

	QLabel* thickness_label = new QLabel(tr("Shell thickness"));


    // Thickness of shell input

    thickness_box = new QDoubleSpinBox();
    thickness_box->setMinimum(0.1);
    thickness_box->setMaximum(2.0);
    thickness_box->setSingleStep(0.05);
    // Will increment the current value with 0.05 (if you use up arrow key) (if you use down arrow key => -0.05)
    thickness_box->setValue(0.2);// Default/begining value


    // Normal threshold for critical points on surface

 	QLabel* normal_thresh_label = new QLabel(tr("Normal threshold"));


     // Thickness of shell input

     normal_thresh_box = new QDoubleSpinBox();
     normal_thresh_box->setMinimum(0.0);
     normal_thresh_box->setMaximum(50.0);
     normal_thresh_box->setSingleStep(0.1);
     // Will increment the current value with 0.1 (if you use up arrow key) (if you use down arrow key => -0.1)
     normal_thresh_box->setValue(5.0);// Default/begining value



    // right controls

	right_controls->addWidget(taubin_label);
    right_controls->addWidget(taubin_step_box);
    right_controls->addWidget(ffr_button);
    right_controls->addWidget(flowstream_meth_box);
    right_controls->addWidget(display_crossfield_cb);
    right_controls->addWidget(display_flowlines_cb);

    right_controls->addWidget(displayLabel);
    right_controls->addWidget(line_input);

    right_controls->addWidget(red_tri_perc_box_label);
    right_controls->addWidget(red_tri_perc_box);

    right_controls->addWidget(thickness_label);
    right_controls->addWidget(thickness_box);

    right_controls->addWidget(normal_thresh_label);
    right_controls->addWidget(normal_thresh_box);

    right_controls->addStretch(10);

    QGroupBox* right_controls_box = new QGroupBox(tr("Right Controls"));
    right_controls_box->setLayout(right_controls);

    right_controls_box->setFixedSize(200, 500);


    QHBoxLayout *main_layout = new QHBoxLayout();

    right_tab = new QTabWidget();

    right_tab->addTab(euler_sliders_box, tr("Euler"));
    right_tab->addTab(right_controls_box, tr("Control"));

    right_tab->setFixedSize(300, 850);


    main_layout->addWidget(screen_tab);
    main_layout->addWidget(right_tab);

    widget->setLayout(main_layout);

    // select surface box

    select_surf_box = new QComboBox();

	connect(select_surf_box, SIGNAL(currentIndexChanged(int)), this, SLOT(selectSurf(int)));


	// infoLabel is not displayed.
	// could be included in ToolBar, for example
	// but not done, because it does not provide interesting information.

    infoLabel = new QLabel(tr("<i>Choose a menu option, or right-click to "
                              "invoke a context menu</i>"));
    infoLabel->setFrameStyle(QFrame::StyledPanel | QFrame::Sunken);
    infoLabel->setAlignment(Qt::AlignCenter);

    screen_main->main_info_label = displayLabel;


    createActions();
    createMenus();

	QToolBar* toolbar = addToolBar("Main");

    toolbar->addAction(openAct);
    toolbar->addAction(printAct);
    toolbar->addSeparator();
    toolbar->addAction(triangularizeAct);
    toolbar->addAction(openscadAct);
    toolbar->addAction(reducetriAct);
    toolbar->addWidget(manifold_sel_box);
    toolbar->addSeparator();
    toolbar->addWidget(select_surf_box);
    toolbar->addSeparator();


    statusProgressBar = new QProgressBar(this);
    statusProgressBar->setFixedWidth(200);
    statusProgressBar->setFixedHeight(12);

    mainStatusProgressBar = statusProgressBar;

    this->statusBar()->addPermanentWidget(statusProgressBar);

    QString message = tr("A context menu is available by right-clicking");
    statusBar()->showMessage(message);

    setWindowTitle(tr("asxp"));

    setMinimumSize(1150, 950);
    resize(1150, 950);

	setInputData("asxp.conf");

	init_complexes();
	taubin_steps = 0;

	manifold_sel = 1;
	manifold_sel_box->setCurrentIndex(manifold_sel);

 }

Window::~Window()
{
	remove_complexes();
}

void Window::processTab(int val)
{
	double euler_phi, euler_theta, euler_psi;

	cout << "processTab val = " << val << endl;
	euler_sliders->tabval = val;

	if (val == 0) {
		euler_phi = screen_main->euler_phi;
		euler_theta = screen_main->euler_theta;
		euler_psi = screen_main->euler_psi;

		select_surf_box->setDisabled(false);
	}

	if (val == 1) {
		euler_phi = screen_main_gl->euler_phi;
		euler_theta = screen_main_gl->euler_theta;
		euler_psi = screen_main_gl->euler_psi;

		select_surf_box->setDisabled(true);
	}

	euler_sliders->setSliders(euler_phi, euler_theta, euler_psi);

}

void Window::contextMenuEvent(QContextMenuEvent *event)
 {
    QMenu menu(this);

    menu.addAction(displaySilhouetteAct);
    menu.addAction(displayColorsAct);
    menu.addAction(displayDotsAct);

    menu.addAction(displayFlowLinesAct);

    menu.exec(event->globalPos());
}


void Window::newFile()
 {
    infoLabel->setText(tr("Invoked <b>File|New</b>"));
    cout << "newFile selected." << endl;
 }

void Window::setInputData(QString fname)
{
	config_obj.set_fname(fname);

	int cnt_blocks = 0;
	config_obj.readSettings(cnt_blocks);

	cout << "cnt_blocks = " << cnt_blocks << endl;

	QList<QString> surf_list;
	config_obj.getSurfBlocks(surf_list);
	QString overview = "";

	select_surf_box->clear();

	for(QList<QString>::iterator it = surf_list.begin(); it != surf_list.end(); ++it){
		select_surf_box->addItem(*it);
	}
}

void Window::open()
 {
    infoLabel->setText(tr("Invoked <b>File|Open</b>"));

    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    ".",
                                                    tr("Configs (*.conf)"));

    if (fileName != "") {
        config_obj.reset();
    	setInputData(fileName);
    }
}

void Window::save()
 {
    infoLabel->setText(tr("Invoked <b>File|Save</b>"));
}

void Window::print()
{
    infoLabel->setText(tr("Invoked <b>File|Print</b>"));

	QPrinter printer;

	QPrintDialog *dialog = new QPrintDialog(&printer, this);

	dialog->setWindowTitle(tr("Print Document"));

	if (dialog->exec() != QDialog::Accepted)
		return;

	printer.setResolution(300);

	//printer.setOutputFormat(QPrinter::PostScriptFormat);

	int printer_width = printer.pageRect().width();
	int printer_height = printer.pageRect().height();

	printf("printer width = %d, printer height = %d\n", printer_width, printer_height );

	//return;

	QPainter painter;
	painter.begin(&printer);



	double xscale = printer_width / double(screen_main->width());
	double yscale = printer_height / double(screen_main->height());

	xscale = 2;
	yscale = 2;

	double scale = qMin(xscale, yscale);

	painter.translate(printer.paperRect().x() + printer.pageRect().width() / 2,
			printer.paperRect().y() + printer.pageRect().height() / 2);

	painter.scale(scale, scale);

	painter.translate(-screen_main->width() / 2, -screen_main->height() / 2);

	screen_main->scale = 0.5;
	screen_main->colmat_valid = false;

	double old_dsep = pp.dsep;

	pp.dsep = 2 * old_dsep;



	screen_main->paint_helper->displayFlowLines = true;
	screen_main->paint_helper->displayCrossField = false;
	screen_main->paint_helper->streamLineColor = 0;

	screen_main->paint_helper->shade_type = -1;

#if 1
	screen_main->render(&painter);
#endif

	screen_main->colmat_valid = true;

	painter.end();

#ifdef WITH_SVG

	QSvgGenerator generator;
	generator.setFileName("print.svg");
	generator.setSize(QSize(1600, 1600));
	generator.setViewBox(QRect(0, 0, 1600, 1600));
	generator.setTitle(tr("SVG output of asxp"));
	generator.setDescription(tr("An SVG drawing created by the SVG Generator "));

	QPainter painter_svg;
	painter_svg.begin(&generator);

	screen_main->render(&painter_svg);

	painter_svg.end();

#endif


	pp.dsep = old_dsep;

	screen_main->scale = 1.0;
	screen_main->colmat_valid = false;



}

void Window::undo()
 {
    infoLabel->setText(tr("Invoked <b>Edit|Undo</b>"));
}

void Window::redo()
 {
    infoLabel->setText(tr("Invoked <b>Edit|Redo</b>"));
}

void Window::cut()
 {
    infoLabel->setText(tr("Invoked <b>Edit|Cut</b>"));
}

void Window::copy()
 {
    infoLabel->setText(tr("Invoked <b>Edit|Copy</b>"));
}

void Window::paste()
 {
    infoLabel->setText(tr("Invoked <b>Edit|Paste</b>"));
}

void Window::bold()
 {
    infoLabel->setText(tr("Invoked <b>Edit|Format|Bold</b>"));
}

void Window::italic()
 {
    infoLabel->setText(tr("Invoked <b>Edit|Format|Italic</b>"));
}

void Window::displaySilhouette()
{
	screen_main->shade_type = 0;

	displayFlowLinesAct->setEnabled(true);
	display_flowlines_cb->setEnabled(true);


	screen_main->refresh();

    infoLabel->setText(tr("Invoked <b>Edit|Format|display Silhouette</b>"));
}

void Window::displayColors()
 {
	screen_main->shade_type = 1;

	screen_main->displayFlowLines = false;

	displayFlowLinesAct->setChecked(false);

	displayFlowLinesAct->setEnabled(false);

	display_flowlines_cb->setChecked(false);
	display_flowlines_cb->setCheckState(Qt::Unchecked);
	display_flowlines_cb->setEnabled(false);

	screen_main->refresh();

    infoLabel->setText(tr("Invoked <b>Edit|Format|display Colors</b>"));
}

void Window::displayDots()
 {
	screen_main->shade_type = 2;

	screen_main->displayFlowLines = false;

	displayFlowLinesAct->setChecked(false);

	displayFlowLinesAct->setEnabled(false);

	display_flowlines_cb->setChecked(false);
	display_flowlines_cb->setCheckState(Qt::Unchecked);
	display_flowlines_cb->setEnabled(false);

	screen_main->refresh();

    infoLabel->setText(tr("Invoked <b>Edit|Format|display Dots</b>"));
}

void Window::displayFlowLines()
{

	cout << "displayFlowLines called...." << endl;

	bool new_disp = !screen_main->displayFlowLines;

	display_flowlines_cb->toggle();

//	screen_main->displayFlowLines = new_disp;
//	screen_main->refresh();

	infoLabel->setText(tr("Invoked <b>Edit|Format|display Flowlines</b>"));

}

void Window::displayFlowLines2(int checked)
{

	cout << "displayFlowLines2 called...." << endl;

	screen_main->displayFlowLines = (checked == Qt::Checked);

	displayFlowLinesAct->setChecked(checked == Qt::Checked);
	screen_main->refresh();
}



void Window::resetFlowfield()
{
	Streamplot* full_plot1 = pp.full_plot1;
	Streamplot* full_plot2 = pp.full_plot2;

	bool do_repaint = true;

	if (full_plot1 != 0) {

		delete pp.full_plot1;
		pp.full_plot1 = 0;

		do_repaint = true;

		screen_main->mousex = 0;
		screen_main->mousey = 0;

	}
	if (full_plot2 != 0) {

		delete pp.full_plot2;
		pp.full_plot2 = 0;

		do_repaint = true;

		screen_main->mousex = 0;
		screen_main->mousey = 0;

	}

	pp.compute_streamfield(gl_win_size, gl_win_size);
	screen_main->colmat_valid = true;

	if (do_repaint)
		screen_main->repaint();

}




void Window::calcTriangulation()
{
	displayLabel->setText(tr("Start triangularize"));

	double normal_thresh = normal_thresh_box->value();

	triangularize(normal_thresh);

	displayLabel->setText(tr("Finished triangularize"));
}

void Window::exportOpenscad()
{
	if (trisurf.vertex_vec.size() > 0) {

		double shell_thickness = thickness_box->value();

		std::cout << "build openscad" << std::endl;

		std::ofstream ofscad("model.scad");

		build_openscad(ofscad, shell_thickness);

		std::cout << "end openscad." << std::endl;

	} else {

		std::cout << "no triangulation: export openSCAD not possible." << std::endl;
	}
}

void Window::reduceTriangulation()
{
	if (trisurf.vertex_vec.size() > 0) {

		double red_perc = double(red_tri_perc_box->value())/100.0;
		double max_dist;

		std::cout << "reduce Triangulation selected: " << red_perc << " ..." << std::endl;

		reduce_gts_surface(trisurf, red_perc, max_dist);

		displayLabel->setText("max_dist = " + QString::number(max_dist));

	} else {

		std::cout << "no triangulation: reduce not possible." << std::endl;

	}
}

void Window::selectSurf(int sel_ind)
{
	int cnt = 0;
	QString f5str;

	cout << "***** selectSurf called: sel_ind = " << sel_ind  << endl;

	for(InputData::SettingsList::iterator it = config_obj.settings_list.begin();
			it != config_obj.settings_list.end(); ++it) {

		if (sel_ind == cnt) {
			f5str = it->smap["poly"];
			int man_sel = 1;
			if (it->smap.count("manifold") > 0) {
				man_sel = it->smap["manifold"].toInt();
			}
			if (0 <= man_sel && man_sel <= 2) {
				manifold_sel_box->setCurrentIndex(man_sel);
			}
			break;
		}
		++cnt;
	}

	pp.init_f5(f5str);

	screen_main->colmat_valid = false;
	screen_main->refresh();

}

void Window::flowstreamMeth(int sel_ind)
{
	if (sel_ind == 0) {
		pp.streamgen_type = 1;

	} else if (sel_ind == 1) {
		pp.streamgen_type = 0;
	}
}

void Window::taubinStepProc(int sel_ind)
{
	taubin_steps = (taubin_step_box->itemText(sel_ind)).toInt();
}

void Window::manifoldProc(int sel_ind)
{
	manifold_sel = sel_ind;
}


void Window::displayCrossfield(int checked)
{
	screen_main->paint_helper->displayCrossField = (checked == Qt::Checked);
	screen_main->repaint();

}

void Window::doLineInput()
{
	QString new_text = line_input->text();

	new_text += "                 ";

	QByteArray ba = new_text.toLatin1();
	const char *buf = ba.data();

	if (strncmp(buf, "tria", 4) == 0) {

		calcTriangulation();

	} else if (strncmp(buf, "cgal", 4) == 0) {

		pp.streamgen_type = 1;

	} else if (strncmp(buf, "own", 3) == 0) {

		pp.streamgen_type = 0;

	} else if (strncmp(buf, "dsep", 4) == 0) {

		double dsep;

		sscanf(buf + 4, "%lf", &dsep);

		pp.dsep = dsep;

		resetFlowfield();
	} else if (strncmp(buf,"color", 5) == 0) {

		int color;

		sscanf(buf + 5, "%d", &color);

		screen_main->paint_helper->streamLineColor = color;

		screen_main->repaint();

	} else if (strncmp(buf, "cfshow", 6) == 0) {

		int showit;

		sscanf(buf + 6, "%d", &showit);

		screen_main->paint_helper->displayCrossField = (showit == 1);

		screen_main->repaint();

	} else if (strncmp(buf, "help", 4) == 0) {

		displayLabel->setText(tr("dsep <float> | color 0/1/2 | cfshow 0/1 | select 0/1/2/3"));

	} else if (strncmp(buf, "select", 6) == 0) {

		QString f5str;

		if (QString(buf + 7).trimmed().length() == 0) {
			QList<QString> surf_list;
			config_obj.getSurfBlocks(surf_list);
			QString overview = "";
			for(QList<QString>::iterator it = surf_list.begin(); it != surf_list.end(); ++it){
				overview += *it;
				overview += " | ";
			}

			displayLabel->setText(overview);

			return;

		} else {
			QString sel_str(buf + 7);
			sel_str = sel_str.trimmed();

			cout << "sel_str = " << sel_str.toLatin1().data() << endl;
			for(InputData::SettingsList::iterator it = config_obj.settings_list.begin();
					it != config_obj.settings_list.end(); ++it) {
				cout << "it->block_name = " << it->block_name.toLatin1().data() << endl;
				if (it->block_name == sel_str) {
					f5str = it->smap["poly"];
					break;
				}
			}

		}

		pp.init_f5(f5str);

		screen_main->colmat_valid = false;
		screen_main->refresh();

	}
}


void Window::center()
 {
    infoLabel->setText(tr("Invoked <b>Edit|Format|Center</b>"));
}

void Window::setLineSpacing()
 {
    infoLabel->setText(tr("Invoked <b>Edit|Format|Set Line Spacing</b>"));
}

void Window::setParagraphSpacing()
 {
    infoLabel->setText(tr("Invoked <b>Edit|Format|Set Paragraph Spacing</b>"));
}

void Window::about()
 {
    infoLabel->setText(tr("Invoked <b>Help|About</b>"));
    QMessageBox::about(this, tr("About Menu"),
            tr("asxp is an <b>Algebraic Surface eXPlorer</b>."));
}

void Window::aboutQt()
 {
    infoLabel->setText(tr("Invoked <b>Help|About Qt</b>"));
}

void Window::createActions()
 {
    newAct = new QAction(tr("&New"), this);
    newAct->setShortcuts(QKeySequence::New);
    newAct->setStatusTip(tr("Create a new file"));
    connect(newAct, SIGNAL(triggered()), this, SLOT(newFile()));

    openAct = new QAction(QIcon("images/open.png"), tr("&Open..."), this);
    openAct->setShortcuts(QKeySequence::Open);
    openAct->setStatusTip(tr("Open an existing file"));
    connect(openAct, SIGNAL(triggered()), this, SLOT(open()));

    saveAct = new QAction(tr("&Save"), this);
    saveAct->setShortcuts(QKeySequence::Save);
    saveAct->setStatusTip(tr("Save the document to disk"));
    connect(saveAct, SIGNAL(triggered()), this, SLOT(save()));

    printAct = new QAction(QIcon("images/print.png"),tr("&Print..."), this);
    printAct->setShortcuts(QKeySequence::Print);
    printAct->setStatusTip(tr("Print the document"));
    connect(printAct, SIGNAL(triggered()), this, SLOT(print()));

    exitAct = new QAction(tr("E&xit"), this);
    exitAct->setShortcuts(QKeySequence::Quit);
    exitAct->setStatusTip(tr("Exit the application"));
    connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));

    undoAct = new QAction(tr("&Undo"), this);
    undoAct->setShortcuts(QKeySequence::Undo);
    undoAct->setStatusTip(tr("Undo the last operation"));
    connect(undoAct, SIGNAL(triggered()), this, SLOT(undo()));

    redoAct = new QAction(tr("&Redo"), this);
    redoAct->setShortcuts(QKeySequence::Redo);
    redoAct->setStatusTip(tr("Redo the last operation"));
    connect(redoAct, SIGNAL(triggered()), this, SLOT(redo()));

    cutAct = new QAction(tr("Cu&t"), this);
    cutAct->setShortcuts(QKeySequence::Cut);
    cutAct->setStatusTip(tr("Cut the current selection's contents to the "
                            "clipboard"));
    connect(cutAct, SIGNAL(triggered()), this, SLOT(cut()));

    copyAct = new QAction(tr("&Copy"), this);
    copyAct->setShortcuts(QKeySequence::Copy);
    copyAct->setStatusTip(tr("Copy the current selection's contents to the "
                             "clipboard"));
    connect(copyAct, SIGNAL(triggered()), this, SLOT(copy()));

    pasteAct = new QAction(tr("&Paste"), this);
    pasteAct->setShortcuts(QKeySequence::Paste);
    pasteAct->setStatusTip(tr("Paste the clipboard's contents into the current "
                              "selection"));
    connect(pasteAct, SIGNAL(triggered()), this, SLOT(paste()));

    boldAct = new QAction(tr("&Bold"), this);
    boldAct->setCheckable(true);
    boldAct->setShortcut(QKeySequence::Bold);
    boldAct->setStatusTip(tr("Make the text bold"));
    connect(boldAct, SIGNAL(triggered()), this, SLOT(bold()));

    QFont boldFont = boldAct->font();
    boldFont.setBold(true);
    boldAct->setFont(boldFont);

    italicAct = new QAction(tr("&Italic"), this);
    italicAct->setCheckable(true);
    italicAct->setShortcut(QKeySequence::Italic);
    italicAct->setStatusTip(tr("Make the text italic"));
    connect(italicAct, SIGNAL(triggered()), this, SLOT(italic()));

    QFont italicFont = italicAct->font();
    italicFont.setItalic(true);
    italicAct->setFont(italicFont);

    setLineSpacingAct = new QAction(tr("Set &Line Spacing..."), this);
    setLineSpacingAct->setStatusTip(tr("Change the gap between the lines of a "
                                       "paragraph"));
    connect(setLineSpacingAct, SIGNAL(triggered()), this, SLOT(setLineSpacing()));

    setParagraphSpacingAct = new QAction(tr("Set &Paragraph Spacing..."), this);
    setLineSpacingAct->setStatusTip(tr("Change the gap between paragraphs"));
    connect(setParagraphSpacingAct, SIGNAL(triggered()),
            this, SLOT(setParagraphSpacing()));

    aboutAct = new QAction(tr("&About"), this);
    aboutAct->setStatusTip(tr("Show the application's About box"));
    connect(aboutAct, SIGNAL(triggered()), this, SLOT(about()));

    aboutQtAct = new QAction(tr("About &Qt"), this);
    aboutQtAct->setStatusTip(tr("Show the Qt library's About box"));
    connect(aboutQtAct, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
    connect(aboutQtAct, SIGNAL(triggered()), this, SLOT(aboutQt()));



    displaySilhouetteAct = new QAction(tr("&Silhouette"), this);
    displaySilhouetteAct->setCheckable(true);
    displaySilhouetteAct->setShortcut(tr("Ctrl+L"));
    displaySilhouetteAct->setStatusTip(tr("Display silhouette type"));
    connect(displaySilhouetteAct, SIGNAL(triggered()), this, SLOT(displaySilhouette()));

    displayColorsAct = new QAction(tr("&Colors"), this);
    displayColorsAct->setCheckable(true);
    displayColorsAct->setShortcut(tr("Ctrl+R"));
    displayColorsAct->setStatusTip(tr("Display colors type"));
    connect(displayColorsAct, SIGNAL(triggered()), this, SLOT(displayColors()));

    displayDotsAct = new QAction(tr("&Dots"), this);
    displayDotsAct->setCheckable(true);
    displayDotsAct->setShortcut(tr("Ctrl+J"));
    displayDotsAct->setStatusTip(tr("Display dots type"));
    connect(displayDotsAct, SIGNAL(triggered()), this, SLOT(displayDots()));

    displayFlowLinesAct = new QAction(tr("&Flowlines"), this);
    displayFlowLinesAct->setCheckable(true);
    displayFlowLinesAct->setShortcut(tr("Ctrl+F"));
    displayFlowLinesAct->setStatusTip(tr("Display flowlines"));
    connect(displayFlowLinesAct, SIGNAL(triggered()), this, SLOT(displayFlowLines()));

    triangularizeAct = new QAction(QIcon("images/triangul.png"),tr("&Triangularize"), this);
    triangularizeAct->setStatusTip(tr("Triangularize"));
    connect(triangularizeAct, SIGNAL(triggered()), this, SLOT(calcTriangulation()));

    openscadAct = new QAction(QIcon("images/scadexport.png"),tr("&export openSCAD"), this);
    openscadAct->setStatusTip(tr("Export openSCAD"));
    connect(openscadAct, SIGNAL(triggered()), this, SLOT(exportOpenscad()));

    reducetriAct = new QAction(QIcon("images/reducetri.png"),tr("&reduce Triangulation"), this);
    reducetriAct->setStatusTip(tr("Reduce triangulation"));
    connect(reducetriAct, SIGNAL(triggered()), this, SLOT(reduceTriangulation()));


    centerAct = new QAction(tr("&Center"), this);
    centerAct->setCheckable(true);
    centerAct->setShortcut(tr("Ctrl+E"));
    centerAct->setStatusTip(tr("Center the selected text"));
    connect(centerAct, SIGNAL(triggered()), this, SLOT(center()));

    alignmentGroup = new QActionGroup(this);
    alignmentGroup->addAction(displaySilhouetteAct);
    alignmentGroup->addAction(displayColorsAct);
    alignmentGroup->addAction(displayDotsAct);
    //alignmentGroup->addAction(centerAct);
    displaySilhouetteAct->setChecked(true);
}


void Window::createMenus()
 {
    fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(newAct);
    fileMenu->addAction(openAct);
    fileMenu->addAction(saveAct);
    fileMenu->addAction(printAct);
    fileMenu->addSeparator();
    fileMenu->addAction(exitAct);

    editMenu = menuBar()->addMenu(tr("&Edit"));
    editMenu->addAction(triangularizeAct);
    editMenu->addAction(openscadAct);
    editMenu->addAction(reducetriAct);

    //editMenu->addAction(undoAct);
    //editMenu->addAction(redoAct);
    //editMenu->addSeparator();
    //editMenu->addAction(cutAct);
    //editMenu->addAction(copyAct);
    //editMenu->addAction(pasteAct);
    //editMenu->addSeparator();

    viewMenu = menuBar()->addMenu(tr("&View"));
    //formatMenu->addAction(boldAct);
    //formatMenu->addAction(italicAct);
    viewMenu->addSeparator()->setText(tr("Display type"));
    viewMenu->addAction(displaySilhouetteAct);
    viewMenu->addAction(displayColorsAct);
    viewMenu->addAction(displayDotsAct);
    viewMenu->addSeparator();
    viewMenu->addAction(displayFlowLinesAct);

    //formatMenu->addAction(centerAct);

    viewMenu->addSeparator();

    //formatMenu->addAction(setLineSpacingAct);
    //formatMenu->addAction(setParagraphSpacingAct);

    helpMenu = menuBar()->addMenu(tr("&Help"));
    helpMenu->addAction(aboutAct);
    helpMenu->addAction(aboutQtAct);

 }

