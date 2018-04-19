
#include <QApplication>

#include <iostream>

#include "mainwindow.h"
#include "asxp.h"
#include "roots.h"
#include "eigen.h"

#include "triangularize.h"

#include "tests.h"

using namespace std;

int main(int argc, char* argv[])
{
	if (argc == 2) {
	} else {
	}
	phong_exponent = 32;

	QApplication app(argc, argv);

	//QApplication::setStyle("motif");

    prepare_binom();

    // cout << "prep binom done" << endl;

    // print_binom();

	//test_unit();

    eigen_test();

    //cin.get();

 	Window graphicWin;
	graphicWin.show();
	return app.exec();
}

