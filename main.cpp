
#include <QApplication>

#include <iostream>
#include <locale>

#include <string>

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

    print_binom();

    if (argc > 1) {

    	if (string(argv[1]) != "no") {

    		test_unit(string(argv[1]));

    	}
    } else {
    	//test_unit("(x+y) * (a + b)");
    }

    eigen_test();

    //cin.get();

 	Window graphicWin;
	graphicWin.show();
	return app.exec();
}

