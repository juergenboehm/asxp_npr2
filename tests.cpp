
#include <iostream>
#include <fstream>
#include <iomanip>

#include <stdlib.h>

#include "roots.h"
#include "poly.h"
#include "parser.h"

#include "expr.h"

#include "tests.h"

using namespace std;

void test_unit_1()
{

	char buf[4096];
	int code = 0;

	Poly5 f(5);
	Poly5 g(5);
	Poly5 f1(5);
	Poly5 f2(5);
	Poly5 f3(5);

	fstream tstfile("test.data");

	while (! tstfile.eof()) {
		tstfile.getline(buf, 4096);
		//cout << buf << endl;

		if (buf[0] == 'L') {
			code = atoi(buf + 1);
			tstfile.getline(buf, 4096);
			if (tstfile.eof()) {
				break;
			}
		}

		if (code == 1) {
			f.read(buf);
		}
		if( code == 2) {
			g.read(buf);
		}
		if (code == 3) { //  f + g
			f1.read(buf);
			//cout << "f1 = " << f1 << endl;
			Poly5 f1c = f1;
			Poly5 fc = f;
			fc.add(g);
			f1c.mul(Poly5(-1,0,0, 5));
			f1c.add(fc);
			cout << "f11 = " << f1c << endl;
		}
		if (code == 4) { // f * g
			f2.read(buf);
			Poly5 f2c = f2;
			Poly5 fc = f;
			fc.mul(g);
			f2c.mul(Poly5(-1,0,0, 5));
			f2c.add(fc);
			cout << "f21 = " << f2c << endl;
		}
		if (code == 5) { // f^3 + g^2
			f3.read(buf);
			Poly5 f3c = f3;
			Poly5 fc = f;
			fc.pow(3);
			Poly5 gc = g;
			gc.pow(2);
			fc.add(gc);
			fc.mul(Poly5(-1,0,0, 5));
			fc.add(f3c);
			cout << "fc = " << fc << endl;
		}

		code = 0;

	}

	cout << "press key." << endl;
	getchar();
}

void test_unit_2()
{
	Poly5 p4(5);
	Poly5 p5(5);
	Poly5 p6(5);

	p4.read("(x+y)^5*z^2*a^3*b^4");

    Poly4 p14(4);
    Poly3 p13(3);
    Poly2 p12(2);
    Poly1 p11(1);

    p4.eval_last(1, p14);
    p14.eval_last(1, p13);
    p13.eval_last(1, p12);
    p12.eval_last(1, p11);

    cout << "(x+1)^5 = " << p11 << endl;

    p4 = Poly5(0,0,0, 5);

    int ret = p4.read("x^3 + y^5 + z^4");

    cout << "ret = " << ret << endl;
    cout << "x^3 + y^5 + z^4 " << p4 << endl;

    cout << "press key." << endl;
    getchar();

    p4 = Poly5(0,0,0, 5);

    ret = p4.read("(a+b)*(x+y)");

    cout << "ret = " << ret << endl;
    cout << "(a+b)*(x+y) = " << p4 << endl;

    cout << "press key." << endl;
    getchar();


	ret = p5.read(
			"-6*x*y*z -3*x^2 -6*x*y -6*x*z -3*y^2 -6 * y * z -3 * z^2 -3 * x^2 * y -3 * x^2 *z"
				"-3 * x* y^2 -3 * x * z^2 -3 * y^2 *z - 3 * y * z^2 - 3 * x -3 * y - 3 * z"
				"");

    cout << "ret = " << ret << endl;
    cout << "clebsch cubic 1 =  " << p5 << endl;

    getchar();

    p6 = Poly5(0,0,0, 5);

    ret = p6.read("x^3+y^3+z^3 + 1 + (-1)*(x+y+z+1)^3");

    cout << "clebsch cubic 2 = " << p6 << endl;

    p6.mul(Poly5(-1,0,0, 5));
    p6.add(p5);

    cout << "difference = " << p6 << endl;

    getchar();

    p6.read("x + a*y + b*z");

    cout << "x + a * y + b * z = " << p6 << endl;

    double substmat[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };

    p6.subst_mat(3, substmat);

    cout << "p6 substituted = " << p6 << endl;

    getchar();

}

void test_unit_3(std::string line)
{
	cout << "Test Lexer:" << endl;

	lexer_test_fun(line);

	cout << "Press key." << endl;

	getchar();


	cout << "Test Parser:" << endl;

	parser_test_fun(line);

	cout << "Press key." << endl;

	getchar();
}

void test_unit(std::string line)
{

	if (1)
	{
		test_expr_1(line);

		cout << "Press key" << endl;

		getchar();
	}


	test_unit_3(line);

	test_unit_1();
	test_unit_2();

    cout << "test root_list" << endl;

    test_root_list();

    cout << "press key." << endl;

    getchar();

}
