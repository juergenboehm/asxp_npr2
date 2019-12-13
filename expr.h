#ifndef __expr_h_
#define __expr_h_

#include <memory>
#include <string>
#include <map>
#include <iostream>

#include <poly.h>


#define EXPR_NODE	0
#define EXPR_SYMBOL	1
#define EXPR_VAL	2


//#define DC(x, y) (std::cout << "Destructor called:" << x << y << std::endl)

#define DC(x, y) if (_alloc_cnt == 0) { std::cout << "Destructor called: alloc_cnt = " << _alloc_cnt << std::endl; }

class Expr;
class ExprVal;
class ExprSymbol;
class ExprNode;

class Ex {

public:

	Ex();

	Ex(const Ex & a);

	Ex & operator=(const Ex & a);

	~Ex();

	Ex(Expr* a);
	Ex(ExprVal* a);
	Ex(ExprSymbol * a);
	Ex(ExprNode * a);

	Ex(double val);
	Ex(std::string str);
	Ex(Ex a, Ex b);


	Ex set_ident(std::string idstr);

	Ex set_double(double val);

	Ex pre_init(Ex a);

	Ex add(Ex a);

	Ex mul(Ex a);

	Ex pow(int d);

	Ex neg();

	std::string to_str() const;

	std::string get_name () const;

	double get_val() const;

	Ex rplaca(Ex e);
	Ex rplacd(Ex e);

	Ex clone();

	Expr* operator->() {return _impl;}

	Expr* _impl;

};


class Expr {

public:

	Expr(int typea): type(typea), ref_cnt(0) { ++_alloc_cnt;
		//std::cout << "Expr(int): type = " << type << " ref_cnt = " << ref_cnt << std::endl;
	};

	virtual ~Expr() {--_alloc_cnt; DC("alloc cnt ", _alloc_cnt);};

	virtual std::string to_str() const = 0;

	virtual Expr* clone() = 0;

	virtual Ex get_car() = 0;
	virtual Ex get_cdr() = 0;

	virtual std::string get_name() const = 0;

	virtual double get_val() const = 0;

	virtual Ex rplaca(Ex e) = 0;
	virtual Ex rplacd(Ex e) = 0;


	int ref_cnt;
	int type;

static int _alloc_cnt;

};



class ExprVal: public Expr {

public:

	ExprVal(double vala): Expr(EXPR_VAL), val(vala) {
		//std::cout << "ExprVal(double): val = " << val << std::endl;
	};

	virtual ~ExprVal() { DC("val = ", val); };

	virtual std::string to_str() const;

	virtual ExprVal* clone();

	virtual Ex get_car();

	virtual Ex get_cdr();

	virtual std::string get_name() const;

	virtual double get_val() const;

	virtual Ex rplaca(Ex e);
	virtual Ex rplacd(Ex e);




	double val;
};

class ExprSymbol: public Expr {

public:

	ExprSymbol(const std::string stra): Expr(EXPR_SYMBOL), str(stra) {
		//std::cout << "ExprSymbol(string): str = " << str << std::endl;
	};

	virtual ~ExprSymbol() { DC("str = ", str);};

	virtual std::string to_str() const;

	virtual ExprSymbol* clone();

	virtual Ex get_car();

	virtual Ex get_cdr();

	virtual std::string get_name() const;

	virtual double get_val() const;

	virtual Ex rplaca(Ex e);
	virtual Ex rplacd(Ex e);



	std::string str;
};


class ExprNode: public Expr {

public:

	ExprNode(Ex a, Ex b): Expr(EXPR_NODE), car(a), cdr(b) {
		//std::cout << "ExprNode(Ex, Ex): ( " << car.to_str() << ". " << cdr.to_str() << " )" << std::endl;
	};

	virtual ~ExprNode() { DC("Node", ""); };

	virtual std::string to_str() const;

	virtual ExprNode* clone();

	virtual Ex get_car();

	virtual Ex get_cdr();

	virtual std::string get_name() const;

	virtual double get_val() const;

	virtual Ex rplaca(Ex e);
	virtual Ex rplacd(Ex e);


	Ex car, cdr;

};




class FactExpr {

public:

	static Ex create(const std::string & str);

	static Ex create(double vala);

	static Ex create(Ex a, Ex b);

	static Ex create(Ex a);

};

class ExEnviron {

public:

	Ex lookup(std::string str) const;

	Ex enter(std::string str, Ex val);

	std::map<std::string, Ex> _env_map;

};

Ex evalsimp(Ex e, const ExEnviron & env, bool deep);



std::ostream & operator<<(std::ostream & os, const Ex & a);

Ex car(Ex e);
Ex cdr(Ex e);

bool is_atom(Ex e);
bool is_val(Ex e);
bool is_symbol(Ex e);
bool is_nil(Ex e);

std::string to_str_infix(Ex e);


Ex operator-(Ex a);

Ex operator+(Ex a, Ex b);
Ex operator-(Ex a, Ex b);
Ex operator*(Ex a, Ex b);
Ex operator/(Ex a, Ex b);

Ex pow(Ex a, Ex b);



int read_expr(std::string expr_str, Ex & perg);
int read_module(std::string expr_str, Ex & perg);

Ex eval_statement_list(Ex & statlis, Poly<double> & retpoly);



void test_expr_1(std::string expr_str);

extern Ex nilEx;



#endif /* __expr_h_ */
