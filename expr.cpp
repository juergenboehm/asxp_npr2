
#include <map>
#include <functional>
#include <sstream>

#include "poly.h"
#include "parserex.h"
#include "parser.h"
#include "expr.h"

#define DBGOUT(x) (std::cout << #x << x << std::endl)

Ex::Ex(): _impl(0) {

}


Ex::Ex(const Ex & a) {
#if 0
	std::cout << "Ex::Ex(Ex): a = " << a.to_str() << std::endl;
#endif

	_impl = a._impl;
	if (_impl) {
		++(_impl->ref_cnt);
	}
}

Ex & Ex::operator=(const Ex & a) {

#if 0
	if (_impl) {
		std::cout << "Ex::operator=:  this =  " << _impl->to_str() << " a = " << a.to_str() << std::endl;
	} else {
		std::cout << "Ex::operator=: _impl = " << _impl << " a = " << a.to_str() << std::endl;
	}
#endif

	if (a._impl == _impl) {
		return *this;
	}

	if (_impl) {
		--(_impl->ref_cnt);
		if (! _impl->ref_cnt) {
			delete _impl;
		}
	}
	_impl = a._impl;
	if (_impl) {
		++(_impl->ref_cnt);
	}
	return *this;
}

Ex::~Ex() {
	if (_impl) {
		--(_impl->ref_cnt);
		if (! (_impl->ref_cnt)) {
			delete _impl;
		}
	}
}


Ex::Ex(Expr* a): _impl(a) {
	++(_impl->ref_cnt);
	//std::cout << "Ex:Ex(Expr*):" << _impl << " rc = " << _impl->ref_cnt << std::endl;

}

Ex::Ex(ExprVal* a): _impl(a) {
	++(_impl->ref_cnt);
	//std::cout << "Ex:Ex(ExprVal*):" << _impl << " rc = " << _impl->ref_cnt << std::endl;
}

Ex::Ex(ExprSymbol * a): _impl(a) {
	++(_impl->ref_cnt);
	//std::cout << "Ex:Ex(ExprSymbol*):" << _impl << " rc = " << _impl->ref_cnt << std::endl;
}

Ex::Ex(ExprNode * a): _impl(a) {
	++(_impl->ref_cnt);
	//std::cout << "Ex:Ex(ExprNode*):" << _impl << " rc = " << _impl->ref_cnt << std::endl;
}

Ex::Ex(double val): _impl(0) {
	*this = FactExpr::create(val);
}

Ex::Ex(std::string str): _impl(0) {
	*this = FactExpr::create(str);
}

Ex::Ex(Ex a, Ex b):_impl(0) {
	*this = FactExpr::create(a, b);
}



Ex Ex::set_ident(std::string idstr) {

	*this = FactExpr::create(idstr);
	return *this;
}

Ex Ex::set_double(double val) {

	*this = FactExpr::create(val);
	return *this;
}

Ex Ex::pre_init(Ex a) {

	//return FactExpr::create(a);
	return *this;
}


Ex Ex::add(Ex a) {

	Ex res = *this + a;
	*this = res;

	return *this;
}

Ex Ex::mul(Ex a) {

	Ex res = *this * a;
	*this = res;

	return *this;
}

Ex Ex::pow(int d) {

	Ex res = ::pow(*this, FactExpr::create(d));
	*this = res;

	return *this;
}

Ex Ex::neg() {

	Ex res = -(*this);
	*this = res;

	return *this;
}

std::string Ex::to_str() const {
	return _impl->to_str();
}

std::string Ex::get_name () const {
	return _impl->get_name();
}

double Ex::get_val() const {
	return _impl->get_val();
}

Ex Ex::rplaca(Ex e)
{
	return _impl->rplaca(e);
}

Ex Ex::rplacd(Ex e)
{
	return _impl->rplacd(e);
}


Ex Ex::clone() {
	Ex res(_impl->clone());
	return res;
}




Ex nilEx = FactExpr::create("nil");


Ex build_unop(std::string unop, Ex a) {

	Ex res = Ex(Ex(unop), Ex(a, nilEx));

	return res;
}


Ex build_op(std::string op, Ex a, Ex b) {

	Ex res = Ex(Ex(op), Ex(a, Ex(b, nilEx)));

	return res;
}

Ex operator-(Ex a) {
	return build_unop("-", a);
}


Ex operator+(Ex a, Ex b) {

	return build_op("+", a, b);
}

Ex operator-(Ex a, Ex b) {

	return build_op("-", a, b);
}

Ex operator*(Ex a, Ex b) {

	return build_op("*", a, b);
}

Ex operator/(Ex a, Ex b) {

	return build_op("/", a, b);
}

Ex pow(Ex a, Ex b) {

	return build_op("^", a, b);
}


int Expr::_alloc_cnt = 0;

ExprVal* ExprVal::clone() {
	return new ExprVal(val);
}

ExprSymbol* ExprSymbol::clone() {
	return new ExprSymbol(str);
}

ExprNode* ExprNode::clone() {
	return new ExprNode(car.clone(), cdr.clone());
}

Ex ExprVal::get_car() {
	return nilEx;
}

Ex ExprVal::get_cdr() {
	return nilEx;
}

Ex ExprSymbol::get_car() {
	return nilEx;
}

Ex ExprSymbol::get_cdr() {
	return nilEx;
}

Ex ExprNode::get_car() {
	return car;
}

Ex ExprNode::get_cdr() {
	return cdr;
}


Ex ExprVal::rplaca(Ex e) {
	return nilEx;
}

Ex ExprVal::rplacd(Ex e) {
	return nilEx;
}

Ex ExprSymbol::rplaca(Ex e) {
	return nilEx;
}

Ex ExprSymbol::rplacd(Ex e) {
	return nilEx;
}

Ex ExprNode::rplaca(Ex e) {
	car = e;
	return Ex(this);
}

Ex ExprNode::rplacd(Ex e) {
	cdr = e;
	return Ex(this);
}






std::string ExprVal::get_name() const {
	return "";
}


std::string ExprSymbol::get_name() const {
	return str;
}


std::string ExprNode::get_name() const {
	return "";
}


double ExprVal::get_val() const {
	return val;
}


double ExprSymbol::get_val() const {
	return 0.0;
}


double ExprNode::get_val() const {
	return 0.0;
}







Ex FactExpr::create(const std::string & str) {
	return Ex((new ExprSymbol(str)));
}

Ex FactExpr::create(double vala) {
	return Ex((new ExprVal(vala)));
}

Ex FactExpr::create(Ex a, Ex b) {
	return Ex((new ExprNode(a, b)));
}


Ex FactExpr::create(Ex a) {
	switch (a->type) {
	case EXPR_SYMBOL: return create(dynamic_cast<ExprSymbol*>(a._impl)->str);
	case EXPR_VAL:	return create(dynamic_cast<ExprVal*>(a._impl)->val);
	case EXPR_NODE: return create((dynamic_cast<ExprNode*>(a._impl)->car),
										(dynamic_cast<ExprNode*>(a._impl)->cdr));
	}
}


Ex ExEnviron::lookup(std::string str) const
{
	auto it = _env_map.find(str);

	if (it == _env_map.end()) {
		return nilEx;
	}

	return it->second;
}

Ex ExEnviron::enter(std::string str, Ex val)
{
	Ex res = lookup(str);

	_env_map[str] = val;

	return res;

}


Ex car(Ex e) {
	return e->get_car();
}

Ex cdr(Ex e) {
	return e->get_cdr();
}

bool is_atom(Ex e) {
	return is_val(e) || is_symbol(e);
}

bool is_val(Ex e) {
	return e->type == EXPR_VAL;
}

bool is_symbol(Ex e) {
	return e->type == EXPR_SYMBOL;
}

bool is_nil(Ex e) {
	return is_symbol(e) && e.get_name() == "nil";
}

std::map<string, bool> is_infi_map = {{"+", true}, {"-", true},
									{"*", true}, {"/", true}, {"^", true}, {"=", true}};

bool is_infix(Ex e) {
	return  is_infi_map.find(e.get_name()) != is_infi_map.end();
}

Ex get_op(Ex e) {
	return car(e);
}

int length(Ex e) {
	int cnt = 0;
	while (! is_nil(e)) {
		e = cdr(e);
		++cnt;
	}
	return cnt;
}

Ex get_operand_1(Ex e) {
	return car(cdr(e));
}

Ex get_operand_2(Ex e) {
	return car(cdr(cdr(e)));
}

int nops(Ex e) {
	return length(cdr(e));
}


std::string to_str_infix(Ex e) {

	if (is_val(e)) {
		return e.to_str();
	} else if (is_symbol(e)) {
		return e.to_str();
	}

	Ex op = car(e);

	if (is_infix(op)) {
		Ex x = car(cdr(e));
		Ex y = cdr(cdr(e));

		if (is_nil(y)) {
			std::string res = to_str_infix(op) + "(" + to_str_infix(x) + ")";
			return res;
		} else {
			y = car(y);
			std::string res = "(" + to_str_infix(x) + ")" + to_str_infix(op) + "(" + to_str_infix(y) + ")";
			return res;
		}

	} else {
		std::string res = to_str_infix(op) + "(";
		Ex argl = cdr(e);
		bool is_first = true;
		while (! is_nil(argl)) {
			if (! is_first) {
				res = res + ",";
			} else {
				is_first = false;
			}
			Ex x = car(argl);
			res = res + to_str_infix(x);
			argl = cdr(argl);
		}
		res = res + ")";
		return res;
	}
	return "";
}


double do_op_2(Ex op, double x, double y)
{
	if (op.get_name() == "+") {
		return x + y;
	} else if (op.get_name() == "-") {
		return x - y;
	} else if (op.get_name() == "*") {
		return x * y;
	} else if (op.get_name() == "/") {
		return x/y;
	} else if (op.get_name() == "^") {
		return ::pow(x,y);
	} else {
		assert(0);
	}
	return 0.0;
}

double do_op_1(Ex op, double x)
{
	if (op.get_name() == "-") {
		return -x;
	} else if (op.get_name() == "sqrt"){
		return sqrt(x);
	} else {
		assert(0);
	}
	return 0.0;
}


Ex evalsimp_op(Ex op, Ex x, Ex y)
{
	if (op.get_name() == "/" && is_val(y) && !is_val(x)) {
		return evalsimp_op(Ex("*"), x, 1/y);
	}
	if (is_val(x) && is_val(y)) {
		return Ex(do_op_2(op, x.get_val(), y.get_val()));
	}
	return Ex(op, Ex(x, Ex(y, nilEx)));
}

Ex evalsimp_op(Ex op, Ex x)
{
	if (is_val(x)) {
		return Ex(do_op_1(op, x.get_val()));
	}
	return Ex(op, Ex(x, nilEx));
}

Ex append(Ex lis, Ex elem)
{
	if (is_nil(lis)) {
		return Ex(elem, nilEx);
	}

	return Ex(car(lis), append(cdr(lis), elem));
}

Ex reverse(Ex lis)
{
	if (is_nil(lis)) {
		return lis;
	}
	return append(reverse(cdr(lis)), car(lis));
}

Ex map_lis(std::function<Ex(Ex)> f, Ex lis)
{

	if (is_nil(lis)) {
		return lis;
	}

	if (is_atom(lis)) {
		return f(lis);
	}

	return Ex(f(car(lis)), map_lis(f, cdr(lis)));
}

Ex evalsubs(Ex e, const ExEnviron & env)
{
	Ex sl = get_operand_1(e);

	if (get_op(sl).get_name() != "subslis" ) {
		return nilEx;
	}

	Ex ee = get_operand_2(e);
	ee = evalsimp(ee, env, true);


	ExEnviron env_aux;

	map_lis([&](Ex ass)->Ex
			{Ex v = get_operand_1(ass);
			 Ex val = get_operand_2(ass);
			 val = evalsimp(val, env, true);
			 env_aux.enter(v.get_name(), val);
			 return nilEx;
			}, sl);

	ee = evalsimp(ee, env_aux, false);

	ee = evalsimp(ee, env, true);

	return ee;

}

Ex evalsimp(Ex e, const ExEnviron & env, bool deep)
{

	//std::cout << "evalsimp: e = " << e << std::endl;

	if (is_nil(e)) {
		return e;
	}

	if (is_val(e)) {
		return e;
	}

	if (is_symbol(e)) {
		Ex res = env.lookup(e.get_name());
		if (is_nil(res)) {
			return e;
		}
		if (deep) {
			res = evalsimp(res, env, deep);
		};
		return res;
	}

	Ex op = car(e);

	int nargs = nops(e);

	if (nargs == 2) {

		Ex x = get_operand_1(e);
		Ex x1;

		Ex y = get_operand_2(e);
		y = evalsimp(y, env, deep).clone();

		if (op.get_name() == "subs") {

			return evalsubs(e, env);

		} else if (op.get_name() == "=") {
			x1 = x;
		} else {
			x1 = evalsimp(x, env, deep).clone();
		}

		return evalsimp_op(op, x1, y);

	} else if (nargs == 1) {

		Ex x = get_operand_1(e);
		x = evalsimp(x, env, deep).clone();

		return evalsimp_op(op, x);
	} else {
		return map_lis([&](Ex e)->Ex{ return evalsimp(e, env, deep);}, e);
	}

}

int ex_to_poly(Ex e, Poly<double> & pol)
{
	if (is_symbol(e)) {
		std::string v = e.get_name();
		read_poly(v, pol);
		return 0;
	}

	if (is_val(e)) {
		Poly<double> pol1(e.get_val(), 0, 0, pol.nvars);
		pol = pol1;
	}

	Ex op = get_op(e);

	int nargs = nops(e);

	Ex e1 = get_operand_1(e);

	Poly<double> pol1(pol.nvars);

	ex_to_poly(e1, pol1);

	if (nargs == 2) {

		Ex e2 = get_operand_2(e);

		if (op.get_name() == "^") {

			if (is_val(e2)) {
				pol1.pow(e2.get_val());
			}

			pol = pol1;

			return 0;
		}

		Poly<double> pol2(pol.nvars);

		ex_to_poly(e2, pol2);

		if (op.get_name() == "+") {
			pol1.add(pol2);
			pol = pol1;
			return 0;
		} else if (op.get_name() == "-") {
			pol1.add(pol2.mul(Poly<double>(-1, 0, 0, pol.nvars)));
			pol = pol1;
			return 0;
		} else if (op.get_name() == "*") {
			pol1.mul(pol2);
			pol = pol1;
			return 0;
		} else {
			assert(0);
		}

	} else if (nargs == 1) {
		if (op.get_name() == "-") {
			pol1.mul(Poly<double>(-1, 0, 0, pol.nvars));
			pol = pol1;
			return 0;
		} else {
			assert(0);
		}
	}
}

Ex eval_statement(Ex & stat, ExEnviron & env, Poly<double> & retpoly)
{
	std::cout << "eval_statement: stat = "  << to_str_infix(stat) << std::endl;

	Ex op = get_op(stat);

	if (op.get_name() == "poly") {

		Ex e = get_operand_1(stat);

		e = evalsimp(e, env, true);

		Poly<double> e_poly(retpoly.nvars);

		ex_to_poly(e, e_poly);

		std::cout << "e_poly = " << e_poly.to_str() << std::endl;

		retpoly = e_poly;

	} else if (op.get_name() == "print") {

		map_lis([&](Ex e)->Ex
				{
					std::string res;

					if (!is_symbol(e)) {
						return nilEx;
					}
					if (is_symbol(e)) {
						res += (e.get_name() + ": ");
					}

					res = res + to_str_infix(evalsimp(e, env, true));

					std::cout << res << std::endl;
					return nilEx;

				}, reverse(cdr(stat)));


	} else if (op.get_name() == "=") {
		Ex v = get_operand_1(stat);
		Ex e = get_operand_2(stat);

		Ex e1 = evalsimp(e, env, true);

		env.enter(v.get_name(), e1);
	}
	return nilEx;
}


Ex eval_statement_list(Ex & statlis, Poly<double> & retpoly)
{
	ExEnviron env;

	if (car(statlis).get_name() != "statlis" ) {
		return nilEx;
	}

	Ex slis = cdr(statlis);

	Ex result = nilEx;

	while (!is_nil(slis)) {

		Ex stat1 = car(slis);

		eval_statement(stat1, env, retpoly);

		slis = cdr(slis);

	}

	return result;
}



ostream & operator<<(ostream & os, const Ex & a) {
	std::string str = a.to_str();
	os << str;
	return os;
}

std::string ExprVal::to_str() const {
	std::stringstream ss;
	ss << val;
	return ss.str();
}

std::string ExprSymbol::to_str() const {
	return str;
}

std::string ExprNode::to_str() const {
	std::string s = "( " + car.to_str() + "  "  + cdr.to_str() + " )";
	return s;
}



int read_expr(std::string expr_str, Ex & perg)
{
	std::unique_ptr<GetchClass> p_getch(new GetchString(expr_str));
	Lexer<Token> lex(ltbl);

	ParserEx parse_expr(lex);

	lex.init(*p_getch, 0);

	Ex res;

	parse_expr.parse_init();
	int err = parse_expr.parse_expression(res);

	assert(err == 0);

	//std::cout << "res = " << res << std::endl;

	if (!err) {
		perg = res;
	};

	return 0;
}


int read_module(std::string expr_str, Ex & perg)
{
	std::unique_ptr<GetchClass> p_getch(new GetchString(expr_str));
	Lexer<Token> lex(ltbl);

	ParserEx parse_expr(lex);

	lex.init(*p_getch, 0);

	Ex res;

	parse_expr.parse_init();
	int err = parse_expr.parse_module(res);

	std::cout << "err = " << err << std::endl;

	assert(err == 0);

	//std::cout << "res = " << res << std::endl;

	if (!err) {
		perg = res;
	};

	return 0;
}



std::string module_test_str =
		"f = \n\nx * y + c * z;\n"
		"g = x^2 * y^3 * a + b * x * (y + z) * (z^2+ x^3)*h(x+y+z^2, w+ 3, w + e);\n"
		"# comment\n"
		"print(f);\n"
		"poly(f,g);\n"
		"\n"
		"\n"
		"f = subs([a = 2, b = 3 + w], f);\n"
		"sl = [a = 2 + y, c = 3 + x^a];\n"
		"\n"
		"g = evalall(f);\n"
		"f = subs(sl, g);\n"
		"# comment 2\n"
		"\n"
		"\n"
		"poly(f,g).";


std::string sextic_str =

"x0 = 1/2 * (sqrt(5) - 1);"
""
"phi = x0/(1-x0);"
""
"print(phi);"
""
"f = 4*(phi^2*x^2-y^2)*(phi^2*y^2-z^2)*(phi^2*z^2-x^2)-(1+2*phi)*(x^2+y^2+z^2-1)^2;"
""
"print(f);"
""
"d = 3;"
""
"f1 = subs([x=x/d], f);"
""
"f2 = subs([x=x/d, y=y/d, z=z/d], f);"
""
"print(f2, f1);"
""
"poly(f2).";


void test_expr_1(std::string expr_str) {

Ex a = FactExpr::create(5.0);
Ex b = FactExpr::create("x");

Ex c = a * b;

std::cout << "c = " << c << std::endl;

Ex d;

read_expr(expr_str, d);

std::cout << "d = " << to_str_infix(d) << std::endl;

Ex e = d.clone();

std::cout << "e = " << to_str_infix(e) << std::endl;

read_expr("a * f(x,y,z+w) + g(h(s+t))", d);

std::cout << "d = " << to_str_infix(d) << std::endl;

read_expr("(1 + 2) * (4 + x) * (5 + 6 + y) + h(x, x+y, x+y + x^2)", d);

std::cout << "d = " << to_str_infix(d) << std::endl;

ExEnviron env1;

env1.enter("x", Ex(10));
env1.enter("y", Ex("x") + Ex(2));

Ex d1 = evalsimp(d, env1, true);

std::cout << "d1 = " << to_str_infix(d1) << std::endl;

d1 = Ex(Ex(1), Ex(2));

std::cout << "d1 = " << d1.to_str() << std::endl;

std::cout << "type d1 = " << d1->type << std::endl;

std::cout << "rplacdr 1: d1 = " <<  d1.rplacd(Ex("x")) << std::endl;

std::cout << "type d1 = " << d1->type << std::endl;

std::cout << "rplacdr 2: d1 = " << d1.to_str()  << std::endl;

Ex mod1;

read_module(sextic_str, mod1);

std::cout << "mod1:" << std::endl;

std::cout << to_str_infix(mod1) << std::endl;

Poly<double> retpoly(5);

assert(is_nil(eval_statement_list(mod1, retpoly)));

std::cout << "retpoly = " << retpoly.to_str() << std::endl;


}




