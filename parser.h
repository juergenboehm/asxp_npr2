
#ifndef __parser_h
#define __parser_h

#include <string>

#include <assert.h>

#include "poly.h"

//#define TRACE


#define SYM_NIL 0

#define SYM_PLUS 10
#define SYM_MINUS 11
#define SYM_TIMES 12
#define SYM_DIVIDE 13
#define SYM_EXPO 14

#define SYM_LPAREN 20
#define SYM_RPAREN 21
#define SYM_SEMICOLON 30

#define SYM_INTEGER 40
#define SYM_DOUBLE 41
#define SYM_NUMBER 42

#define SYM_VARIABLE 43

#define SYM_EOF	50

#define SYM_ERR	100

#define ERR_NUM_AFTER_EXPO 10
#define ERR_NO_RPAREN 11
#define ERR_BASE_MALFORMED	12



class GetchClass {

public:

	virtual ~GetchClass(){};

	virtual char operator()(int & is_end) = 0;

	virtual void advance() = 0;

};


class GetchString: public GetchClass {

public:

		GetchString(std::string sa): s(sa), pos(0) {};
		virtual ~GetchString() {};

		virtual char operator()(int & is_end) {
			if (pos == s.length()) {
				is_end = 1;
				return 255;
			}
			is_end = 0;
			return s[pos];
		}

		virtual void advance() {
			if (pos < s.length()) {
				++pos;
			}
		};

	private:

		std::string s;
		size_t pos;


};

class Token {

public:

	Token(): sym_code(0), ident(""), vald(0.0), vali(0) {};

	Token(int sym_codea, std::string identa, double valda, int valia):
		sym_code(sym_codea), ident(identa), vald(valda), vali(valia) {};

	Token(int sym_codea, std::string identa):
		sym_code(sym_codea), ident(identa), vald(0.0), vali(0) {};

	int sym_code;

	std::string ident;
	double vald;
	int vali;

};

template<class Tok>
class Lexer {

public:

	typedef Tok (*token_fun)(std::string);

	struct lex_table_line_s {

		int state;
		std::string init_str;
		int new_state;
		int action;

		token_fun tok_fun;

	};

	typedef struct lex_table_line_s lex_table_line_t;


	Lexer(vector<lex_table_line_t> & lex_tbla): lex_tbl(lex_tbla), pgetch(0), lex_state(0), text_buf("") {};

	void init(GetchClass & getcha, int lex_inita) {
		pgetch = &getcha;
		lex_state = lex_inita;
		text_buf = "";
	}

	Tok getsym();

private:

	vector<lex_table_line_t> & lex_tbl;
	GetchClass* pgetch;

	int lex_state;

	std::string text_buf;


};


void parser_test_fun(std::string line);
void lexer_test_fun(std::string line);

// grammar is
//
// expression = [+/-] term ( +|- term)*
// term = factor ( * factor)*
// factor = base [^number]
// base = variable | (expression) | number
//


template<class T>
class ParserExpr {

public:

	ParserExpr(Lexer<Token> & lexa): lex(lexa) {};

	int parse_init() { akt_sym = lex.getsym(); return 0; }

	int parse_expression(T & pvalue);
	int parse_term(T & value);
	int parse_factor(T & pvalue);
	int parse_base(T & value);

	Lexer<Token>& lex;

	Token akt_sym;

};

#ifdef TEST_PARSER

class DoubleExpr {

public:

	DoubleExpr(): val(0.0) {};

	DoubleExpr & pre_init(DoubleExpr & a) { return *this; };

	DoubleExpr & set_ident(std::string s ) {
		val = 0;
		if (s == "x") {
			val = 2;
		}
		if (s == "y") {
			val = 5;
		}
		return *this;
	};

	DoubleExpr & set_double(double vala) {
		val = vala;
		return *this;
	}

	DoubleExpr & add(DoubleExpr & a) {
		val += a.val;
		return *this;
	};

	DoubleExpr & mul(DoubleExpr & a) {
		val *= a.val;
		return *this;
	};

	DoubleExpr & pow(int d) {
		std::cout << "val = " << val << " pow: d =" << d << std::endl;
		val = exp(d * log(val));
		return *this;
	}

	DoubleExpr & neg() {
		val = -val;
		return *this;
	}

	double val;

};

#endif

class PolyExpr {

public:

	PolyExpr(): vars(), poly() {};

	PolyExpr & set_ident(std::string idstr) {
		for(size_t i = 0; i < vars.size(); ++i) {
			if (idstr == vars[i]) {
				poly = Poly<double>(1, i, 1, vars.size());
				return *this;
			}
		}
		poly = Poly<double>(1, 0, 0, vars.size());
		return *this;
	}

	PolyExpr& set_double(double val) {
		poly = Poly<double>(val, 0, 0, vars.size());
		return *this;
	}

	PolyExpr & pre_init(PolyExpr & a) {
		vars = a.vars;
		return *this;
	}

	void set_vars(std::string varstr) {
		for(size_t i = 0; i < varstr.size(); ++i) {
			vars.push_back(std::string(1,varstr[i]));
		}
	}

	PolyExpr & add(PolyExpr & a) {
		poly.add(a.poly);
		return *this;
	}

	PolyExpr & mul(PolyExpr & a) {
		poly.mul(a.poly);
		return *this;
	}

	PolyExpr & pow(int d) {
		poly.pow(d);
		return *this;
	}

	PolyExpr & neg() {
		poly.mul(Poly<double>(-1, 0, 0, vars.size()));
		return *this;
	}


	vector<std::string> vars;

	Poly<double> poly;
};


extern vector<Lexer<Token>::lex_table_line_t> ltbl;

template< class C>
int read_poly(char* poly_str, Poly<C> & perg)
{
	GetchClass* p_getch = new GetchString(std::string(poly_str));
	Lexer<Token> lex(ltbl);

	ParserExpr<PolyExpr> parse_poly(lex);

	lex.init(*p_getch, 0);

	PolyExpr res;

	std::string varstr = "xyzabcdefg";

	varstr = varstr.substr(0, perg.nvars);

	res.set_vars(varstr);

	parse_poly.parse_init();
	int err = parse_poly.parse_expression(res);

	assert(err == 0);

	if (!err) {
		perg = res.poly;
	};

	delete p_getch;

	return 0;
}





#endif
