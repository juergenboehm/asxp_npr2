
#ifndef __parser_h_
#define __parser_h_

#include <string>
#include <memory>

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
#define SYM_COMMA 31
#define SYM_EQUAL 32
#define SYM_LBRACKET 33
#define SYM_RBRACKET 34
#define SYM_PERIOD 35

#define SYM_INTEGER 40
#define SYM_DOUBLE 41
#define SYM_NUMBER 42

#define SYM_VARIABLE 43

#define SYM_EOF	50

#define SYM_ERR	100

#define ERR_NUM_AFTER_EXPO 10
#define ERR_NO_RPAREN 11
#define ERR_BASE_MALFORMED	12
#define ERR_SYMBOL_EXPECTED	13
#define ERR_EQUALSIGN_EXPECTED	14
#define ERR_NO_LPAREN	15
#define ERR_NO_RBRACKET	16
#define ERR_BAD_STATEMENT_DELIM 17



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
			//std::cout << s[pos];
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

template< class T>
int ParserExpr<T>::parse_expression(T & pvalue)
{

	T ptemp;

	ptemp.pre_init(pvalue);

	bool first_minus = false;

	if (akt_sym.sym_code == SYM_PLUS) {
		akt_sym = lex.getsym();
	} else if (akt_sym.sym_code == SYM_MINUS) {
		akt_sym = lex.getsym();
		first_minus = true;
	};

	int err = parse_term(ptemp);

	if (err) {
		return err;
	}

	if (first_minus) {
		ptemp.neg();
	}


	do {

		if (akt_sym.sym_code == SYM_PLUS || akt_sym.sym_code == SYM_MINUS) {

			bool is_minus = (akt_sym.sym_code == SYM_MINUS);

			T ptemp1;

			ptemp1.pre_init(pvalue);

			akt_sym = lex.getsym();

			int err = parse_term(ptemp1);

			if (err) {
				return err;
			}

			if (is_minus) {
				ptemp1.neg();
			}

			ptemp.add(ptemp1);

		} else {

			pvalue = ptemp;

			return 0;
		}

	} while (1);

	pvalue = ptemp;

	return 0;
}


template<class T>
int ParserExpr<T>::parse_term(T & pvalue)
{

	T ptemp;

	ptemp.pre_init(pvalue);

	int err = parse_factor(ptemp);

	if (err) {
		return err;
	}

	do {

		if (akt_sym.sym_code == SYM_TIMES) {

			akt_sym = lex.getsym();

			T ptemp1;

			ptemp1.pre_init(pvalue);

			int err = parse_factor(ptemp1);

			if (err) {
				return err;
			}

			ptemp.mul(ptemp1);

		} else {
			pvalue = ptemp;
			return 0;
		}

	} while (1);

	pvalue = ptemp;

	return 0;
}

template<class T>
int ParserExpr<T>::parse_factor(T & pvalue)
{

	T ptemp;

	ptemp.pre_init(pvalue);

	int err = parse_base(ptemp);

	if (err) {
		return err;
	}

	if (akt_sym.sym_code == SYM_EXPO) {

		akt_sym = lex.getsym();

		if (akt_sym.sym_code == SYM_INTEGER) {

			ptemp.pow(akt_sym.vali);

			akt_sym = lex.getsym();


			pvalue = ptemp;
			return 0;

		} else {

			assert(0);
			return ERR_NUM_AFTER_EXPO;
		}
	} else {

		pvalue = ptemp;

		return 0;

	}
}



template<class T>
int ParserExpr<T>::parse_base(T & pvalue)
{


	if (akt_sym.sym_code == SYM_INTEGER || akt_sym.sym_code == SYM_DOUBLE) {

		T pval;
		pval.pre_init(pvalue);
		pval.set_double(akt_sym.vald);
		pvalue = pval;

		akt_sym = lex.getsym();

		return 0;

	} else if (akt_sym.sym_code == SYM_LPAREN) {

		akt_sym = lex.getsym();

		int err = parse_expression(pvalue);

		if (err) {
			return err;
		}

		if (akt_sym.sym_code == SYM_RPAREN) {

			akt_sym = lex.getsym();

			return 0;

		} else {

			assert(0);
			return ERR_NO_RPAREN;
		}
	} else if (akt_sym.sym_code == SYM_VARIABLE) {

		T pval;
		pval.pre_init(pvalue);
		pval.set_ident(akt_sym.ident);
		pvalue = pval;

		akt_sym = lex.getsym();

		return 0;

	} else {
		assert(0);
		return ERR_BASE_MALFORMED;
	}
	return 0;
}


//#define TEST_PARSER

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
int read_poly(std::string poly_str, Poly<C> & perg)
{
	std::unique_ptr<GetchClass> p_getch(new GetchString(std::string(poly_str)));
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

	return 0;
}





#endif
