
#include <iostream>
#include <sstream>
#include <string>

#include "stdlib.h"
#include "string.h"

#include "parser.h"


#define INIT_STATE	0
#define DIGIT_STATE1	1
#define DIGIT_STATE2	2
#define DIGIT_STATE3	3
#define DIGIT_STATE4	4

#define COMMENT_STATE	6

#define SYMCHAR_STATE1	7

#define WHITESPACE_STATE 8

#define IDENT_STATE 9

//#define NO_ACT	0
#define ACCUM_ACT	(1 << 0)
#define RETURN_ACT	(1 << 1)
#define NUM_ERR_ACT (1 << 2)
#define ERRSYM_ACT	(1 << 3)
#define CONSUME_ACT	(1 << 4)


#define ALPHA_STR "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
#define DIGIT_STR	"0123456789"
#define SYMCHAR_STR "+-*/^()"


Token tokenize_ident(std::string idstr)
{
	return Token(SYM_VARIABLE, idstr, 0.0, 0);
}

Token tokenize_symchar(std::string idstr)
{
	int sym = 0;
	switch (idstr[0]) {
	case '+':	sym = SYM_PLUS;
				break;
	case '-': 	sym = SYM_MINUS;
				break;
	case '*': 	sym = SYM_TIMES;
				break;
	case '/':	sym = SYM_DIVIDE;
				break;
	case '^':	sym = SYM_EXPO;
				break;
	case '(':	sym = SYM_LPAREN;
				break;
	case ')':	sym = SYM_RPAREN;
				break;
	default:	sym = 0;
	}

	return Token(sym, idstr, 0.0, 0);

}

Token tokenize_int(std::string idstr)
{
	int vali = std::stoi(idstr);
	return Token(SYM_INTEGER, idstr, vali, vali);

}

Token tokenize_double(std::string idstr)
{
	std::istringstream is(idstr);
	double vald;
	is >> vald;
	return Token(SYM_DOUBLE, idstr, vald, 0);

}

Token tokenize_eof(std::string idstr)
{
	return Token(SYM_EOF, idstr, 0.0, 0);
}


vector<Lexer<Token>::lex_table_line_t> ltbl =
{
{ INIT_STATE, DIGIT_STR, DIGIT_STATE1, ACCUM_ACT, 0 },
{ INIT_STATE, "\xff", INIT_STATE, RETURN_ACT, tokenize_eof },
{ INIT_STATE, SYMCHAR_STR, INIT_STATE, RETURN_ACT, tokenize_symchar },
{ INIT_STATE, " \t\n", INIT_STATE, CONSUME_ACT, 0 },
{ INIT_STATE, "#", COMMENT_STATE, CONSUME_ACT, 0 },
{ INIT_STATE, ALPHA_STR, IDENT_STATE, ACCUM_ACT, 0 },
{ INIT_STATE, "", INIT_STATE, ERRSYM_ACT, 0 },
{ IDENT_STATE, ALPHA_STR, IDENT_STATE, ACCUM_ACT, 0 },
{ IDENT_STATE, DIGIT_STR, IDENT_STATE, ACCUM_ACT, 0 },
{ IDENT_STATE, "_", IDENT_STATE, ACCUM_ACT, 0 },
{ IDENT_STATE, "", INIT_STATE, RETURN_ACT, tokenize_ident },
{ COMMENT_STATE, "\n", INIT_STATE, CONSUME_ACT, 0 },
{ COMMENT_STATE, "", COMMENT_STATE, CONSUME_ACT, 0 },
{ DIGIT_STATE1, DIGIT_STR, DIGIT_STATE1, ACCUM_ACT, 0 },
{ DIGIT_STATE1, ".", DIGIT_STATE2, ACCUM_ACT, 0 },
{ DIGIT_STATE1, "eE", DIGIT_STATE3, ACCUM_ACT, 0 },
{ DIGIT_STATE1, "", INIT_STATE, RETURN_ACT, tokenize_int },
{ DIGIT_STATE2, DIGIT_STR, DIGIT_STATE2, ACCUM_ACT, 0 },
{ DIGIT_STATE2, "eE", DIGIT_STATE3, ACCUM_ACT, 0 },
{ DIGIT_STATE2, "", INIT_STATE, RETURN_ACT, tokenize_double },
{ DIGIT_STATE3, "+-", DIGIT_STATE4, ACCUM_ACT, 0 },
{ DIGIT_STATE3, DIGIT_STR, DIGIT_STATE4, ACCUM_ACT, 0 },
{ DIGIT_STATE3, "", INIT_STATE, NUM_ERR_ACT, 0 },
{ DIGIT_STATE4, DIGIT_STR, DIGIT_STATE4, ACCUM_ACT, 0 },
{ DIGIT_STATE4, "", INIT_STATE, RETURN_ACT, tokenize_double }
};

template<class Tok>
Tok Lexer<Tok>::getsym()
{

	while (1) {

		int is_end;
		char c = (*pgetch)(is_end);

		//std::cout << "c = " << (int)(c) << std::endl;

		size_t posi;

		for(posi = 0; posi < lex_tbl.size(); ++posi) {
			if (lex_tbl[posi].state == lex_state) {
				break;
			}
		}

		if (posi == lex_tbl.size()) {
			std::cout << "getsym: state not found." << std::endl;
			return Tok(SYM_ERR, "err: state not found.");
		}

		bool ok = false;

		do {

			bool is_in_init_str = (lex_tbl[posi].init_str.find(c) != string::npos);
			bool is_fall_of = (lex_tbl[posi].init_str.size() == 0);

			//std::cout << "init_str test = " << lex_tbl[posi].init_str << std::endl;

			if (is_in_init_str || is_fall_of) {
				if (lex_tbl[posi].action & ACCUM_ACT) {
					lex_state = lex_tbl[posi].new_state;

					text_buf += c;
					pgetch->advance();
					ok = true;

					break;
				}
				if (lex_tbl[posi].action & RETURN_ACT) {
					lex_state = lex_tbl[posi].new_state;

					std::string text_bufa = text_buf;
					text_buf = "";

					if (!is_fall_of) {
						text_bufa += c;
						pgetch->advance();
					}

					return (*lex_tbl[posi].tok_fun)(text_bufa);
				}
				if (lex_tbl[posi].action & CONSUME_ACT) {
					lex_state = lex_tbl[posi].new_state;

					pgetch->advance();

					ok = true;

					break;
				}

			} else {
				++posi;
			}

		} while ((posi < lex_tbl.size()) && (lex_tbl[posi].state == lex_state));

		if (! ok) {
			break;
		}
	}

	std::cout << "getsym: parse error: fallen out of main loop" << std::endl;

	return Tok(SYM_ERR, "err: fall out of main loop");

}


void lexer_test_fun(std::string line)
{
	std::cout << "line = " << line << std::endl;

	GetchClass* p = new GetchString(line);
	Lexer<Token> lex(ltbl);

	lex.init(*p, 0);

	while (1) {
		Token tok = lex.getsym();

		std::cout << "sym_code = " << tok.sym_code << " ident = " << tok.ident <<
				" vald = " << tok.vald << " vali = " << tok.vali << std::endl;

		if (tok.sym_code == SYM_EOF) {
			break;
		}
	}

	delete p;
}

void parser_test_fun(std::string line)
{
	std::cout << "line = " << line << std::endl;

	GetchClass* p = new GetchString(line);
	Lexer<Token> lex(ltbl);

	ParserExpr<PolyExpr> parse_poly(lex);

	lex.init(*p, 0);

	PolyExpr res;

	res.set_vars("xyzab");

	parse_poly.parse_init();
	parse_poly.parse_expression(res);

	std::cout << "result = " << res.poly << std::endl;

	delete p;
}


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












