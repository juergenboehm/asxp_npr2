
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
#define SYMCHAR_STR "+-*/^();,=[]."


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
	case ';':	sym = SYM_SEMICOLON;
				break;
	case ',': 	sym = SYM_COMMA;
				break;
	case '=':	sym = SYM_EQUAL;
				break;
	case '[':	sym = SYM_LBRACKET;
				break;
	case ']':	sym = SYM_RBRACKET;
				break;
	case '.':	sym = SYM_PERIOD;
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

	std::unique_ptr<GetchClass> p_getch (new GetchString(line));
	Lexer<Token> lex(ltbl);

	ParserExpr<PolyExpr> parse_poly(lex);

	lex.init(*p_getch, 0);

	PolyExpr res;

	res.set_vars("xyzab");

	parse_poly.parse_init();
	parse_poly.parse_expression(res);

	std::cout << "result = " << res.poly << std::endl;

}














