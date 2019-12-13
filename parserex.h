#ifndef __parserex_h_
#define __parserex_h_

#include "expr.h"
#include "parser.h"


class ParserEx {

public:

	ParserEx(Lexer<Token> & lexa): lex(lexa) {};

	int parse_init() { akt_sym = lex.getsym(); return 0; }

	int parse_expression(Ex & pvalue);
	int parse_term(Ex & value);
	int parse_factor(Ex & pvalue);
	int parse_base(Ex & value);

	int parse_arglist(Ex & pvalue);

	int parse_funexpr(Ex & pvalue);

	int parse_subslis(Ex & pvalue);

	int parse_statement_subslis(Ex pvar, Ex & pvalue);

	int parse_statement(Ex & pvalue);

	int parse_statement_list(Ex & pvalue);

	int parse_module(Ex & pvalue);

	Lexer<Token>& lex;

	Token akt_sym;

};





#endif /* __parserex_h_ */
