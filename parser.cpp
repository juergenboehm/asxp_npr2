
#include <iostream>

#include "stdlib.h"
#include "string.h"

#include "parser.h"


int next_symbol(char* p_start, char * & p_sym_beg, char* & p_sym_next, int & symbol,
		long & numval, int & varindex )
{
	char *p = p_start;
	char *q;
	long val;

	val = 0;
	numval = 0;
	varindex = 0;
	symbol = SYM_NIL;

	while (*p == ' ')
		++p;

#ifdef TRACE
	cerr << "next_symbol p = " << p << endl;

	cerr << "trying opcode" << endl;
#endif

	switch (*p) {

	case '+': symbol = SYM_PLUS;
				break;
	case '-': symbol = SYM_MINUS;
				break;
	case '*': symbol = SYM_TIMES;
				break;
	case '^': symbol = SYM_EXPO;
				break;
	case '(': symbol = SYM_LPAREN;
				break;
	case ')': symbol = SYM_RPAREN;
				break;
	case ';': symbol = SYM_SEMICOLON;
				break;
	default: symbol = SYM_NIL;
				break;
	}
	if (symbol != SYM_NIL) {
#ifdef TRACE
		cout << "opcode found: symbol = " << symbol << endl;
#endif
		p_sym_beg = p;
		++p;
		p_sym_next = p;
		return 0;
	}

	val = strtol(p, &q, 10);

	if (q != p) {
		symbol = SYM_NUMBER;
		p_sym_beg = p;
		p_sym_next = q;
		numval = val;
#ifdef TRACE
		cout << "number found: val = " << numval << endl;
#endif
		return 0;
	}

#ifdef TRACE
	cerr << "trying variable" << endl;
#endif

	char* var_str = "xyzabcdef";

	char* p_vpos = strchr(var_str, *p);
	if (p_vpos != NULL) {

#ifdef TRACE
		cerr << "variable_case p = " << p << endl;
#endif
		varindex = p_vpos - var_str;
#ifdef TRACE
		cerr << "varindex = " << varindex << endl;
#endif
		symbol = SYM_VARIABLE;

		p_sym_beg = p;
		p_sym_next = ++p;
		return 0;
	} else {
		return 3;
	}
#ifdef TRACE
	cerr << "all cases failed" << endl;
#endif
	return 1;

}

