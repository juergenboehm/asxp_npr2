
#ifndef __parser_h
#define __parser_h

#include <assert.h>

#include "poly.h"

//#define TRACE


#define SYM_NIL 0

#define SYM_PLUS 10
#define SYM_MINUS 11
#define SYM_TIMES 12
#define SYM_EXPO 13

#define SYM_LPAREN 20
#define SYM_RPAREN 21
#define SYM_SEMICOLON 30

#define SYM_NUMBER 40
#define SYM_VARIABLE 41

#define ERR_NUM_AFTER_EXPO 10
#define ERR_NO_RPAREN 11

int next_symbol(char* p_start, char * & p_sym_beg, char* & p_sym_next, int & symbol,
		long & numval, int & varindex );

template<class C>
int parse_expression(char* p, char* & p_next, Poly<C> & pvalue);

template<class C>
int parse_term(char* p, char* & p_next, Poly<C> & pvalue);

template<class C>
int parse_factor(char* p, char* & p_next, Poly<C> &  pvalue);

template<class C>
int parse_base(char* p, char* & p_next, Poly<C> & pvalue);


// grammar is
//
// expression = [+/-] term ( +|- term)*
// term = factor ( * factor)*
// factor = base [^number]
// base = variable | (expression) | number
//

template< class C>
int read_poly(char* poly_str, Poly<C> & perg)
{
	char *p_next;
	parse_expression(poly_str, p_next, perg);

	return 0;
}

template<class C>
int parse_expression(char* p, char* & p_next, Poly<C>  & pvalue)
{
	int nvars = pvalue.nvars;

	Poly<C> ptemp(nvars);
	char* pnexttemp;
	char* pbeg;
	int symbol;
	long numval;
	int varindex;
	bool first_minus = false;

#ifdef TRACE
	cout << "Enter expression:" << p << endl;
#endif

	next_symbol(p, pbeg, pnexttemp, symbol, numval, varindex);

	if (symbol == SYM_PLUS) {
		p = pnexttemp;
	} else if (symbol == SYM_MINUS) {
		p = pnexttemp;
		first_minus = true;
	};

	parse_term(p, p_next, ptemp);

	if (first_minus) {
		ptemp.mul(Poly<C>(-1,0,0, nvars));
	}

#ifdef TRACE
	cout << "term start" << ptemp << endl;
#endif

	p = p_next;

	do {

		next_symbol(p, pbeg, pnexttemp, symbol, numval, varindex);

		if (symbol == SYM_PLUS || symbol == SYM_MINUS) {
			p = pnexttemp;
			Poly<C> ptemp1(nvars);
			parse_term(p, pnexttemp, ptemp1);
			if (symbol == SYM_MINUS) {
				ptemp1.mul(Poly<C>(-1,0,0, nvars));
			}
#ifdef TRACE
			cout << "term follow = " << ptemp1 << endl;
#endif
			ptemp.add(ptemp1);
			p = pnexttemp;
		} else {
			p_next = p;
			pvalue = ptemp;
			return 0;
		}

	} while (1);

	return 0;
}


template<class C>
int parse_term(char* p, char* & p_next, Poly<C> & pvalue)
{
	int nvars = pvalue.nvars;

	Poly<C> ptemp(nvars);
	char* pnexttemp;
	char* pbeg;
	int symbol;
	long numval;
	int varindex;

#ifdef TRACE
	cout << "Enter term:" << p << endl;
#endif

	parse_factor(p, p_next, ptemp);

#ifdef TRACE
	cout << "factor start:" << ptemp << endl;
#endif

	p = p_next;

	do {

		next_symbol(p, pbeg, pnexttemp, symbol, numval, varindex);

		if (symbol == SYM_TIMES) {
			p = pnexttemp;
			Poly<C> ptemp1(nvars);
			parse_factor(p, pnexttemp, ptemp1);

#ifdef TRACE
			cout << "factor follow: " << ptemp1 << endl;
#endif

			ptemp.mul(ptemp1);
			p = pnexttemp;


		} else {
			p_next = p;
			pvalue = ptemp;
			return 0;
		}

	} while (1);

	return 0;
}

template<class C>
int parse_factor(char* p, char* &p_next, Poly<C> & pvalue)
{
	int nvars = pvalue.nvars;

	Poly<C> ptemp(nvars);
	char* pnexttemp;
	char* pbeg;
	int symbol;
	long numval;
	int varindex;

#ifdef TRACE
	cout << "Enter factor:" << p << endl;
#endif

	parse_base(p, pnexttemp, ptemp);

#ifdef TRACE
	cout << "base start:" << ptemp << endl;
#endif

	p = pnexttemp;

	next_symbol(p, pbeg, pnexttemp, symbol, numval, varindex);

	if (symbol == SYM_EXPO) {
		p = pnexttemp;
		next_symbol(p, pbeg, pnexttemp, symbol, numval, varindex);
		if (symbol == SYM_NUMBER) {
			ptemp.pow(numval);

#ifdef TRACE
			cout << "varindex = " << varindex << " expo = " << numval << endl;
			cout << "++++++ base full = " << ptemp << endl;
#endif
			p_next = pnexttemp;
			pvalue = ptemp;
			return 0;
		} else {
			assert(0);
			return ERR_NUM_AFTER_EXPO;
		}
	} else {
		p_next = p;
		pvalue = ptemp;
		return 0;
	}
}



template<class C>
int parse_base(char* p, char* & p_next, Poly<C> & pvalue)
{
	char* pbeg;
	char* p_nexttemp;
	int symbol;
	long numval;
	int varindex;

	int nvars = pvalue.nvars;

#ifdef TRACE
	cout << "Enter base:" << p << endl;
#endif

	next_symbol(p, pbeg, p_nexttemp, symbol, numval, varindex);

	if (symbol == SYM_NUMBER) {
		Poly<C> pval(numval, 0, 0, nvars);
		pvalue = pval;
		p_next = p_nexttemp;
		return 0;
	} else if (symbol == SYM_LPAREN) {
		p = p_nexttemp;
		parse_expression(p, p_nexttemp, pvalue);
		p = p_nexttemp;
		next_symbol(p, pbeg, p_nexttemp, symbol, numval, varindex);
		if (symbol == SYM_RPAREN) {
			p_next = p_nexttemp;
			return 0;
		} else {
			assert(0);
			return ERR_NO_RPAREN;
		}
	} else if (symbol == SYM_VARIABLE) {
		Poly<C> pval(1, varindex, 1, nvars);
		p_next = p_nexttemp;
		pvalue = pval;
		return 0;
	} else {
		assert(0);
	}
	return 0;
}



#endif
