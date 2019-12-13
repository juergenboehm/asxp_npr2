
#include <map>

#include "expr.h"


#include "parserex.h"



int ParserEx::parse_expression(Ex & pvalue)
{

	Ex ptemp;

	//ptemp.pre_init(pvalue);

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
		ptemp = -ptemp;
	}


	do {

		if (akt_sym.sym_code == SYM_PLUS || akt_sym.sym_code == SYM_MINUS) {

			bool is_minus = (akt_sym.sym_code == SYM_MINUS);

			Ex ptemp1;

			//ptemp1.pre_init(pvalue);

			akt_sym = lex.getsym();

			int err = parse_term(ptemp1);

			if (err) {
				return err;
			}

			if (is_minus) {

				ptemp = ptemp - ptemp1;

			} else {

				ptemp = ptemp + ptemp1;
			}

		} else {

			pvalue = ptemp;

			//std::cout << "parse_expression: return = " << to_str_infix(pvalue) << std::endl;

			return 0;
		}

	} while (1);

	pvalue = ptemp;

	//std::cout << "parse_expression: return = " << to_str_infix(pvalue) << std::endl;

	return 0;
}


int ParserEx::parse_term(Ex & pvalue)
{

	Ex ptemp;

	int err = parse_factor(ptemp);

	if (err) {
		return err;
	}

	do {

		if (akt_sym.sym_code == SYM_TIMES) {

			akt_sym = lex.getsym();

			Ex ptemp1;

			int err = parse_factor(ptemp1);

			if (err) {
				return err;
			}

			ptemp = ptemp * ptemp1;

		} else if (akt_sym.sym_code == SYM_DIVIDE) {

			akt_sym = lex.getsym();

			Ex ptemp1;

			int err = parse_factor(ptemp1);

			if (err) {
				return err;
			}

			ptemp = ptemp / ptemp1;

		} else {
			pvalue = ptemp;
			return 0;
		}

	} while (1);

	pvalue = ptemp;

	return 0;
}


int ParserEx::parse_factor(Ex & pvalue)
{

	Ex ptemp;

	int err = parse_base(ptemp);

	if (err) {
		return err;
	}

	if (akt_sym.sym_code == SYM_EXPO) {

		akt_sym = lex.getsym();

		if (akt_sym.sym_code == SYM_INTEGER) {

			ptemp = pow(ptemp, akt_sym.vali);

			akt_sym = lex.getsym();


			pvalue = ptemp;
			return 0;

		} else {

			Ex ptemp1;
			int err = parse_base(ptemp1);

			if (err) {
				return err;
			}

			pvalue = pow(ptemp, ptemp1);

			return 0;
		}
	} else {

		pvalue = ptemp;

		return 0;

	}
}


int ParserEx::parse_base(Ex & pvalue)
{


	if (akt_sym.sym_code == SYM_INTEGER || akt_sym.sym_code == SYM_DOUBLE) {

		Ex pval;

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
	} else if (akt_sym.sym_code == SYM_LBRACKET) {

		akt_sym = lex.getsym();

		Ex p1;

		int err = parse_subslis(p1);

		if (err) {
			return err;
		}

		if (akt_sym.sym_code == SYM_RBRACKET) {
			akt_sym = lex.getsym();

			pvalue = Ex(Ex("subslis"), p1);

			return 0;
		} else {

			return ERR_NO_RBRACKET;

		}

	} else if (akt_sym.sym_code == SYM_VARIABLE) {

		Ex pval;

		pval.set_ident(akt_sym.ident);
		pvalue = pval;

		akt_sym = lex.getsym();

		if (akt_sym.sym_code == SYM_LPAREN) {

			akt_sym = lex.getsym();

			Ex pargl;

			int err = parse_arglist(pargl);

			if (err) {
				return err;
			}

			pvalue = Ex(pval, pargl);

			if (akt_sym.sym_code == SYM_RPAREN) {
				akt_sym = lex.getsym();

			} else {

				assert(0);
				return ERR_NO_RPAREN;
			}

		}

		return 0;

	} else {
		assert(0);
		return ERR_BASE_MALFORMED;
	}
	return 0;
}

int ParserEx::parse_arglist(Ex & pvalue) {

	Ex pfirst;

	int err = parse_expression(pfirst);

	if (err) {
		return err;
	}

	if (akt_sym.sym_code == SYM_COMMA) {
		akt_sym = lex.getsym();

		Ex prest;

		int err = parse_arglist(prest);

		if (err) {
			return err;
		}

		pvalue = Ex(pfirst, prest);

		return 0;
	} else {

		pvalue = Ex(pfirst, nilEx);

		return 0;
	}

	return 0;
}

int ParserEx::parse_funexpr(Ex & pvalue) {

	Ex pval;

	pval.set_ident(akt_sym.ident);
	pvalue = pval;

	akt_sym = lex.getsym();

	if (akt_sym.sym_code == SYM_LPAREN) {

		akt_sym = lex.getsym();

		Ex pargl;

		int err = parse_arglist(pargl);

		if (err) {
			return err;
		}

		pvalue = Ex(pval, pargl);

		if (akt_sym.sym_code == SYM_RPAREN) {
			akt_sym = lex.getsym();

		} else {

			assert(0);
			return ERR_NO_RPAREN;
		}

	} else {
		return ERR_NO_LPAREN;
	}

	return 0;

}


int ParserEx::parse_subslis(Ex & pvalue) {

	Ex p1;

	int err = parse_expression(p1);

	if (err) {
		return err;
	}

	if (! is_symbol(p1)) {
		err = ERR_SYMBOL_EXPECTED;
		return err;
	}

	if (akt_sym.sym_code != SYM_EQUAL) {
		err = ERR_EQUALSIGN_EXPECTED;
		return err;
	}

	akt_sym = lex.getsym();

	Ex p12;

	err = parse_expression(p12);

	if (err) {
		return err;
	}

	Ex prest;

	if (akt_sym.sym_code == SYM_COMMA) {

		akt_sym = lex.getsym();

		err = parse_subslis(prest);

		if (err) {
			return err;
		}

		pvalue = Ex(Ex(Ex("="),Ex(p1, Ex(p12, nilEx))), prest);

		return 0;
	} else {

		pvalue = Ex(Ex(Ex("="), Ex(p1, Ex(p12,nilEx))), nilEx);

		return 0;
	}

	return 0;
}

std::map<string, bool> leadfun_map = { {"print", true}, {"poly", true}};
std::map<string, bool> masterfun_map = {{"subs", true}, {"evalall", true}};

bool is_leadfun(std::string str) {
	return (leadfun_map.find(str) != leadfun_map.end());
}

bool is_masterfun(std::string str) {
	return masterfun_map.find(str) != masterfun_map.end();
}


int ParserEx::parse_statement_subslis(Ex pvar, Ex & pvalue) {

	Ex p1;

	int err = parse_subslis(p1);

	if (err) {
		return err;
	}

	pvalue = Ex(Ex("="), Ex(pvar, Ex(Ex(Ex("subslis"), p1), nilEx)));

	return 0;

}


int ParserEx::parse_statement(Ex & pvalue) {

	if (akt_sym.sym_code == SYM_VARIABLE) {

		// statement line begin
		// must be an identifier

		std::string str = akt_sym.ident;

		if (is_leadfun(str)) {

			// is it poly or print

			int err = parse_funexpr(pvalue);

			if (err) {
				return err;
			}

			return 0;

		}

		// must be an assignment
		// pvar is assigned variable

		Ex pvar(str);

		akt_sym = lex.getsym();

		if (akt_sym.sym_code != SYM_EQUAL) {
			return ERR_EQUALSIGN_EXPECTED;
		}

		akt_sym = lex.getsym();

		// is rhs a [..] list?

		if (akt_sym.sym_code == SYM_LBRACKET) {

			akt_sym = lex.getsym();

			int err = parse_statement_subslis(pvar, pvalue);

			if (err) {
				return err;
			}

			if (akt_sym.sym_code == SYM_RBRACKET) {
				akt_sym = lex.getsym();
				return 0;
			}

			return ERR_NO_RBRACKET;

		}

		// does right hand side start with an identifier?

		if (akt_sym.sym_code == SYM_VARIABLE) {

			// yes!

			std::string str = akt_sym.ident;

			if (is_masterfun(str)) {

				// identifier is a "masterfun", subs or evalall

				Ex p1;

				int err = parse_funexpr(p1);

				if (err) {
					return err;
				}

				pvalue = Ex(Ex("="), Ex(pvar, Ex(p1, nilEx)));

				return 0;

			} else {

				// identifier leads in an ordinary expression

				Ex p1;

				int err = parse_expression(p1);

				if (err) {
					return err;
				}

				pvalue = Ex(Ex("="), Ex(pvar, Ex(p1, nilEx)));

				return 0;
			}
		} else {

			// is right hand side an expression not starting with
			// an identifier

			Ex p1;

			int err = parse_expression(p1);

			if (err) {
				return err;
			}

			pvalue = Ex(Ex("="), Ex(pvar, Ex(p1, nilEx)));

			return 0;
		}

	} else {
		return ERR_SYMBOL_EXPECTED;
	}

}

int ParserEx::parse_statement_list(Ex & pvalue) {

	Ex p1 = nilEx;
	int err = parse_statement(p1);

	if (akt_sym.sym_code == SYM_SEMICOLON) {

		akt_sym = lex.getsym();
		Ex prest;

		err = parse_statement_list(prest);

		if (err) {
			return err;
		}

		pvalue = Ex(p1, prest);

		return 0;
	} else if (akt_sym.sym_code == SYM_PERIOD) {

		akt_sym = lex.getsym();

		pvalue = Ex(p1, nilEx);

		return 0;
	}

	return ERR_BAD_STATEMENT_DELIM;
}

int ParserEx::parse_module(Ex & pvalue) {

	Ex pstatlis;

	int err = parse_statement_list(pstatlis);

	if (err) {
		return err;
	}

	pvalue = Ex(Ex("statlis"), pstatlis);

	return 0;
}





