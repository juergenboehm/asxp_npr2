
#ifndef __poly_h_
#define __poly_h_

#include <iostream>

#include <vector>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define CHECK_TERMS(x) check_vector_ordered(x)
//#define CHECK_TERMS(x)

using namespace std;


template< typename A >
struct Select {};

template< class C, int nvars >
class Poly {

public:

	class Term;

	Poly<C, nvars>();
	Poly<C, nvars>(const C a0, int i, int exp);
	Poly<C, nvars>(const Term & t0);
	virtual ~Poly<C,nvars>();

	int degree();

	int read(char* poly_str);

	Poly<C,nvars> & add(const Poly<C,nvars> & p);
	Poly<C, nvars> & term_mult(const Term & tmult, Poly<C, nvars> & q );
	Poly<C,nvars> & mul(const Poly<C,nvars> & p);
	Poly<C,nvars> & pow(int e);


	Poly<C,nvars> & diff(int i);

	C eval(const C* vals) const;
	Poly<C, nvars-1> & eval_last(const C & val, Poly<C, nvars-1> & perg) const;
	int get_z_coef_fast(C* sl, C* coefs) const;

	Poly<C,nvars> & subst(const Poly<C,nvars>* subst_lis);

	Poly<C,nvars> & subst_mat(int ndim, C* subst_mat);

	char* sprint(char* buf) const;


public:
	class Term {
	public:
		Term();
		Term(C a, int i, int expo);
		Term(const Term & tright) {
			// cout << "Copy Construct Term" << endl;
			coef = tright.coef;
			memcpy(expos, tright.expos, nvars*sizeof(int));
		}

		Term & operator = (const Term & tright) {
			//cout << "Assigning Term" << endl;
			coef = tright.coef;
			memcpy(expos,tright.expos, nvars*sizeof(int));
			return *this;
		}

		virtual ~Term();

		int degree() const;
		bool operator < (const Term & t2) const;
		void diff(int i);
		void subst(const Poly<C, nvars> *subst_lis, Poly<C, nvars> & perg);

		template< typename C1 >
		int sprintcoefImpl(Select< C1 >, char* buf, bool with_sign) const;

		int sprintcoefImpl(Select< double >, char* buf, bool with_sign) const;

		int sprintcoef(char* buf, bool with_sign) const;

		char* sprint(char* buf) const;

		template< typename C1 >
		C powxImpl(Select<C1>, C val, int x) const;

		double powxImpl(Select<double>, double val, int x) const;

		C powx(C val, int x) const;

		bool same_monom(const Term & t1);

		void mult(const Term & t1);

		C eval(const C* vals) const;

		int sum_z_coef_fast(C* sl, C* coefs) const;

		typename Poly<C, nvars-1>::Term eval_last(const C & val) const;

	public:
		C coef;
		int expos[nvars];

	};

	typedef vector< Term > termVec;
	typedef typename vector< Term >::iterator vecIt;
	typedef typename vector< Term >::const_iterator vecItconst;


	termVec terms;

};

template< class C, int nvars>
ostream & operator << (ostream & os, const Poly<C, nvars> & p );



template< class C, int nvars >
Poly<C,nvars>::Poly()
{

}

template< class C, int nvars >
Poly<C,nvars>::Poly(const C a, int i, int expo)
{
	//cerr << "Entering constructor" << endl;
	terms.push_back(Term(a, i, expo));
	//cerr << "terms.size = " << terms.size() << endl;
}

template<class C, int nvars>
Poly<C, nvars>::Poly(const Term & t0)
{
	terms.push_back(t0);
}


template< class C, int nvars >
Poly<C, nvars>::~Poly() {
}

template< class C, int nvars >
int Poly<C, nvars>::degree()
{
	int deg_max = -1;
	for(vecIt it = terms.begin(); it < terms.end(); ++it) {
		int deg_it = it->degree();
		deg_max = ::max(deg_max, deg_it);
	}
	return deg_max;
}


template< class C, int nvars>
int read_poly(char* poly_str, Poly<C, nvars> & perg);


template< class C, int nvars>
int Poly<C, nvars>::read(char* poly_str) {

	int ret = read_poly<C, nvars>(poly_str, *this);
	return ret;
}

template< typename T>
bool check_vector_ordered(vector<T> & termsx)
{
	for(int i = 0; i < int(termsx.size()) - 1; ++i) {
		assert(termsx[i] < termsx[i+1]);
	}
	return true;
}


template< class C, int nvars >
Poly<C,nvars> & Poly<C, nvars>::add(const Poly<C,nvars> & p)
{
	size_t i = 0;
	size_t j = 0;

	vector<Term> new_terms;

	while (i < terms.size() && j < p.terms.size()) {
		const Term & ti = terms[i];
		const Term & tj = p.terms[j];
		if (ti < tj) {
			if (ti.coef != 0) {
				new_terms.push_back(ti);
			}
			++i;
		} else if (tj < ti) {
			if (tj.coef != 0) {
				new_terms.push_back(tj);
			}
			j++;
		} else {
			Term tsum;
			tsum.coef = ti.coef + tj.coef;
			if (tsum.coef != 0) {
				memcpy(tsum.expos, ti.expos, sizeof(int) * nvars);
				new_terms.push_back(tsum);
			}
			++i;
			++j;
		}
	}
	while (i < terms.size()) {
		if (terms[i].coef != 0) {
			new_terms.push_back(terms[i]);;
		}
		++i;
	}
	while (j < p.terms.size()) {
		if (p.terms[j].coef != 0) {
			new_terms.push_back(p.terms[j]);
		}
		++j;
	}

	CHECK_TERMS(new_terms);

	terms = new_terms;

	return *this;

}

template< class C, int nvars>
Poly<C, nvars> & Poly<C,nvars>::term_mult(const Term & tmult, Poly<C, nvars> & q )
{
	q.terms = terms;
	size_t i = 0;
	while (i < q.terms.size()) {
		q.terms[i].mult(tmult);
		++i;
	}

	CHECK_TERMS(q.terms);

	return *this;
}

template< class C, int nvars >
Poly<C,nvars> & Poly<C, nvars>::mul(const Poly<C,nvars> & p)
{
	size_t j = 0;
	Poly<C, nvars> p_sum_total;
	Poly<C, nvars> q;

	while (j < p.terms.size()) {
		term_mult(p.terms[j], q);
		p_sum_total.add(q);
		++j;
	}

	CHECK_TERMS(p_sum_total.terms);

	*this = p_sum_total;
	return *this;
}

template< class C, int nvars >
Poly<C,nvars> & Poly<C, nvars>::pow(int e)
{
	if (e == 0) {
		Poly<C, nvars> px(1, 0, 0);
		*this = px;
		return *this;
	}

	Poly<C, nvars> px = *this;
	for(int i = 0; i < e - 1; ++i) {
		mul(px);
	}
	return *this;
}


template< class C, int nvars >
Poly<C,nvars> & Poly<C, nvars>::diff(int i)
{

	vecIt j = terms.begin();

	while (j < terms.end()) {
		Term & ti = *j;
		if (ti.expos[i] > 0) {
			ti.coef = ti.expos[i] * ti.coef;
			ti.expos[i]--;
			++j;
		} else {
			terms.erase(j);
		}
	}

	return *this;
}

template< class C, int nvars >
C Poly<C, nvars>::eval(const C* vals) const
{
	C sum = 0;
	for(size_t i = 0; i < terms.size(); ++i) {
		sum += terms[i].eval(vals);
	}
	return sum;
}

template< class C, int nvars >
Poly<C, nvars-1> & Poly<C, nvars>::eval_last(const C & val, Poly<C, nvars-1> & perg) const
{
	perg = Poly<C,nvars-1>(0,0,0);
	for(vecItconst it = terms.begin(); it < terms.end(); ++it) {
		typename Poly<C, nvars-1>::Term t1 = it->eval_last(val);
		perg.add(Poly<C,nvars-1>(t1));
	}
	return perg;
}

template< class C, int nvars >
int Poly<C,nvars>::get_z_coef_fast(C* sl, C* coefs) const
{
	for(vecItconst it = terms.begin(); it < terms.end(); ++it) {
		it->sum_z_coef_fast(sl, coefs);
	}
	return 0;
}



template< class C, int nvars >
Poly<C,nvars> & Poly<C, nvars>::subst(const Poly<C,nvars>* subst_list)
{
	vecIt it = terms.begin();
	Poly<C, nvars> psum;
	while (it < terms.end()) {
		Poly<C, nvars> p_it;
		it->subst(subst_list, p_it);
		psum.add(p_it);
		//cerr << "psum = " << psum << endl;
		++it;
	}
	*this = psum;
	return *this;
}

template< class C, int nvars >
Poly<C,nvars> & Poly<C, nvars>::subst_mat(int ndim, C* subst_mat)
{
	Poly<C, nvars> psubsl[nvars];

	for(int i = 0; i < ndim; ++i) {
		Poly<C, nvars> psum;

		for(int j = 0; j < ndim; ++j) {
			int ind_ij = i * ndim + j;

			Poly<C,nvars> pacterm(subst_mat[ind_ij], j, 1);
			psum.add(pacterm);

		}
		psubsl[i] = psum;
	}
	for(int i = ndim; i < nvars; ++i) {
		psubsl[i] = Poly<C,nvars>(1,i,1);
	}
	return subst(psubsl);
}


template< class C, int nvars >
char* Poly<C, nvars>::sprint(char* buf) const
{
	char* p = buf;
	int tlen = terms.size();
	int nc = 0;
	nc = sprintf(p, "(%d):", tlen);
	p += nc;
	if (tlen > 0) {
		if (terms[0].coef < 0) {
			nc = sprintf(p, "-");
			p += nc;
		}
	}
	for(int i = 0; i < tlen - 1; ++i) {
		p = terms[i].sprint(p);
		char op = terms[i+1].coef < 0 ? '-' : '+';
		nc = sprintf(p, " %c ", op);
		p += nc;
	}
	if (tlen > 0) {
		p = terms[tlen-1].sprint(p);
	} else {
		Term a(0,0,0);
		nc = a.sprintcoef(p, 0);
		p += nc;
	}
	return p;
}

template< class C, int nvars>
ostream & operator << (ostream & os, const Poly<C, nvars> & p ) {
	char buf[2048];
	p.sprint(buf);
	os << buf;
	return os;
}


template< class C, int nvars >
Poly<C,nvars>::Term::Term() {
	coef = 0;
	for(int i = 0; i < nvars; ++i) {
		expos[i] = 0;
	}
}

template< class C, int nvars >
Poly<C,nvars>::Term::Term(C a0, int i, int expo) {
	coef = a0;
	for(int j = 0; j < nvars; ++j) {
		expos[j] = 0;
	}
	expos[i] = expo;
	//cerr << "(" << i << "," << expo << ")" << endl;
}


template< class C, int nvars >
Poly<C, nvars>::Term::~Term() {

}

template< class C, int nvars >
int Poly<C, nvars>::Term::degree () const
{

	int d = 0;
	for(int i = 0; i < nvars; ++i) {
		d += expos[i];
	}
	return d;
}

template< class C, int nvars >
bool Poly<C, nvars>::Term::operator < (const Poly<C, nvars>::Term & t) const {

	int d = degree();
	int d1 = t.degree();

	if (d < d1) {
		return true;
	}
	if (d > d1) {
		return false;
	}
	for(int i = 0; i < nvars; ++i) {
		if (expos[i] == t.expos[i]) {
			continue;
		} else {
			return expos[i] < t.expos[i];
		}
	}
	return false;
}


template< class C, int nvars >
void Poly<C, nvars>::Term::diff(int i) {
	coef *= expos[i];
	expos[i]--;
}

template< class C, int nvars >
void Poly<C, nvars>::Term::subst(const Poly<C, nvars> *subst_lis, Poly<C, nvars> & perg) {

	Poly<C, nvars>p(coef,0,0);
	for(int i = 0; i < nvars; ++i) {
		Poly<C,nvars> pvari = subst_lis[i];
		pvari.pow(expos[i]);
		p.mul(pvari);
	}
	perg = p;
}

template<class C, int nvars>
template< typename C1 >
int Poly<C, nvars>::Term::sprintcoefImpl(Select< C1 >, char* buf, bool with_sign) const
{
	return 0;
}

template<class C, int nvars>
int Poly<C,nvars>::Term::sprintcoefImpl(Select< double >, char* buf, bool with_sign) const
{
	double cf = with_sign ? coef : fabs(coef);
	return sprintf(buf, "%.2f", cf);
}

template<class C, int nvars>
int Poly<C,nvars>::Term::sprintcoef(char* buf, bool with_sign) const
{
	return sprintcoefImpl(Select<C>(), buf, with_sign);
}

template<class C, int nvars>
template< typename C1 >
C Poly<C,nvars>::Term::powxImpl(Select<C1>, C val, int x) const
{
	return val;
}

template<class C, int nvars>
double Poly<C,nvars>::Term::powxImpl(Select<double>, double val, int x) const
{
	return ::pow(val, x);
}

template<class C, int nvars>
C Poly<C,nvars>::Term::powx(C val, int x) const
{
	return powxImpl(Select<C>(), val, x);
}



template< class C, int nvars >
char* Poly<C, nvars>::Term::sprint(char* buf) const
{
	const char* vars = "xyzabcdefgh";
    char* pout = buf;
	int nchar = sprintcoef(pout, false);
	pout += nchar;
	for(int i = 0; i < nvars - 1; ++i) {
		if (expos[i] == 1) {
			nchar = sprintf(pout, " * %c", vars[i]);
			pout += nchar;
		} else if (expos[i] > 1) {
			nchar = sprintf(pout, " * %c^%d", vars[i], expos[i]);
			pout += nchar;
		}
	}
	if (expos[nvars-1] > 1) {
		nchar = sprintf(pout, "* %c^%d", vars[nvars-1], expos[nvars-1]);
		pout += nchar;
	} else if (expos[nvars-1] > 0) {
		nchar = sprintf(pout, "* %c", vars[nvars-1]);
		pout += nchar;
	}
	return pout;
}

template< class C, int nvars >
bool Poly<C, nvars>::Term::same_monom(const Term & t1) {
	return memcmp(expos, t1.expos, sizeof(int) * nvars ) == 0;
}


template< class C, int nvars >
void Poly<C, nvars>::Term::mult(const Term & t1) {
	for(int i = 0; i < nvars; ++i) {
		expos[i] += t1.expos[i];
	}
	coef *= t1.coef;
}

template< class C, int nvars >
C Poly<C, nvars>::Term::eval(const C* vals) const {
	C val = coef;
	for(int i = 0; i < nvars; ++i) {
		val *= powx(vals[i], expos[i]);
	}
	return val;
}

template< class C, int nvars >
typename Poly<C, nvars-1>::Term Poly<C, nvars>::Term::eval_last(const C & val) const {
	typename Poly<C,nvars-1>::Term t1;
	for(int i = 0; i < nvars - 1; ++i) {
		t1.expos[i] = expos[i];
	}
	t1.coef = coef * powx(val, expos[nvars-1]);
	return t1;
}

template< class C, int nvars >
int Poly<C,nvars>::Term::sum_z_coef_fast(C* sl, C* coefs) const
{
	C prod = coef;
	for(int i = 0; i < nvars - 1; ++i ) {
		prod *= powx(sl[i], expos[i]);
	}
	coefs[expos[nvars-1]] += prod;
	return 0;
}


typedef Poly<double, 1> Poly1;
typedef Poly<double, 2> Poly2;
typedef Poly<double, 3> Poly3;
typedef Poly<double, 4> Poly4;
typedef Poly<double, 5> Poly5;



#endif // __poly_h_
