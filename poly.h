
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

template< class C>
class Poly {

public:

	class Term;

	Poly<C>();
	Poly<C>(int nvarsa);
	Poly<C>(const C a0, int i, int exp, int nvarsa);
	Poly<C>(const Term & t0);
	virtual ~Poly<C>();

	Poly<C>(const Poly<C> &) = default;

	Poly<C> & operator=(const Poly<C> & rhs);

	int degree();

	int read(char* poly_str);

	Poly<C> & add(const Poly<C> & p);
	Poly<C> & term_mult(const Term & tmult, Poly<C> & q );
	Poly<C> & mul(const Poly<C> & p);
	Poly<C> & pow(int e);


	Poly<C> & diff(int i);

	C eval(const C* vals) const;
	Poly<C> & eval_last(const C & val, Poly<C> & perg) const;
	int get_z_coef_fast(C* sl, C* coefs) const;

	Poly<C> & subst(const Poly<C>* subst_lis);

	Poly<C> & subst_mat(int ndim, C* subst_mat);

	char* sprint(char* buf) const;


public:
	class Term {
	public:
		Term();
		Term(int nvarsa);

		Term(C a, int i, int expo, int nvarsa);
		Term(const Term & tright) {
			// cout << "Copy Construct Term" << endl;
			coef = tright.coef;
			expos = tright.expos;
		}

		Term & operator = (const Term & tright) {
			//cout << "Assigning Term" << endl;
			coef = tright.coef;
			expos = tright.expos;
			return *this;
		}

		virtual ~Term();

		int degree() const;
		bool operator < (const Term & t2) const;
		void diff(int i);
		void subst(const Poly<C> *subst_lis, Poly<C> & perg);

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

		typename Poly<C>::Term eval_last(const C & val) const;

	public:
		C coef;
		//int expos[nvars];


		vector<int> expos;

	};

	typedef vector< Term > termVec;
	typedef typename vector< Term >::iterator vecIt;
	typedef typename vector< Term >::const_iterator vecItconst;


	termVec terms;

	int nvars;

};

template< class C>
ostream & operator << (ostream & os, const Poly<C> & p );


template< class C>
Poly<C>::Poly(): nvars(0), terms()
{

}

template< class C>
Poly<C>::Poly(int nvarsa) : nvars(nvarsa), terms()
{

}

template< class C>
Poly<C>::Poly(const C a, int i, int expo, int nvarsa): nvars(nvarsa)
{
	//cerr << "Entering constructor" << endl;
	terms.push_back(Term(a, i, expo, nvarsa));
	//cerr << "terms.size = " << terms.size() << endl;
}

template<class C>
Poly<C>::Poly(const Term & t0): nvars(t0.expos.size())
{
	terms.push_back(t0);
}


template< class C>
Poly<C>::~Poly() {
}

template< class C>
Poly<C> & Poly<C>::operator=(const Poly<C> & rhs)
{
	nvars = rhs.nvars;
	terms = rhs.terms;

	return *this;
}


template< class C>
int Poly<C>::degree()
{
	int deg_max = -1;
	for(vecIt it = terms.begin(); it < terms.end(); ++it) {
		int deg_it = it->degree();
		deg_max = ::max(deg_max, deg_it);
	}
	return deg_max;
}


template< class C>
int read_poly(char* poly_str, Poly<C> & perg);


template< class C>
int Poly<C>::read(char* poly_str) {

	int ret = read_poly<C>(poly_str, *this);
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


template< class C>
Poly<C> & Poly<C>::add(const Poly<C> & p)
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
				tsum.expos = ti.expos;
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

template< class C>
Poly<C> & Poly<C>::term_mult(const Term & tmult, Poly<C> & q )
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

template< class C>
Poly<C> & Poly<C>::mul(const Poly<C> & p)
{
	size_t j = 0;
	Poly<C> p_sum_total(nvars);
	Poly<C> q(nvars);

	while (j < p.terms.size()) {
		term_mult(p.terms[j], q);
		p_sum_total.add(q);
		++j;
	}

	CHECK_TERMS(p_sum_total.terms);

	*this = p_sum_total;
	return *this;
}

template< class C>
Poly<C> & Poly<C>::pow(int e)
{
	if (e == 0) {
		Poly<C> px(1, 0, 0, nvars);
		*this = px;
		return *this;
	}

	Poly<C> px = *this;
	for(int i = 0; i < e - 1; ++i) {
		mul(px);
	}
	return *this;
}


template< class C>
Poly<C> & Poly<C>::diff(int i)
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

template< class C>
C Poly<C>::eval(const C* vals) const
{
	C sum = 0;
	for(size_t i = 0; i < terms.size(); ++i) {
		sum += terms[i].eval(vals);
	}
	return sum;
}

template< class C>
Poly<C> & Poly<C>::eval_last(const C & val, Poly<C> & perg) const
{
	perg = Poly<C>(0,0,0, nvars-1);
	for(vecItconst it = terms.begin(); it < terms.end(); ++it) {
		typename Poly<C>::Term t1 = it->eval_last(val);
		perg.add(Poly<C>(t1));
	}
	return perg;
}

template< class C>
int Poly<C>::get_z_coef_fast(C* sl, C* coefs) const
{
	for(vecItconst it = terms.begin(); it < terms.end(); ++it) {
		it->sum_z_coef_fast(sl, coefs);
	}
	return 0;
}



template< class C>
Poly<C> & Poly<C>::subst(const Poly<C>* subst_list)
{
	vecIt it = terms.begin();
	Poly<C> psum(nvars);
	while (it < terms.end()) {
		Poly<C> p_it(nvars);
		it->subst(subst_list, p_it);
		psum.add(p_it);
		//cerr << "psum = " << psum << endl;
		++it;
	}
	*this = psum;
	return *this;
}


template< class C>
Poly<C> & Poly<C>::subst_mat(int ndim, C* subst_mat)
{
	Poly<C> psubsl[nvars];

	for(int i = 0; i < ndim; ++i) {
		Poly<C> psum(nvars);

		for(int j = 0; j < ndim; ++j) {
			int ind_ij = i * ndim + j;

			Poly<C> pacterm(subst_mat[ind_ij], j, 1, nvars);
			psum.add(pacterm);

		}
		psubsl[i] = psum;
	}
	for(int i = ndim; i < nvars; ++i) {
		psubsl[i] = Poly<C>(1,i,1, nvars);
	}
	return subst(psubsl);
}



template< class C>
char* Poly<C>::sprint(char* buf) const
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
		Term a(0,0,0, nvars);
		nc = a.sprintcoef(p, 0);
		p += nc;
	}
	return p;
}

template< class C>
ostream & operator << (ostream & os, const Poly<C> & p ) {
	char buf[2048];
	p.sprint(buf);
	os << buf;
	return os;
}


template< class C>
Poly<C>::Term::Term() {
	coef = 0;
//	for(int i = 0; i < nvars; ++i) {
//		expos[i] = 0;
//	}
}

template< class C>
Poly<C>::Term::Term(int nvarsa) {
	coef = 0;
	expos.resize(nvarsa);
	for(int i = 0; i < nvarsa; ++i) {
		expos[i] = 0;
	}
}


template< class C>
Poly<C>::Term::Term(C a0, int i, int expo, int nvarsa) {
	coef = a0;
	expos.resize(nvarsa);
	for(int j = 0; j < nvarsa; ++j) {
		expos[j] = 0;
	}
	expos[i] = expo;
	//cerr << "(" << i << "," << expo << ")" << endl;
}


template< class C>
Poly<C>::Term::~Term() {

}

template< class C>
int Poly<C>::Term::degree () const
{
	int nvars = expos.size();

	int d = 0;
	for(int i = 0; i < nvars; ++i) {
		d += expos[i];
	}
	return d;
}

template< class C>
bool Poly<C>::Term::operator < (const Poly<C>::Term & t) const {

	int nvars = expos.size();

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


template< class C>
void Poly<C>::Term::diff(int i) {
	coef *= expos[i];
	expos[i]--;
}

template< class C>
void Poly<C>::Term::subst(const Poly<C> *subst_lis, Poly<C> & perg) {

	int nvars = expos.size();

	Poly<C>p(coef,0,0, nvars);
	for(int i = 0; i < nvars; ++i) {
		Poly<C> pvari = subst_lis[i];
		pvari.pow(expos[i]);
		p.mul(pvari);
	}
	perg = p;
}

template<class C>
template< typename C1 >
int Poly<C>::Term::sprintcoefImpl(Select< C1 >, char* buf, bool with_sign) const
{
	return 0;
}

template<class C>
int Poly<C>::Term::sprintcoefImpl(Select< double >, char* buf, bool with_sign) const
{
	double cf = with_sign ? coef : fabs(coef);
	return sprintf(buf, "%.2f", cf);
}

template<class C>
int Poly<C>::Term::sprintcoef(char* buf, bool with_sign) const
{
	return sprintcoefImpl(Select<C>(), buf, with_sign);
}

template<class C>
template< typename C1 >
C Poly<C>::Term::powxImpl(Select<C1>, C val, int x) const
{
	return val;
}

template<class C>
double Poly<C>::Term::powxImpl(Select<double>, double val, int x) const
{
	return ::pow(val, x);
}

template<class C>
C Poly<C>::Term::powx(C val, int x) const
{
	return powxImpl(Select<C>(), val, x);
}



template< class C>
char* Poly<C>::Term::sprint(char* buf) const
{
	int nvars = expos.size();

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

template< class C>
bool Poly<C>::Term::same_monom(const Term & t1) {
	bool erg = true;
	for(size_t i = 0; i < expos.size(); ++i) {
		if (expos[i] != t1.expos[i]) {
			erg = false;
			break;
		}
	}
	return erg;
}


template< class C>
void Poly<C>::Term::mult(const Term & t1) {
	int nvars = expos.size();

	for(int i = 0; i < nvars; ++i) {
		expos[i] += t1.expos[i];
	}
	coef *= t1.coef;
}

template< class C>
C Poly<C>::Term::eval(const C* vals) const {
	int nvars = expos.size();

	C val = coef;
	for(int i = 0; i < nvars; ++i) {
		val *= powx(vals[i], expos[i]);
	}
	return val;
}

template< class C>
typename Poly<C>::Term Poly<C>::Term::eval_last(const C & val) const {

	int nvars = expos.size();

	typename Poly<C>::Term t1(nvars-1);
	for(int i = 0; i < nvars - 1; ++i) {
		t1.expos[i] = expos[i];
	}
	t1.coef = coef * powx(val, expos[nvars-1]);
	return t1;
}

template< class C>
int Poly<C>::Term::sum_z_coef_fast(C* sl, C* coefs) const
{
	int nvars = expos.size();

	C prod = coef;
	for(int i = 0; i < nvars - 1; ++i ) {
		prod *= powx(sl[i], expos[i]);
	}
	coefs[expos[nvars-1]] += prod;
	return 0;
}


typedef Poly<double> Poly1;
typedef Poly<double> Poly2;
typedef Poly<double> Poly3;
typedef Poly<double> Poly4;
typedef Poly<double> Poly5;



#endif // __poly_h_
