
#ifndef __defined_libdai_factorsp_h
#define __defined_libdai_factorsp_h

#include <iostream>
#include <functional>
#include <cmath>
#include <dai/prob.h>
#include <dai/varset.h>
#include <dai/index.h>
#include <dai/util.h>
#include <dai/spvector.h>

namespace dai {

template <typename T, typename spvector_type>
class TFactor {
private:
	VarSet _vs;
	TProb<T,spvector_type> _p;
public:


	TFactor ( T p = 1) : _vs(), _p(1,p) {}

	TFactor( const Var &v ) : _vs(v), _p(v.states()) {}

	TFactor( const VarSet& vars ) : _vs(vars), _p() {
		_p = TProb<T,spvector_type>( BigInt_size_t( _vs.nrStates() ) );
	}

	TFactor( const VarSet& vars, T p ) : _vs(vars), _p() {
		_p = TProb<T,spvector_type>( BigInt_size_t( _vs.nrStates() ), p );
	}


	template<typename S>
	TFactor( const VarSet& vars, const std::vector<S> &x ) : _vs(vars), _p() {
		DAI_ASSERT( x.size() == vars.nrStates() );
		_p = TProb<T,spvector_type>( x.begin(), x.end(), x.size() );
	}


	TFactor( const VarSet& vars, const T* p ) : _vs(vars), _p() {
		size_t N = BigInt_size_t( _vs.nrStates() );
		_p = TProb<T,spvector_type>( p, p + N, N );
	}

	TFactor( const VarSet& vars, const TProb<T,spvector_type> &p ) : _vs(vars), _p(p) {
		DAI_ASSERT( _vs.nrStates() == _p.size() );
	}

	TFactor( const std::vector<Var> &vars, const std::vector<T> &p ) : _vs(vars.begin(), vars.end(), vars.size()), _p(p.size()) {
		BigInt nrStates = 1;
		for( size_t i = 0; i < vars.size(); i++ )
			nrStates *= vars[i].states();
		DAI_ASSERT( nrStates == p.size() );
		Permute permindex(vars);
		for( size_t li = 0; li < p.size(); ++li )
			_p.set( permindex.convertLinearIndex(li), p[li] );
	}



	void set( size_t i, T val ) { _p.set( i, val ); }

	T get( size_t i ) const { return _p[i]; }



	const TProb<T,spvector_type>& p() const { return _p; }

	TProb<T,spvector_type>& p() { return _p; }

	T operator[] (size_t i) const { return _p[i]; }

	const VarSet& vars() const { return _vs; }

	VarSet& vars() { return _vs; }


	size_t nrStates() const { return _p.size(); }

	T entropy() const { return _p.entropy(); }

	T max() const { return _p.max(); }

	T min() const { return _p.min(); }

	T sum() const { return _p.sum(); }

	T sumAbs() const { return _p.sumAbs(); }

	T maxAbs() const { return _p.maxAbs(); }

	bool hasNaNs() const { return _p.hasNaNs(); }

	bool hasNegatives() const { return _p.hasNegatives(); }

	T strength( const Var &i, const Var &j ) const;

	bool operator==( const TFactor<T,spvector_type>& y ) const {
		return (_vs == y._vs) && (_p == y._p);
	}



	TFactor<T,spvector_type> operator- () const {
		// Note: the alternative (shorter) way of implementing this,
				// return TFactor<T,spvector_type>( _vs, _p.abs() );
		// is slower because it invokes the copy constructor of TProb<T>
		TFactor<T,spvector_type> x;
		x._vs = _vs;
		x._p = -_p;
		return x;
	}

	TFactor<T,spvector_type> abs() const {
		TFactor<T,spvector_type> x;
		x._vs = _vs;
		x._p = _p.abs();
		return x;
	}

	TFactor<T,spvector_type> exp() const {
		TFactor<T,spvector_type> x;
		x._vs = _vs;
		x._p = _p.exp();
		return x;
	}


	TFactor<T,spvector_type> log(bool zero=false) const {
		TFactor<T,spvector_type> x;
		x._vs = _vs;
		x._p = _p.log(zero);
		return x;
	}


	TFactor<T,spvector_type> inverse(bool zero=true) const {
		TFactor<T,spvector_type> x;
		x._vs = _vs;
		x._p = _p.inverse(zero);
		return x;
	}


	TFactor<T,spvector_type> normalized( ProbNormType norm=NORMPROB ) const {
		TFactor<T,spvector_type> x;
		x._vs = _vs;
		x._p = _p.normalized( norm );
		return x;
	}



	TFactor<T,spvector_type>& randomize() { _p.randomize(); return *this; }

	TFactor<T,spvector_type>& setUniform() { _p.setUniform(); return *this; }

	TFactor<T,spvector_type>& takeAbs() { _p.takeAbs(); return *this; }

	TFactor<T,spvector_type>& takeExp() { _p.takeExp(); return *this; }


	TFactor<T,spvector_type>& takeLog( bool zero = false ) { _p.takeLog(zero); return *this; }


	T normalize( ProbNormType norm=NORMPROB ) { return _p.normalize( norm ); }



	TFactor<T,spvector_type>& fill (T x) { _p.fill( x ); return *this; }

	TFactor<T,spvector_type>& operator+= (T x) { _p += x; return *this; }

	TFactor<T,spvector_type>& operator-= (T x) { _p -= x; return *this; }

	TFactor<T,spvector_type>& operator*= (T x) { _p *= x; return *this; }

	TFactor<T,spvector_type>& operator/= (T x) { _p /= x; return *this; }

	TFactor<T,spvector_type>& operator^= (T x) { _p ^= x; return *this; }



	TFactor<T,spvector_type> operator+ (T x) const {
		// Note: the alternative (shorter) way of implementing this,
				// TFactor<T,spvector_type> result(*this);
		// result._p += x;
		// is slower because it invokes the copy constructor of TFactor<T,spvector_type>
		TFactor<T,spvector_type> result;
		result._vs = _vs;
		result._p = p() + x;
		return result;
	}

	TFactor<T,spvector_type> operator- (T x) const {
		TFactor<T,spvector_type> result;
		result._vs = _vs;
		result._p = p() - x;
		return result;
	}

	TFactor<T,spvector_type> operator* (T x) const {
		TFactor<T,spvector_type> result;
		result._vs = _vs;
		result._p = p() * x;
		return result;
	}

	TFactor<T,spvector_type> operator/ (T x) const {
		TFactor<T,spvector_type> result;
		result._vs = _vs;
		result._p = p() / x;
		return result;
	}

	TFactor<T,spvector_type> operator^ (T x) const {
		TFactor<T,spvector_type> result;
		result._vs = _vs;
		result._p = p() ^ x;
		return result;
	}




	template<typename binOp> TFactor<T,spvector_type>& binaryOp( const TFactor<T,spvector_type> &g, binOp op ) {
		if( _vs == g._vs ){ // optimize special case

//			std::cout << "binaryOp: Same variable sets\n";
			_p.pwBinaryOp( g._p, op );
		}
		else {

//			std::cout << "binaryOp: Different variable sets\n";

			*this = pointwiseOp( *this, g, op );
		}
		return *this;
 	}


	TFactor<T,spvector_type>& operator+= (const TFactor<T,spvector_type>& g) {
//		std::cout<< this->vars() << " += " << g.vars() << "\n";

		if( _vs == g._vs ){ // optimize special case
//			std::cout << "Same variable sets\n";

			_p.pwBinaryOpAddSub_eff( g._p, std::plus<T>() );
		}
		else{
//			std::cout << "Different variable sets\n";

			return binaryOp( g, std::plus<T>() );
		}

		return *this;
	}


	TFactor<T,spvector_type>& operator-= (const TFactor<T,spvector_type>& g) {
//		std::cout<< this->vars() << " -= " << g.vars() << "\n";

		if( _vs == g._vs ){ // optimize special case
//			std::cout << "Same variable sets\n";

			_p.pwBinaryOpAddSub_eff( g._p, std::minus<T>() );
		}
		else{
//			std::cout << "Different variable sets\n";

			return binaryOp( g, std::minus<T>() );
		}

		return *this;
	}


	TFactor<T,spvector_type>& operator*= (const TFactor<T,spvector_type>& g) {
		// Note that the following implementation is slow, because it doesn't exploit the special case of value 0
		// return binaryOp( g, std::multiplies<T>() );

//		std::cout<< this->vars() << " *= " << g.vars() << "\n";

		if( _vs == g._vs ){ // optimize special case
//			std::cout << "Same variable sets\n";
			_p.pwBinaryOpMultDiv_eff( g._p, std::multiplies<T>() );
		}
		else{

//			std::cout << "Different variable sets\n";

			*this = pointwiseOp( *this, g, std::multiplies<T>(), p().def() == (T)0 && g.p().def() == (T)0);
		}
		return *this;
	}


	TFactor<T,spvector_type>& operator/= (const TFactor<T,spvector_type>& g) {
//		std::cout<< this->vars() << " /= " << g.vars() << "\n";

		if( _vs == g._vs ){ // optimize special case
//			std::cout << "Same variable sets\n";
			_p.pwBinaryOpMultDiv_eff( g._p, fo_divides0<T>() );
		}
		else{

//	    	std::cout << "Different variable sets\n";
			return binaryOp( g, fo_divides0<T>() );
		}
		return *this;
	}

	template<typename binOp> TFactor<T,spvector_type> binaryTr( const TFactor<T,spvector_type> &g, binOp op ) const {
		// OPTIMIZE ME
		return pointwiseOp( *this, g, op );

	}


	TFactor<T,spvector_type> operator+ (const TFactor<T,spvector_type>& g) const {
//		std::cout<< this->vars() << " + " << g.vars() << "\n";

		if( _vs == g._vs ){ // optimize special case

//			std::cout << "Same variable sets\n";

			TFactor<T,spvector_type> result( vars() );
			result.p() = p().pwBinaryTrAddSub_eff( g.p(), std::plus<T>() );

			return result;
		}
		else{

//			std::cout << "Different variable sets\n";
			return binaryTr(g,std::plus<T>());
		}
	}


	TFactor<T,spvector_type> operator- (const TFactor<T,spvector_type>& g) const {
//		std::cout<< this->vars() << " - " << g.vars() << "\n";

		if( _vs == g._vs ){ // optimize special case

//			std::cout << "Same variable sets\n";

			TFactor<T,spvector_type> result( vars() );
			result.p() = p().pwBinaryTrAddSub_eff( g.p(), std::minus<T>() );

			return result;
		}
		else{

//			std::cout << "Different variable sets\n";
			return binaryTr(g,std::minus<T>());
		}
	}


	TFactor<T,spvector_type> operator* (const TFactor<T,spvector_type>& g) const {
//		std::cout<< this->vars() << " * " << g.vars() << "\n";

		if( _vs == g._vs ){ // optimize special case

//			std::cout << "Same variable sets\n";

			TFactor<T,spvector_type> result( vars() );
			result.p() = p().pwBinaryTrMultDiv_eff( g.p(), std::multiplies<T>() );

			return result;
		}
		else{

//			std::cout << "Different variable sets\n";
			return pointwiseOp( *this, g, std::multiplies<T>(), p().def() == (T)0 && g.p().def() == (T)0);
		}
	}


	TFactor<T,spvector_type> operator/ (const TFactor<T,spvector_type>& g) const {
//		std::cout<< this->vars() << " / " << g.vars() << "\n";



		if( _vs == g._vs ){ // optimize special case

//			std::cout << "Same variable sets\n";

			TFactor<T,spvector_type> result( vars() );
			result.p() = p().pwBinaryTrMultDiv_eff( g.p(), fo_divides0<T>() );

			return result;
		}
		else{

//			std::cout << "Different variable sets\n";
			return binaryTr(g,fo_divides0<T>());
		}


	}




	TFactor<T,spvector_type> slice( const VarSet& vars, size_t varsState ) const;


	TFactor<T,spvector_type> embed(const VarSet & vars) const {
		DAI_ASSERT( vars >> _vs );
		if( _vs == vars )
			return *this;
		else
			return (*this) * TFactor<T,spvector_type>(vars / _vs, (T)1);
	}

	TFactor<T,spvector_type> marginal(const VarSet &vars, bool normed=true) const;

	TFactor<T,spvector_type> maxMarginal(const VarSet &vars, bool normed=true) const;
};


template<typename T, typename spvector_type> TFactor<T,spvector_type> TFactor<T,spvector_type>::slice( const VarSet& vars, size_t varsState ) const {
	DAI_ASSERT( vars << _vs );
	VarSet varsrem = _vs / vars;

	TFactor<T,spvector_type> result( varsrem, p().def() );
	for( typename TProb<T,spvector_type>::const_iterator it = p().begin(); it != p().end(); it++ ) {
		State s( _vs, it->first );
		size_t vars_s = BigInt_size_t( s( vars ) );
		if( vars_s == varsState )
			result.set( BigInt_size_t( s(varsrem) ), it->second );
	}

	/* SLOW BECAUSE IT ITERATES OVER ALL VALUES */
	// OPTIMIZE ME
	/* TFactor<T,spvector_type> res( varsrem, T(0) );
  IndexFor i_vars (vars, _vs);
  IndexFor i_varsrem (varsrem, _vs);
  for( size_t i = 0; i < nrStates(); i++, ++i_vars, ++i_varsrem )
  if( (size_t)i_vars == varsState )
 res.set( i_varsrem, _p[i] );

  if( !((result.p() <= res.p()) && (res.p() <= result.p())) ) {
  std::cerr << result << std::endl;
  std::cerr << res << std::endl;
  DAI_ASSERT( ((result.p() <= res.p()) && (res.p() <= result.p())) );
  }*/

	return result;
}


template<typename T, typename spvector_type> TFactor<T,spvector_type> TFactor<T,spvector_type>::marginal(const VarSet &vars, bool normed) const {
	VarSet res_vars = vars & _vs;

//	DAI_ASSERT( !isnan(p().def()) );
//	VarSet rem(_vs / res_vars);

	TFactor<T,spvector_type> result( res_vars, (T)0);

	//allocate memory
	result.p().p().reserve(std::min( p().p().nrNonDef(), res_vars.nrStates().get_ui()) );

	for( typename TProb<T,spvector_type>::const_iterator it = p().begin(); it != p().end(); it++ ) {
		State s(_vs, it->first);
		size_t res_vars_s = BigInt_size_t( s( res_vars) );
		result.set( res_vars_s, result[res_vars_s] + it->second );
	}
//	std::cout << "capacity: " << result.p().p().nonDef().capacity() << "\n";
//	std::cout << "size: " << result.p().p().nonDef().size() << "\n";

	if( normed )
		result.normalize( NORMPROB );

	return result;
}


template<typename T, typename spvector_type> TFactor<T,spvector_type> TFactor<T,spvector_type>::maxMarginal(const VarSet &vars, bool normed) const {
	VarSet res_vars = vars & _vs;

	VarSet rem(_vs / res_vars);
	TFactor<T,spvector_type> result( res_vars, p().def() );
	for( typename TProb<T,spvector_type>::const_iterator it = p().begin(); it != p().end(); it++ ) {
		State s( _vs, it->first );
		size_t res_vars_s = BigInt_size_t( s( res_vars ) );
		if( it->second > result[res_vars_s] )
			result.set( res_vars_s, it->second );
	}

	/* SLOW BECAUSE IT ITERATES OVER ALL VALUES */
	/*
  TFactor<T,spvector_type> res( res_vars, 0.);
 IndexFor i_res( res_vars, _vs );
  for( size_t i = 0; i < _p.size(); i++, ++i_res )
  if( _p[i] > res._p[i_res] )
  res.set( i_res, _p[i] );

  if( !((result.p() <= res.p()) && (res.p() <= result.p())) ) {
  std::cerr << result << std::endl;
  std::cerr << res << std::endl;
  DAI_ASSERT( ((result.p() <= res.p()) && (res.p() <= result.p())) );
  }
	 */
	if( normed )
		result.normalize( NORMPROB );

	return result;
}

template<typename T, typename spvector_type> T TFactor<T,spvector_type>::strength( const Var &i, const Var &j ) const {
	DAI_DEBASSERT( _vs.contains( i ) );
	DAI_DEBASSERT( _vs.contains( j ) );
	DAI_DEBASSERT( i != j );
	VarSet ij(i, j);

	T max = 0.0;
	for( size_t alpha1 = 0; alpha1 < i.states(); alpha1++ )
		for( size_t alpha2 = 0; alpha2 < i.states(); alpha2++ )
			if( alpha2 != alpha1 )
				for( size_t beta1 = 0; beta1 < j.states(); beta1++ )
					for( size_t beta2 = 0; beta2 < j.states(); beta2++ )
						if( beta2 != beta1 ) {
							size_t as = 1, bs = 1;
							if( i < j )
								bs = i.states();
							else
								as = j.states();
							T f1 = slice( ij, alpha1 * as + beta1 * bs ).p().divide( slice( ij, alpha2 * as + beta1 * bs ).p() ).max();
							T f2 = slice( ij, alpha2 * as + beta2 * bs ).p().divide( slice( ij, alpha1 * as + beta2 * bs ).p() ).max();
							T f = f1 * f2;
							if( f > max )
								max = f;
						}
}


template<typename T, typename spvector_type, typename binaryOp> TFactor<T,spvector_type> pointwiseOp( const TFactor<T,spvector_type> &f, const TFactor<T,spvector_type> &g, binaryOp op, bool fast=false ) {

//	std::cout << "f = " << f << "\n\n";
//	std::cout << "g = " << g << "\n\n";

	//	if( f.vars() == g.vars() ) { // optimizate special case
	//
	//		std::cout << "Should not be here: factor.h, pointwiseOp()\n";
	//		exit(1);
	//
	//		TFactor<T,spvector_type> result( f.vars() );
	//		result.p() = f.p().pwBinaryTr( g.p(), op );
	//		return result;
	//	} else {

//	std::cout << "pointwiseOp: Different variable set\n";

	// Union of variables
	VarSet un( f.vars() | g.vars() );
	// Intersection of variables
	VarSet is( f.vars() & g.vars() );
	// Result factor
	TFactor<T,spvector_type> result( un, op( f.p().def(), g.p().def() ) );

//	std::cout << "container size before: "<< result.p().p().nrNonDef()<<"\n";

	if( fast ) {

//		std::cout << "\nMult Div \n";
		//allocate memory
		result.p().p().reserve(std::max(f.p().p().nrNonDef(), g.p().p().nrNonDef()));

		// For all non-default states of f and all non-default states of g
		for( typename TProb<T,spvector_type>::const_iterator fit = f.p().begin(); fit != f.p().end(); fit++ ) {
			// calculate state of f
			State fs( f.vars(), fit->first );

			for( typename TProb<T,spvector_type>::const_iterator git = g.p().begin(); git != g.p().end(); git++ ) {
				// calculate state of g
				State gs( g.vars(), git->first );

				// check whether these states are compatible
				bool compatible = true;

				for( typename VarSet::const_iterator v = is.begin(); v != is.end() && compatible; v++ ){

					if( fs(*v) != gs(*v) )
						compatible = false;
				}

				if( compatible ) {
					State fgs = fs;
					fgs.insert( gs.begin(), gs.end() );
					result.set( BigInt_size_t(fgs(un)), op( fit->second, git->second ) );
				}

//				std::cout << "container size: "<< result.p().p().nonDef().capacity()<<"\n";
			}
		}
	} else {

//		std::cout << "\nAdd Sub \n";
		//allocate memory
		result.p().p().reserve( f.p().p().nrNonDef() + g.p().p().nrNonDef() );

		// For all non-default states of f and all states of g
		for( typename TProb<T,spvector_type>::const_iterator fit = f.p().begin(); fit != f.p().end(); fit++ ) {
			State fs( f.vars(), fit->first );

			for( State g_minus_f_s(g.vars() / f.vars()); g_minus_f_s.valid(); g_minus_f_s++ ) {
				State fgs = g_minus_f_s.get();
				fgs.insert( fs.begin(), fs.end() );
				result.set( BigInt_size_t(fgs(un)), op( fit->second, g[BigInt_size_t(fgs(g.vars()))] ) );
//				std::cout << "container size: "<< result.p().p().nonDef().capacity()<<"\n";
			}

		}

//		std::cout << "\n result1 " << result << "\n";

		// For all states of f and all non-default states of g
		for( typename TProb<T,spvector_type>::const_iterator git = g.p().begin(); git != g.p().end(); git++ ) {
			State gs( g.vars(), git->first );

			for( State f_minus_g_s(f.vars() / g.vars()); f_minus_g_s.valid(); f_minus_g_s++ ) {
				State fgs = f_minus_g_s.get();
				fgs.insert( gs.begin(), gs.end() );
				result.set( BigInt_size_t(fgs(un)), op( f[BigInt_size_t(fgs(f.vars()))], git->second ) );
//				std::cout << "container size: "<< result.p().p().nonDef().capacity()<<"\n";
			}
		}
//		std::cout << "\n result2 " << result << "\n";
	}

//	std::cout << "container size after: "<< result.p().p().nrNonDef()<<"\n";
//
//	std::cout << "\n result " << result << "\n";
//	std::cout << "=========================================================\n\n";
	return result;

}



template<typename T, typename spvector_type> std::ostream& operator<< (std::ostream& os, const TFactor<T,spvector_type>& f) {
	// os << "(" << f.vars() << ", " << f.p() << ")";
	os << "[" << f.vars() << ", (";
	os << f.p();
	//for( size_t i = 0; i < f.nrStates(); i++ )
	//	os << (i == 0 ? "" : ", ") << f[i];
	os << ")]";
	return os;
}



template<typename T, typename spvector_type> T dist( const TFactor<T,spvector_type> &f, const TFactor<T,spvector_type> &g, ProbDistType dt ) {
	if( f.vars().empty() || g.vars().empty() )
		return -1;
	else {
		DAI_DEBASSERT( f.vars() == g.vars() );
		return dist( f.p(), g.p(), dt );
	}
}



template<typename T, typename spvector_type> TFactor<T,spvector_type> max( const TFactor<T,spvector_type> &f, const TFactor<T,spvector_type> &g ) {
	DAI_ASSERT( f.vars() == g.vars() );
	return TFactor<T,spvector_type>( f.vars(), max( f.p(), g.p() ) );
}



template<typename T, typename spvector_type> TFactor<T,spvector_type> min( const TFactor<T,spvector_type> &f, const TFactor<T,spvector_type> &g ) {
	DAI_ASSERT( f.vars() == g.vars() );
	return TFactor<T,spvector_type>( f.vars(), min( f.p(), g.p() ) );
}



template<typename T, typename spvector_type> T MutualInfo(const TFactor<T,spvector_type> &f) {
	DAI_ASSERT( f.vars().size() ==  2);
	VarSet::const_iterator it = f.vars().begin();
	Var i = *it; it++; Var j = *it;
	TFactor<T,spvector_type> projection = f.marginal(i) * f.marginal(j);
	return dist( f.normalized(), projection, DISTKL );
}


/// Represents a factor with values of type dai::Real.
typedef TFactor<Real, spvector<Real> > Factor;


Factor createFactorIsing( const Var &x, Real h );

Factor createFactorIsing( const Var &x1, const Var &x2, Real J );

Factor createFactorExpGauss( const VarSet &vs, Real beta );

Factor createFactorPotts( const Var &x1, const Var &x2, Real J );

Factor createFactorDelta( const Var &v, size_t state );

Factor createFactorDelta( const VarSet& vs, size_t state );

} // end of namespace dai


#endif
