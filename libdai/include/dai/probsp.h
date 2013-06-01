
#ifndef __defined_libdai_probsp_h
#define __defined_libdai_probsp_h


#include <cmath>
#include <vector>
#include <ostream>
#include <algorithm>
#include <numeric>
#include <functional>
#include <dai/util.h>
#include <dai/exceptions.h>
#include <dai/fo.h>
#include <dai/spvector.h>


namespace dai {



template <typename T, typename spvector_type>
class TProb {
public:
	typedef spvector_type container_type;
	typedef TProb<T,spvector_type> this_type;

private:
	container_type _p;

public:


	TProb() : _p() {}

	explicit TProb( size_t n ) : _p( n, (T)1 / n ) {}

	explicit TProb( size_t n, T p ) : _p( n, p ) {}

	//prob vector of size n able to store cap elements
	explicit TProb( size_t n, size_t cap ) : _p( n, cap ) {}

	template <typename TIterator>
	TProb( TIterator begin, TIterator end, size_t sizeHint) : _p( begin, end, sizeHint) {}


	template <typename S>
	TProb( const std::vector<S> &v) : _p( v, v.size()) {}

	typedef typename container_type::const_iterator const_iterator;
	typedef typename container_type::iterator iterator;
	typedef typename container_type::const_reverse_iterator const_reverse_iterator;
	typedef typename container_type::reverse_iterator reverse_iterator;



	iterator begin() { return _p.begin(); }
	const_iterator begin() const { return _p.begin(); }

	iterator end() { return _p.end(); }
	const_iterator end() const { return _p.end(); }

	reverse_iterator rbegin() { return _p.rbegin(); }
	const_reverse_iterator rbegin() const { return _p.rbegin(); }

	reverse_iterator rend() { return _p.rend(); }
	const_reverse_iterator rend() const { return _p.rend(); }


	void resize( size_t sz ) {
		_p.resize( sz );
	}



	T get( size_t i ) const { return _p[i]; }

	void set( size_t i, T val ) { _p.set( i, val ); }



	const container_type& p() const { return _p; }

	container_type& p() { return _p; }

	T operator[]( size_t i ) const { return get(i); }

	size_t size() const { return _p.size(); }

	size_t nrDef() const { return _p.nrDef(); }

	size_t nrNonDef() const { return _p.nrNonDef(); }

	T def() const { return _p.def(); }

	void setDef( T def ) { _p.setDef( def ); }


	template<typename unOp> T accumulateSum( T init, unOp op ) const {
		T t = op(init);
		for( const_iterator it = begin(); it != end(); it++ )
			t += op(it->second);
		t += nrDef() * op(def());
		return t;
	}


	template<typename unOp> T accumulateMax( T init, unOp op, bool minimize ) const {
		T t = op(init);
		if( minimize ) {
			for( const_iterator it = begin(); it != end(); it++ )
				t = std::min( t, op(it->second) );
			if( nrDef() )
				t = std::min( t, op(def()) );
		} else {
			for( const_iterator it = begin(); it != end(); it++ )
				t = std::max( t, op(it->second) );
			if( nrDef() )
				t = std::max( t, op(def()) );
		}
		return t;
	}

	T entropy() const { return -accumulateSum( (T)0, fo_plog0p<T>() ); }

	T max() const { return accumulateMax( (T)(-INFINITY), fo_id<T>(), false ); }

	T min() const { return accumulateMax( (T)INFINITY, fo_id<T>(), true ); }

	T sum() const { return accumulateSum( (T)0, fo_id<T>() ); }

	T sumAbs() const { return accumulateSum( (T)0, fo_abs<T>() ); }

	T maxAbs() const { return accumulateMax( (T)0, fo_abs<T>(), false ); }

	bool hasNaNs() const {
		if( isnan( def() ) && nrDef() )
			return true;
		else {
			bool foundnan = false;
			for( const_iterator it = begin(); it != end(); it++ )
				if( isnan( it->second ) ) {
					foundnan = true;
					break;
				}
			return foundnan;
		}
	}

	bool hasNegatives() const {
		if( (def() < 0) && nrDef() )
			return true;
		else {
			bool foundnegative = false;
			for( const_iterator it = begin(); it != end(); it++ )
				if( it->second < 0) {
					foundnegative = true;
					break;
				}
			return foundnegative;
		}
	}

	std::pair<size_t,T> argmax() const {
		T max;
		size_t arg;
		DAI_ASSERT( size() );
		if( nrDef() == size() ) {
			max = def();
			arg = 0;
		} else if( nrDef() > 0) {
			max = begin()->second;
			arg = begin()->first;
			size_t i = 0;
			size_t argdef = 0;
			for( const_iterator it = begin(); it != end(); it++ ) {
				if( it->second > max ) {
					max = it->second;
					arg = it->first;
				}
				if( it->first != i )
					argdef = i;
				i = it->first + 1;
			}
			if( def() > max ) {
				max = def();
				arg = argdef;
			}
		} else {
			max = begin()->second;
			arg = begin()->first;
			for( const_iterator it = begin(); it != end(); it++ )
				if( it->second > max ) {
					max = it->second;
					arg = it->first;
				}
		}
		return std::make_pair( arg, max );
	}

	size_t draw() {
		Real x = rnd_uniform() * sum();
		T s = 0;
		for( size_t i = 0; i < size(); i++ ) {
			s += get(i);
			if( s > x )
				return i;
		}
		return( size() - 1);
	}


	bool operator<( const this_type& q ) const {
		DAI_DEBASSERT( size() == q.size() );
		for( size_t i = 0; i < size(); i++ ) {
			T a = get(i);
			T b = q.get(i);
			if( a > b )
				return false;
			if( a < b )
				return true;
		}
		return false;
	}

	bool operator==( const this_type& q ) const {
		if( size() != q.size() )
			return false;
		return _p == q._p;
	}



	template<typename unaryOp> this_type pwUnaryTr( unaryOp op ) const {
		this_type r;
		r.setDef( op( def() ) );
		r._p.resize( size() );
		for( const_iterator it = begin(); it != end(); it++ ) {
			T new_val = op( it->second );
			if( new_val != r.def() )
				r._p.push_back( it->first, new_val );
		}
		return r;
	}

	this_type operator- () const { return pwUnaryTr( std::negate<T>() ); }

	this_type abs() const { return pwUnaryTr( fo_abs<T>() ); }

	this_type exp() const { return pwUnaryTr( fo_exp<T>() ); }


	this_type log(bool zero=false) const {
		if( zero )
			return pwUnaryTr( fo_log0<T>() );
		else
			return pwUnaryTr( fo_log<T>() );
	}


	this_type inverse(bool zero=true) const {
		if( zero )
			return pwUnaryTr( fo_inv0<T>() );
		else
			return pwUnaryTr( fo_inv<T>() );
	}


	this_type normalized( ProbNormType norm = dai::NORMPROB ) const {
		T Z = 0;
		if( norm == dai::NORMPROB )
			Z = sum();
		else if( norm == dai::NORMLINF )
			Z = maxAbs();
		if( Z == (T)0) {
			DAI_THROW(NOT_NORMALIZABLE);
			return *this;
		} else
			return pwUnaryTr( std::bind2nd( std::divides<T>(), Z ) );
	}



	template<typename unaryOp> this_type& pwUnaryOp( unaryOp op ) {
		setDef( op( def() ) );
		for( iterator it = begin(); it != end(); ) {
			T new_val = op( it->second );
			if( new_val != def() ) {
				it->second = new_val;
				it++;
			} else
				it = _p.erase( it );
		}
		return *this;
	}

	this_type& randomize() {
		setDef(0);
		for( size_t i = 0; i < size(); i++ )
			set( i, (T)rnd_uniform() );
		return *this;
	}

	this_type& setUniform () {
		setDef(0);
		for( size_t i = 0; i < size(); i++ )
			set( i, (T)1/size() );
		return *this;
	}

	this_type& takeAbs() { return pwUnaryOp( fo_abs<T>() ); }

	this_type& takeExp() { return pwUnaryOp( fo_exp<T>() ); }


	this_type& takeLog(bool zero=false) {
		if( zero ) {
			return pwUnaryOp( fo_log0<T>() );
		} else
			return pwUnaryOp( fo_log<T>() );
	}


	T normalize( ProbNormType norm = dai::NORMPROB ) {
		T Z = 0;
		if( norm == dai::NORMPROB )
			Z = sum();
		else if( norm == dai::NORMLINF )
			Z = maxAbs();
		if( Z == (T)0)
			DAI_THROW(NOT_NORMALIZABLE);
		else
			*this /= Z;
		return Z;
	}



	this_type& fill( T x ) {
		for( size_t i = 0; i < size(); i++ )   set( i, x );
		return *this;
	}

	this_type& operator+= (T x) {
		if( x != 0)
			return pwUnaryOp( std::bind2nd( std::plus<T>(), x ) );
		else
			return *this;
	}

	this_type& operator-= (T x) {
		if( x != 0)
			return pwUnaryOp( std::bind2nd( std::minus<T>(), x ) );
		else
			return *this;
	}

	this_type& operator*= (T x) {
		if( x !=  1)
			return pwUnaryOp( std::bind2nd( std::multiplies<T>(), x ) );
		else
			return *this;
	}

	this_type& operator/= (T x) {
		if( x !=  1)
			return pwUnaryOp( std::bind2nd( fo_divides0<T>(), x ) );
		else
			return *this;
	}

	this_type& operator^= (T x) {
		if( x != (T)1)
			return pwUnaryOp( std::bind2nd( fo_pow<T>(), x) );
		else
			return *this;
	}


	this_type operator+ (T x) const { return pwUnaryTr( std::bind2nd( std::plus<T>(), x ) ); }

	this_type operator- (T x) const { return pwUnaryTr( std::bind2nd( std::minus<T>(), x ) ); }

	this_type operator* (T x) const { return pwUnaryTr( std::bind2nd( std::multiplies<T>(), x ) ); }

	this_type operator/ (T x) const { return pwUnaryTr( std::bind2nd( fo_divides0<T>(), x ) ); }

	this_type operator^ (T x) const { return pwUnaryTr( std::bind2nd( fo_pow<T>(), x ) ); }




	template<typename binaryOp> this_type& pwBinaryOp( const this_type &q, binaryOp op ) {

//		std::cout << "\n pwBinaryOp\n";
//		std::cout << "*this = " << *this << "\n\n";
//		std::cout << "q = " << q << "\n";

		DAI_DEBASSERT( size() == q.size() );
		this_type p(*this);

		setDef( op( p.def(), q.def() ) );

		for( typename this_type::const_iterator it = p.begin(); it != p.end(); it++ ) {
			T new_val = op( it->second, q[it->first] );
			set( it->first, new_val );
		}

//		std::cout << "\n *this1 " << *this << "\n";

		 //these are the same, why do they have them?
		for( typename this_type::const_iterator it = q.begin(); it != q.end(); it++ ) {
			T new_val = op( p[it->first], it->second );
			set( it->first, new_val );
		}

//		std::cout << "\n *this2 " << *this << "\n";
//		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++\n";
		return *this;
	}




	template<typename binaryOp> this_type& pwBinaryOpAddSub_eff( const this_type &q, binaryOp op ) {

//		std::cout << "\n pwBinaryOpAddSub_efficient\n";
//		std::cout << "*this = " << *this << "\n\n";
//		std::cout << "q = " << q << "\n";

		DAI_DEBASSERT( size() == q.size() );

		this_type p(*this);

//		(over) allocate memory for elements in the result
		this->p().union_set(p.p(), q.p());

//		std::cout << "container size before: "<< this->p().nrNonDef()<<"\n";

		T new_val;
		for( typename this_type::const_iterator it = begin(); it != end(); it++ ) {
			new_val = op( p[it->first], q[it->first] );
			set( it->first, new_val );
		}

//		std::cout << "container size after: "<< this->p().nrNonDef()<<"\n";

//		std::cout << "\n *this " << *this << "\n";
//		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++\n";
		return *this;
	}



	template<typename binaryOp> this_type& pwBinaryOpMultDiv_eff( const this_type &q, binaryOp op ) {

//		std::cout << "\n pwBinaryOpMultDiv_efficient\n";
//		std::cout << "*this = " << *this << "\n\n";
//		std::cout << "q = " << q << "\n";

		DAI_DEBASSERT( size() == q.size() );

		this_type p(*this);

		//(over) allocate memory for elements in the result
		this->p().intersect(p.p(), q.p());

//		std::cout << "container size before: "<< this->p().nrNonDef()<<"\n";

		T new_val;
		for( typename this_type::const_iterator it = begin(); it != end(); it++ ) {
			new_val = op( p[it->first], q[it->first] );
			set( it->first, new_val );
		}

//		std::cout << "container size after: "<< this->p().nrNonDef()<<"\n";
//
//		std::cout << "\n *this " << *this << "\n";
//		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++\n";
		return *this;
	}





	this_type& operator+= (const this_type & q) { return pwBinaryOp( q, std::plus<T>() ); }


	this_type& operator-= (const this_type & q) { return pwBinaryOp( q, std::minus<T>() ); }


	this_type& operator*= (const this_type & q) { return pwBinaryOp( q, std::multiplies<T>() ); }


	this_type& operator/= (const this_type & q) { return pwBinaryOp( q, fo_divides0<T>() ); }


	this_type& divide (const this_type & q) { return pwBinaryOp( q, std::divides<T>() ); }


	this_type& operator^= (const this_type & q) { return pwBinaryOp( q, fo_pow<T>() ); }




	template<typename binaryOp> this_type pwBinaryTr( const this_type &q, binaryOp op ) const {
//		std::cout << "\n pwBinaryTr\n";
//		std::cout << "*this = " << *this << "\n\n";
//		std::cout << "q = " << q << "\n";

		DAI_DEBASSERT( size() == q.size() );
		this_type result;

		result.setDef( op( def(), q.def() ) );

		result._p.resize( size() );

		for( typename this_type::const_iterator it = begin(); it != end(); it++ ) {
			T new_val = op( it->second, q[it->first] );

			if( new_val != result.def() )
				result._p.push_back(it->first, new_val);
		}

//		std::cout << "\n result1 " << result << "\n";

		for( typename this_type::const_iterator it = q.begin(); it != q.end(); it++ ) {
			T new_val = op( get(it->first), it->second );
			if( new_val != result.def() )
				result.set( it->first, new_val );
		}

//		std::cout << "\n result2 " << result << "\n";
//		std::cout << "-----------------------------------------------\n";

		return result;
	}


	template<typename binaryOp> this_type pwBinaryTrAddSub_eff( const this_type &q, binaryOp op ) const {
//		std::cout << "\n pwBinaryTrAddSub_efficient\n";
//		std::cout << "*this = " << *this << "\n\n";
//		std::cout << "q = " << q << "\n";

		DAI_DEBASSERT( size() == q.size() );

		this_type result;
		result._p.resize( size() );
		result._p.union_set(this->p(), q.p());

//		std::cout << "container size before: "<< result._p.nrNonDef()<<"\n";

		T new_val;
		for( typename this_type::const_iterator it = result.begin(); it != result.end(); it++ ) {
			new_val = op( operator [](it->first), q[it->first] );
			result.set( it->first, new_val );
		}

//		std::cout << "container size after: "<< result._p.nrNonDef()<<"\n";
//
//		std::cout << "\n result " << result << "\n";
//		std::cout << "-----------------------------------------------\n";

		return result;
	}





	template<typename binaryOp> this_type pwBinaryTrMultDiv_eff( const this_type &q, binaryOp op ) const {
//		std::cout << "\n pwBinaryTrMultDiv_efficient\n";
//		std::cout << "*this = " << *this << "\n\n";
//		std::cout << "q = " << q << "\n";

		DAI_DEBASSERT( size() == q.size() );

		this_type result;
		result._p.resize( size() );
		result._p.intersect(this->p(), q.p());

//		std::cout << "container size before: "<< result._p.nrNonDef()<<"\n";

		T new_val;
		for( typename this_type::const_iterator it = result.begin(); it != result.end(); it++ ) {
			new_val = op( operator [](it->first), q[it->first] );
			result.set( it->first, new_val );
		}

//		std::cout << "container size after: "<< result._p.nrNonDef()<<"\n";
//
//		std::cout << "\n result " << result << "\n";
//		std::cout << "-----------------------------------------------\n";

		return result;
	}






	this_type operator+ ( const this_type& q ) const { return pwBinaryTr( q, std::plus<T>() ); }


	this_type operator- ( const this_type& q ) const { return pwBinaryTr( q, std::minus<T>() ); }


	this_type operator* ( const this_type &q ) const { return pwBinaryTr( q, std::multiplies<T>() ); }


	this_type operator/ ( const this_type &q ) const { return pwBinaryTr( q, fo_divides0<T>() ); }


	this_type divided_by( const this_type &q ) const { return pwBinaryTr( q, std::divides<T>() ); }


	this_type operator^ ( const this_type &q ) const { return pwBinaryTr( q, fo_pow<T>() ); }


	template<typename binOp1, typename binOp2> T innerProduct( const this_type &q, T init, binOp1 binaryOp1, binOp2 binaryOp2 ) const {
		DAI_DEBASSERT( size() == q.size() );
		// OPTIMIZE ME
		T result = init;
		for( size_t i = 0; i < size(); i++ )
			result = binaryOp1( result, binaryOp2( get(i), q.get(i) ) );
		return result;
	}
};



template<typename T, typename spvector_type> T dist( const TProb<T,spvector_type>& p, const TProb<T,spvector_type>& q, ProbDistType dt ) {
	switch( dt ) {
	case DISTL1:
		return p.innerProduct( q, (T)0, std::plus<T>(), fo_absdiff<T>() );
	case DISTLINF:
		return p.innerProduct( q, (T)0, fo_max<T>(), fo_absdiff<T>() );
	case DISTTV:
		return p.innerProduct( q, (T)0, std::plus<T>(), fo_absdiff<T>() ) / 2;
	case DISTKL:
		return p.innerProduct( q, (T)0, std::plus<T>(), fo_KL<T>() );
	case DISTHEL:
		return p.innerProduct( q, (T)0, std::plus<T>(), fo_Hellinger<T>() ) / 2;
	default:
		DAI_THROW(UNKNOWN_ENUM_VALUE);
		return INFINITY;
	}
}



template<typename T, typename spvector_type> std::ostream& operator<< (std::ostream& os, const TProb<T,spvector_type>& p) {
	os << p.p();
	return os;
}



template<typename T, typename spvector_type> TProb<T,spvector_type> min( const TProb<T,spvector_type> &a, const TProb<T,spvector_type> &b ) {
	return a.pwBinaryTr( b, fo_min<T>() );
}



template<typename T, typename spvector_type> TProb<T,spvector_type> max( const TProb<T,spvector_type> &a, const TProb<T,spvector_type> &b ) {
	return a.pwBinaryTr( b, fo_max<T>() );
}


/// Represents a vector with entries of type dai::Real.
typedef TProb<Real, spvector<Real> > Prob;


} // end of namespace dai


#endif
