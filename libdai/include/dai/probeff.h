
#ifndef __defined_libdai_probeff_h
#define __defined_libdai_probeff_h


#include <cmath>
#include <vector>
#include <ostream>
#include <algorithm>
#include <numeric>
#include <functional>
#include <dai/util.h>
#include <dai/exceptions.h>
#include <dai/fo.h>


namespace dai {

template <typename T>
class TProbEff {
    public:
        /// Type of data structure used for storing the values
        typedef std::vector<std::vector<T> > container_type;

        /// Shorthand
        typedef TProbEff<T> this_type;

    private:
        /// The data structure that stores the values
        container_type _p;

    public:
    /// \name Constructors and destructors
    //@{
        /// Default constructor (constructs empty vector)
        TProbEff() : _p() {}

        //construct matrix of size n1 x n2
        explicit TProbEff( size_t n1, size_t n2 ) : _p(n1, std::vector<T>(n2,1)) {}


        /// Constant iterator over the elements
        typedef typename container_type::const_iterator const_iterator;
        /// Iterator over the elements
        typedef typename container_type::iterator iterator;
        /// Constant reverse iterator over the elements
        typedef typename container_type::const_reverse_iterator const_reverse_iterator;
        /// Reverse iterator over the elements
        typedef typename container_type::reverse_iterator reverse_iterator;


    /// \name Iterator interface
    //@{
        /// Returns iterator that points to the first element
        iterator begin() { return _p.begin(); }
        /// Returns constant iterator that points to the first element
        const_iterator begin() const { return _p.begin(); }

        /// Returns iterator that points beyond the last element
        iterator end() { return _p.end(); }
        /// Returns constant iterator that points beyond the last element
        const_iterator end() const { return _p.end(); }

        /// Returns reverse iterator that points to the last element
        reverse_iterator rbegin() { return _p.rbegin(); }
        /// Returns constant reverse iterator that points to the last element
        const_reverse_iterator rbegin() const { return _p.rbegin(); }

        /// Returns reverse iterator that points beyond the first element
        reverse_iterator rend() { return _p.rend(); }
        /// Returns constant reverse iterator that points beyond the first element
        const_reverse_iterator rend() const { return _p.rend(); }
    //@}

    /// \name Get/set individual entries
    //@{
        /// Gets \a i 'th entry
        T get( size_t i, size_t j ) const {
            return _p[i][j];
        }

        /// Sets \a i 'th entry to \a val
        void set( size_t i, size_t j, T val ) {
            DAI_DEBASSERT( i < _p.size() );
            DAI_DEBASSERT( j < _p[i].size() );
            _p[i][j] = val;
        }

    /// \name Queries
    //@{
        /// Returns a const reference to the wrapped container
        const container_type& p() const { return _p; }

        /// Returns a reference to the wrapped container
        container_type& p() { return _p; }

        /// Accumulate all values (similar to std::accumulate) by summing
        /** The following calculation is done:
         *  \code
         *  T t = op(init);
         *  for( const_iterator it = begin(); it != end(); it++ )
         *      t += op(*it);
         *  return t;
         *  \endcode
         */
        template<typename unOp> T accumulateSum( T init, unOp op ) const {
            T t = op(init);
            for( const_iterator it = begin(); it != end(); it++ )
                t += op(*it);
            return t;
        }

        /// Accumulate all values (similar to std::accumulate) by maximization/minimization
        /** The following calculation is done (with "max" replaced by "min" if \a minimize == \c true):
         *  \code
         *  T t = op(init);
         *  for( const_iterator it = begin(); it != end(); it++ )
         *      t = std::max( t, op(*it) );
         *  return t;
         *  \endcode
         */
        template<typename unOp> T accumulateMax( T init, unOp op, bool minimize ) const {
            T t = op(init);
            if( minimize ) {
                for( const_iterator it = begin(); it != end(); it++ )
                    t = std::min( t, op(*it) );
            } else {
                for( const_iterator it = begin(); it != end(); it++ )
                    t = std::max( t, op(*it) );
            }
            return t;
        }

        /// Returns the Shannon entropy of \c *this, \f$-\sum_i p_i \log p_i\f$
        T entropy() const { return -accumulateSum( (T)0, fo_plog0p<T>() ); }

        /// Returns maximum value of all entries
        T max() const { return accumulateMax( (T)(-INFINITY), fo_id<T>(), false ); }

        /// Returns minimum value of all entries
        T min() const { return accumulateMax( (T)INFINITY, fo_id<T>(), true ); }

        /// Returns sum of all entries
        T sum() const { return accumulateSum( (T)0, fo_id<T>() ); }

        /// Return sum of absolute value of all entries
        T sumAbs() const { return accumulateSum( (T)0, fo_abs<T>() ); }

        /// Returns maximum absolute value of all entries
        T maxAbs() const { return accumulateMax( (T)0, fo_abs<T>(), false ); }

        /// Returns \c true if one or more entries are NaN
        bool hasNaNs() const {
            bool foundnan = false;
            for( const_iterator x = _p.begin(); x != _p.end(); x++ )
                if( dai::isnan( *x ) ) {
                    foundnan = true;
                    break;
                }
            return foundnan;
        }

        /// Returns \c true if one or more entries are negative
        bool hasNegatives() const {
            return (std::find_if( _p.begin(), _p.end(), std::bind2nd( std::less<T>(), (T)0 ) ) != _p.end());
        }

        /// Returns a pair consisting of the index of the maximum value and the maximum value itself
        std::pair<size_t,T> argmax() const {
            T max = _p[0];
            size_t arg = 0;
            for( size_t i = 1; i < size(); i++ ) {
              if( _p[i] > max ) {
                max = _p[i];
                arg = i;
              }
            }
            return std::make_pair( arg, max );
        }

        /// Returns a random index, according to the (normalized) distribution described by *this
        size_t draw() {
            Real x = rnd_uniform() * sum();
            T s = 0;
            for( size_t i = 0; i < size(); i++ ) {
                s += get(i);
                if( s > x )
                    return i;
            }
            return( size() - 1 );
        }

        /// Lexicographical comparison
        /** \pre <tt>this->size() == q.size()</tt>
         */
        bool operator<( const this_type& q ) const {
            DAI_DEBASSERT( size() == q.size() );
            return lexicographical_compare( begin(), end(), q.begin(), q.end() );
        }

        /// Comparison
        bool operator==( const this_type& q ) const {
            if( size() != q.size() )
                return false;
            return p() == q.p();
        }
    //@}

    /// \name Unary transformations
    //@{
        /// Returns the result of applying operation \a op pointwise on \c *this
        template<typename unaryOp> this_type pwUnaryTr( unaryOp op ) const {
            this_type r;
            r._p.reserve( size() );
            std::transform( _p.begin(), _p.end(), back_inserter( r._p ), op );
            return r;
        }

        /// Returns negative of \c *this
        this_type operator- () const { return pwUnaryTr( std::negate<T>() ); }

        /// Returns pointwise absolute value
        this_type abs() const { return pwUnaryTr( fo_abs<T>() ); }

        /// Returns pointwise exponent
        this_type exp() const { return pwUnaryTr( fo_exp<T>() ); }

        /// Returns pointwise logarithm
        /** If \a zero == \c true, uses <tt>log(0)==0</tt>; otherwise, <tt>log(0)==-Inf</tt>.
         */
        this_type log(bool zero=false) const {
            if( zero )
                return pwUnaryTr( fo_log0<T>() );
            else
                return pwUnaryTr( fo_log<T>() );
        }

        /// Returns pointwise inverse
        /** If \a zero == \c true, uses <tt>1/0==0</tt>; otherwise, <tt>1/0==Inf</tt>.
         */
        this_type inverse(bool zero=true) const {
            if( zero )
                return pwUnaryTr( fo_inv0<T>() );
            else
                return pwUnaryTr( fo_inv<T>() );
        }

        /// Returns normalized copy of \c *this, using the specified norm
        /** \throw NOT_NORMALIZABLE if the norm is zero
         */
        this_type normalized( ProbNormType norm = dai::NORMPROB ) const {
            T Z = 0;
            if( norm == dai::NORMPROB )
                Z = sum();
            else if( norm == dai::NORMLINF )
                Z = maxAbs();
            if( Z == (T)0 ) {
                DAI_THROW(NOT_NORMALIZABLE);
                return *this;
            } else
                return pwUnaryTr( std::bind2nd( std::divides<T>(), Z ) );
        }
    //@}

    /// \name Unary operations
    //@{
        /// Applies unary operation \a op pointwise
        template<typename unaryOp> this_type& pwUnaryOp( unaryOp op ) {
            std::transform( _p.begin(), _p.end(), _p.begin(), op );
            return *this;
        }

        /// Draws all entries i.i.d. from a uniform distribution on [0,1)
        this_type& randomize() {
            std::generate( _p.begin(), _p.end(), rnd_uniform );
            return *this;
        }

        /// Sets all entries to \f$1/n\f$ where \a n is the length of the vector
        this_type& setUniform () {
            fill( (T)1 / size() );
            return *this;
        }

        /// Applies absolute value pointwise
        this_type& takeAbs() { return pwUnaryOp( fo_abs<T>() ); }

        /// Applies exponent pointwise
        this_type& takeExp() { return pwUnaryOp( fo_exp<T>() ); }

        /// Applies logarithm pointwise
        /** If \a zero == \c true, uses <tt>log(0)==0</tt>; otherwise, <tt>log(0)==-Inf</tt>.
         */
        this_type& takeLog(bool zero=false) {
            if( zero ) {
                return pwUnaryOp( fo_log0<T>() );
            } else
                return pwUnaryOp( fo_log<T>() );
        }

        /// Normalizes vector using the specified norm
        /** \throw NOT_NORMALIZABLE if the norm is zero
         */
        T normalize( ProbNormType norm=dai::NORMPROB ) {
            T Z = 0;
            if( norm == dai::NORMPROB )
                Z = sum();
            else if( norm == dai::NORMLINF )
                Z = maxAbs();

            if( Z == (T)0 )
                //DAI_THROW(NOT_NORMALIZABLE);
            	//if Z is zero, it means all elements are 0
            	//we can let the normalized distrib to be all 0 as well
            	*this /= (T)1;
            else
                *this /= Z;
            return Z;
        }
    //@}

    /// \name Operations with scalars
    //@{
        /// Sets all entries to \a x
        this_type& fill( T x ) {
            std::fill( _p.begin(), _p.end(), x );
            return *this;
        }

        /// Adds scalar \a x to each entry
        this_type& operator+= (T x) {
            if( x != 0 )
                return pwUnaryOp( std::bind2nd( std::plus<T>(), x ) );
            else
                return *this;
        }

        /// Subtracts scalar \a x from each entry
        this_type& operator-= (T x) {
            if( x != 0 )
                return pwUnaryOp( std::bind2nd( std::minus<T>(), x ) );
            else
                return *this;
        }

        /// Multiplies each entry with scalar \a x
        this_type& operator*= (T x) {
            if( x != 1 )
                return pwUnaryOp( std::bind2nd( std::multiplies<T>(), x ) );
            else
                return *this;
        }

        /// Divides each entry by scalar \a x, where division by 0 yields 0
        this_type& operator/= (T x) {
            if( x != 1 )
                return pwUnaryOp( std::bind2nd( fo_divides0<T>(), x ) );
            else
                return *this;
        }

        /// Raises entries to the power \a x
        this_type& operator^= (T x) {
            if( x != (T)1 )
                return pwUnaryOp( std::bind2nd( fo_pow<T>(), x) );
            else
                return *this;
        }
    //@}

    /// \name Transformations with scalars
    //@{
        /// Returns sum of \c *this and scalar \a x
        this_type operator+ (T x) const { return pwUnaryTr( std::bind2nd( std::plus<T>(), x ) ); }

        /// Returns difference of \c *this and scalar \a x
        this_type operator- (T x) const { return pwUnaryTr( std::bind2nd( std::minus<T>(), x ) ); }

        /// Returns product of \c *this with scalar \a x
        this_type operator* (T x) const { return pwUnaryTr( std::bind2nd( std::multiplies<T>(), x ) ); }

        /// Returns quotient of \c *this and scalar \a x, where division by 0 yields 0
        this_type operator/ (T x) const { return pwUnaryTr( std::bind2nd( fo_divides0<T>(), x ) ); }

        /// Returns \c *this raised to the power \a x
        this_type operator^ (T x) const { return pwUnaryTr( std::bind2nd( fo_pow<T>(), x ) ); }
    //@}

    /// \name Operations with other equally-sized vectors
    //@{
        /// Applies binary operation pointwise on two vectors
        /** \tparam binaryOp Type of function object that accepts two arguments of type \a T and outputs a type \a T
         *  \param q Right operand
         *  \param op Operation of type \a binaryOp
         */
        template<typename binaryOp> this_type& pwBinaryOp( const this_type &q, binaryOp op ) {
            DAI_DEBASSERT( size() == q.size() );
            std::transform( _p.begin(), _p.end(), q._p.begin(), _p.begin(), op );
            return *this;
        }

        /// Pointwise addition with \a q
        /** \pre <tt>this->size() == q.size()</tt>
         */
        this_type& operator+= (const this_type & q) { return pwBinaryOp( q, std::plus<T>() ); }

        /// Pointwise subtraction of \a q
        /** \pre <tt>this->size() == q.size()</tt>
         */
        this_type& operator-= (const this_type & q) { return pwBinaryOp( q, std::minus<T>() ); }

        /// Pointwise multiplication with \a q
        /** \pre <tt>this->size() == q.size()</tt>
         */
        this_type& operator*= (const this_type & q) { return pwBinaryOp( q, std::multiplies<T>() ); }

        /// Pointwise division by \a q, where division by 0 yields 0
        /** \pre <tt>this->size() == q.size()</tt>
         *  \see divide(const TProb<T> &)
         */
        this_type& operator/= (const this_type & q) { return pwBinaryOp( q, fo_divides0<T>() ); }

        /// Pointwise division by \a q, where division by 0 yields +Inf
        /** \pre <tt>this->size() == q.size()</tt>
         *  \see operator/=(const TProb<T> &)
         */
        this_type& divide (const this_type & q) { return pwBinaryOp( q, std::divides<T>() ); }

        /// Pointwise power
        /** \pre <tt>this->size() == q.size()</tt>
         */
        this_type& operator^= (const this_type & q) { return pwBinaryOp( q, fo_pow<T>() ); }
    //@}

    /// \name Transformations with other equally-sized vectors
    //@{
        /// Returns the result of applying binary operation \a op pointwise on \c *this and \a q
        /** \tparam binaryOp Type of function object that accepts two arguments of type \a T and outputs a type \a T
         *  \param q Right operand
         *  \param op Operation of type \a binaryOp
         */
        template<typename binaryOp> this_type pwBinaryTr( const this_type &q, binaryOp op ) const {
            DAI_DEBASSERT( size() == q.size() );
            TProb<T> r;
            r._p.reserve( size() );
            std::transform( _p.begin(), _p.end(), q._p.begin(), back_inserter( r._p ), op );
            return r;
        }

        /// Returns sum of \c *this and \a q
        /** \pre <tt>this->size() == q.size()</tt>
         */
        this_type operator+ ( const this_type& q ) const { return pwBinaryTr( q, std::plus<T>() ); }

        /// Return \c *this minus \a q
        /** \pre <tt>this->size() == q.size()</tt>
         */
        this_type operator- ( const this_type& q ) const { return pwBinaryTr( q, std::minus<T>() ); }

        /// Return product of \c *this with \a q
        /** \pre <tt>this->size() == q.size()</tt>
         */
        this_type operator* ( const this_type &q ) const { return pwBinaryTr( q, std::multiplies<T>() ); }

        /// Returns quotient of \c *this with \a q, where division by 0 yields 0
        /** \pre <tt>this->size() == q.size()</tt>
         *  \see divided_by(const TProb<T> &)
         */
        this_type operator/ ( const this_type &q ) const { return pwBinaryTr( q, fo_divides0<T>() ); }

        /// Pointwise division by \a q, where division by 0 yields +Inf
        /** \pre <tt>this->size() == q.size()</tt>
         *  \see operator/(const TProb<T> &)
         */
        this_type divided_by( const this_type &q ) const { return pwBinaryTr( q, std::divides<T>() ); }

        /// Returns \c *this to the power \a q
        /** \pre <tt>this->size() == q.size()</tt>
         */
        this_type operator^ ( const this_type &q ) const { return pwBinaryTr( q, fo_pow<T>() ); }
    //@}

        /// Performs a generalized inner product, similar to std::inner_product
        /** \pre <tt>this->size() == q.size()</tt>
         */
        template<typename binOp1, typename binOp2> T innerProduct( const this_type &q, T init, binOp1 binaryOp1, binOp2 binaryOp2 ) const {
            DAI_DEBASSERT( size() == q.size() );
            return std::inner_product( begin(), end(), q.begin(), init, binaryOp1, binaryOp2 );
        }
};


/// Returns distance between \a p and \a q, measured using distance measure \a dt
/** \relates TProb
 *  \pre <tt>this->size() == q.size()</tt>
 */
template<typename T> T dist( const TProb<T> &p, const TProb<T> &q, ProbDistType dt ) {
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


/// Writes a TProb<T> to an output stream
/** \relates TProb
 */
template<typename T> std::ostream& operator<< (std::ostream& os, const TProb<T>& p) {
    os << "(";
    for( size_t i = 0; i < p.size(); i++ )
        os << ((i != 0) ? ", " : "") << p.get(i);
    os << ")";
    return os;
}


/// Returns the pointwise minimum of \a a and \a b
/** \relates TProb
 *  \pre <tt>this->size() == q.size()</tt>
 */
template<typename T> TProb<T> min( const TProb<T> &a, const TProb<T> &b ) {
    return a.pwBinaryTr( b, fo_min<T>() );
}


/// Returns the pointwise maximum of \a a and \a b
/** \relates TProb
 *  \pre <tt>this->size() == q.size()</tt>
 */
template<typename T> TProb<T> max( const TProb<T> &a, const TProb<T> &b ) {
    return a.pwBinaryTr( b, fo_max<T>() );
}


/// Represents a vector with entries of type dai::Real.
typedef TProb<Real> Prob;


} // end of namespace dai


#endif