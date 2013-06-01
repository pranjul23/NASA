
#ifndef __defined_libdai_spvector_boost_h
#define __defined_libdai_spvector_boost_h


#include <iostream>
#include <ostream>
#include <algorithm>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace dai {

using namespace boost::numeric::ublas;

template <typename T>
class spvector_boost {
public:

	typedef compressed_vector<T> nondefaults_type;

private:
	nondefaults_type _p;
	size_t _size;
	T _def;
public:


	spvector_boost() : _p(), _size(0){}

	explicit spvector_boost( size_t n, T p ) : _p(n,n), _size(n){
		if(p != (T)0){
			for(size_t i=0; i<n; i++){
				_p[i]=p;
			}
		}
	}


	template <typename TIterator>
	spvector_boost( TIterator begin, TIterator end, size_t sizeHint) : _p(sizeHint, sizeHint), _size(0){
		size_t iter = 0;
		for( TIterator it = begin; it != end; it++, iter++ ){
			if( *it != _def )
				push_back( iter, *it );
		}
		_size = iter;

		compressed_vector<T> gg;
		gg.nnz() .value_data()
	}


	template <typename S>
	spvector_boost( const std::vector<S> &v, size_t sizeHint) : _p(sizeHint, sizeHint), _size(v.size()), _def((T)0) {
		for( size_t i = 0; i < v.size(); i++ )
			if( v[i] != _def )
				push_back( i, v[i] );
	}



	typedef typename nondefaults_type::const_iterator const_iterator;
	typedef typename nondefaults_type::iterator iterator;
	typedef typename nondefaults_type::const_reverse_iterator const_reverse_iterator;
	typedef typename nondefaults_type::reverse_iterator reverse_iterator;

	iterator begin() { return _p.begin(); }
	const_iterator begin() const { return _p.begin(); }

	iterator end() { return _p.end(); }
	const_iterator end() const { return _p.end(); }

	reverse_iterator rbegin() { return _p.rbegin(); }
	const_reverse_iterator rbegin() const { return _p.rbegin(); }

	reverse_iterator rend() { return _p.rend(); }
	const_reverse_iterator rend() const { return _p.rend(); }



	void reserve( size_t n ) {
		_p.reserve(n);
	}

	size_t size() const { return _size; }

	void resize( size_t n ) {
		_size = n;
		if( _p.size() ) {
			for( iterator it = begin(); it != end(); ) {
				if( it->first >= n )
					_p.erase( it++ );
				else
					it++;
			}
		}
	}

	iterator erase( iterator position ) {
		iterator next = position;
		next++;
		_p.erase_element( position );
		return next;
	}

	void push_back( size_t idx, T val ) { _p[idx] = val; }

	void setDefault( T def ) {
		DAI_ASSERT( def == (T)0 );
		_def = def;
	}

	T getDefault() const { return _def; }

	void clearNonDef() { _p.clear(); }

	void set( size_t i, T val ) {
		DAI_ASSERT( i < _size );
		if( val != _def )
			_p[i] = val;
		else
			_p.erase( i );
	}

	T get( size_t i ) const {
		DAI_ASSERT( i < _size );
		const_iterator it = _p.find( i );
		if( it != _p.end() )
			return it->second;
		else
			return _def;
	}

	T operator[]( size_t i ) const { return get(i); }

	bool operator==( const spvector_map<T>& q ) const {
		if( size() != q.size() )
			return false;
		// OPTIMIZE
		for( size_t i = 0; i < size(); i++ )
			if( !(get(i) == q.get(i)) )
				return false;
		return true;
	}

	size_t nrNonDef() const { return _p.size(); }

	size_t nrDef() const { return _size - _p.size(); }

	T def() const {
		return _def;
	}

	void setDef( T def ) {
		DAI_ASSERT( def == (T)0 );
		_def = def;
	}

	const nondefaults_type & nonDef() const { return _p; }


	//intersect p and q and store result in *this:


	void intersect(const spvector_map<T>& p, const spvector_map<T>& q){
		_p.clear();
		std::set_intersection(p.begin(), p.end(), q.begin(), q.end(), std::inserter(_p,  _p.begin()),  _p.value_comp());
	}

	void union_set(const spvector_map<T>& p, const spvector_map<T>& q){
		std::set_union(p.begin(), p.end(), q.begin(), q.end(), std::inserter(_p,  _p.begin()),  _p.value_comp());
	}

};


template<class T>
std::ostream& operator << (std::ostream& os, const spvector_map<T> &x) {
	os << "(";
	os << "size:" << x.size();
	os << ", def:" << x.def();
	for( typename spvector_map<T>::const_iterator it = x.begin(); it != x.end(); it++ )
		os << ", " << it->first << ":" << it->second;
	os << ")";
	return os;
}


} // end of namespace dai


#endif
