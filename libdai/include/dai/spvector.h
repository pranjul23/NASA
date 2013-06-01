#ifndef __defined_libdai_spvector_h
#define __defined_libdai_spvector_h


#include <ostream>
#include <vector>
#include <algorithm>
#include <iostream>



template <typename T1, typename T2>
struct first_less : public std::binary_function<const std::pair<T1,T2> &, const std::pair<T1,T2> &, bool> {
	bool operator()( const std::pair<T1,T2> &a, const std::pair<T1,T2> &b ) { return a.first < b.first; }
};


template <typename T>
class spvector {
public:
	typedef std::pair<size_t, T> nondefault_type;

	typedef std::vector<nondefault_type> nondefaults_type;

private:

	nondefaults_type _p;
	size_t _size;
	T _def;

public:


	spvector() : _p(), _size(0), _def((T)0) {}

	explicit spvector( size_t n, T p ) : _size(n), _def((T)0) {

		if(p != (T)0){
			_p.reserve(n);

			for(size_t i=0; i<n; i++){
				push_back(i, p);
			}
		}
	}


	explicit spvector( size_t n, size_t cap ) : _def((T)0) {
		_size = n;
		_p.reserve(cap);
	}



	template <typename TIterator>
	spvector( TIterator begin, TIterator end, size_t sizeHint) : _p(), _size(0), _def((T)0) {
		if( sizeHint )
			reserve( sizeHint );
		size_t iter = 0;
		for( TIterator it = begin; it != end; it++, iter++ ){
			if( *it != _def )
				push_back( iter, *it );
		}
		_size = iter;
	}


	template <typename S>
	spvector( const std::vector<S> &v, size_t sizeHint) : _p(), _size(v.size()), _def((T)0) {
		if( sizeHint )
			reserve( sizeHint );
		for( size_t i = 0; i < v.size(); i++ )
			if( v[i] != def )
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



	void reserve( size_t n ) { _p.reserve( n ); }

	size_t size() const { return _size; }

	void resize( size_t n ) {
		_size = n;
		if( _p.size() ) {
			if( _p.back().first >= n ) {
				iterator it = lower_bound( _p.begin(), _p.end(), std::make_pair(n, T()), first_less<size_t, T>() );
				_p.resize( distance( _p.begin(), it ) );
			}
		}
	}

	iterator erase( iterator position ) {
		return _p.erase( position );
	}

	void push_back( size_t idx, T val ) { _p.push_back( std::make_pair( idx, val ) ); }

	void setDefault( T def ) {
		DAI_ASSERT( def == (T)0 );
		_def = def;
	}

	T getDefault() const { return _def; }

	void clearNonDef() { _p.clear(); }

	void set( size_t i, T val ) {
		DAI_ASSERT( i < _size );
		iterator it = lower_bound( _p.begin(), _p.end(), std::make_pair(i, T()), first_less<size_t, T>() );
		if( (it != _p.end()) && (it->first == i) ) {
			// nondefault value already present
			if( val == _def )
				_p.erase( it );
			else
				it->second = val;
		} else {
			// no nondefault value present yet
			if( val == _def )
				; // do nothing
			else
				_p.insert( it, std::make_pair(i, val) );
		}
	}

	T get( size_t i ) const {
		DAI_ASSERT( i < _size );
		const_iterator it = lower_bound( _p.begin(), _p.end(), std::make_pair(i, T()), first_less<size_t, T>() );
		if( (it != _p.end()) && (it->first == i) )
			return it->second;
		else
			return _def;
	}

	T operator[]( size_t i ) const { return get(i); }

	bool operator==( const spvector<T>& q ) const {
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

	T def() const { return _def; }

	void setDef( T def ) {
		DAI_ASSERT( def == (T)0 );
		_def = def;
	}

	const nondefaults_type & nonDef() const { return _p; }


	void intersect(const spvector<T>& p, const spvector<T>& q){
		_p.clear();
		_p.reserve( std::min(p.nrNonDef(), q.nrNonDef()) );

		std::set_intersection(p.begin(), p.end(), q.begin(), q.end(), std::inserter(_p, _p.begin()), first_less<size_t, T>());
	}

	void union_set(const spvector<T>& p, const spvector<T>& q){
		_p.clear();
		_p.reserve( p.nrNonDef()+q.nrNonDef() );

		std::set_union(p.begin(), p.end(), q.begin(), q.end(),  std::inserter(_p, _p.begin()),  first_less<size_t, T>());
	}

};


template<class T>
std::ostream& operator<< (std::ostream& os, const spvector<T> &x) {
	os << "(";
	os << "size:" << x.size();
	os << ", def:" << x.def();
	for( typename spvector<T>::const_iterator it = x.begin(); it != x.end(); it++ )
		os << ", " << it->first << ":" << it->second;
	os << ")";
	return os;
}

#endif
