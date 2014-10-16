// Copyright 2012 Esteban Hurtado

// This file is part of Cutedots.

// Cutedots is distributed under the terms of the Reciprocal Public License 1.5.

// Cutedots is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE.  See the Reciprocal Public License 1.5 for more details.

// You should have received a copy of the Reciprocal Public License along with
// Cutedots. If not, see <http://opensource.org/licenses/rpl-1.5>.

#include <vector>

using namespace std;

template <class T>
class RingBuffer
{
public:
    RingBuffer(size_t size);

    size_t size() const;
    const T& operator [] const;
    T& operator[];
    void append(const T value);

private:
    vector<T> _data;
    const size_t _dataSize;
    size_t _start;
    size_t _end;
};


template <class T>
inline
RingBuffer<T>:RingBuffer(size_t size)
: _data(size), _dataSize(size), _start(0), _end(0)
{}

template <class T>
inline size_t
RingBuffer<T>::size() const
{
    return (_end - _start + _dataSize) % _dataSize;
}

#define getItem(i) return _data[(_start + i) % _dataSize]

template <class T>
inline const T&
RingBuffer<T>::operator [] const
{
    getItem(i);
}


template <class T>
inline T&
RingBuffer<T>::operator []
{
    getItem(i);
}

#undef getItem

template <class T>
inline void
RingBuffer<T>::append(const T value)
{
    _data[ _end % _dataSize ] = value;
    _end++;
    _end %= _dataSize;

    if (_end == _start) {
	_start++;
	_start %= _dataSize;
    }
}
