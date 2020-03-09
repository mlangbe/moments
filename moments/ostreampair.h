/*
    output stream being forwarded to two other output streams.

	This file is part of the source code used for the calculation of the moment invariants
	as described in the dissertation of Max Langbein
	https://nbn-resolving.org/urn:nbn:de:hbz:386-kluedo-38558

	Copyright (C) 2020 TU Kaiserslautern,
	              Prof.Dr. Hans Hagen (AG Computergraphik und HCI) (hagen@cs.uni-kl.de)

	Author: Max Langbein (m_langbe@cs.uni-kl.de)

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef OSTREAM_PAIR_H
#define OSTREAM_PAIR_H

/*
this file defines 
basic_ostreampair<T> and ostreampair,
which is an ostream that forwards
anything written to it
to two other ostreams.
*/


/*
* the streambuffer that does the work.
*/
template<class T>
class ostreampairbuf:public basic_streambuf<T>
{
	
	basic_ostream<T>&a;	
	basic_ostream<T>&b;

	//a small buffer is needed because 
	//the virtual function overflow shouldn't be called too often.
	static const int nbuf= 32;
	T buf[nbuf];


	void flushme()
	{
		for(T * s = pbase();s!=pptr();++s)
		{
			a.put(*s);b.put(*s);
		}
		
		setp(buf,epptr());	
	}



	virtual int_type overflow(int_type c)
	{
		
		flushme();

		a.put(traits_type::to_char_type(c));
		b.put(traits_type::to_char_type(c));

		return traits_type::not_eof(c);
	}

	virtual int sync()
	{
		flushme();
		a.flush();
		b.flush();
		return 0;
	}

public :
	ostreampairbuf(basic_ostream<T> & a,
				 basic_ostream<T> & b )
		:a(a),b(b)
	{		
		setp(buf,buf+nbuf);
	}

};








/*
a basic_ostream that forwards
anything written to it
to two other ostreams.
*/
template<class T=char>
class basic_ostreampair :public basic_ostream<T>
{
	ostreampairbuf<T> s;

public:

	basic_ostreampair(basic_ostream<T> & a,
				     basic_ostream<T> & b )
		:s(a,b),basic_ostream<T>(&s)
	{
	}
};

typedef basic_ostreampair<char> ostreampair;

#endif
