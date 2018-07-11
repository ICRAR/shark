//
// Various utilities for shark
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2017
// Copyright by UWA (in the framework of the ICRAR)
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307  USA
//

#include <algorithm>
#include <cctype>
#include <functional>
#include <locale>
#include <string>
#include <sstream>

#include <errno.h>
#include <string.h>

#include "utils.h"

using namespace std;

namespace shark {

vector<string> tokenize(const string &s, const string &delims)
{
	string::size_type lastPos = s.find_first_not_of(delims, 0);
	string::size_type pos     = s.find_first_of(delims, lastPos);

	vector<string> tokens;
	while (string::npos != pos || string::npos != lastPos) {
		tokens.push_back(s.substr(lastPos, pos - lastPos));
		lastPos = s.find_first_not_of(delims, pos);
		pos = s.find_first_of(delims, lastPos);
	}
	return tokens;
}

// trim from start
static inline string &ltrim(std::string &s) {
	s.erase(s.begin(), find_if(s.begin(), s.end(),
	        not1(ptr_fun<int, int>(isspace))));
	return s;
}

// trim from end
static inline string &rtrim(std::string &s) {
	s.erase(find_if(s.rbegin(), s.rend(),
	        not1(ptr_fun<int, int>(isspace))).base(), s.end());
	return s;
}

// trim from both ends
void trim(std::string &s) {
	ltrim(rtrim(s));
}

void lower(string &s)
{
	transform(s.begin(), s.end(), s.begin(), ::tolower);
}

std::string lower(const std::string &s)
{
	std::string low(s);
	lower(low);
	return low;
}

void upper(string &s)
{
	transform(s.begin(), s.end(), s.begin(), ::toupper);
}

ifstream open_file(const string &name)
{
	ifstream f(name);
	if( !f ) {
		ostringstream os;
		os << "Error when opening file '" << name << "': " << strerror(errno);
		throw runtime_error(os.str());
	}
	return f;
}

bool empty_or_comment(const std::string &s) {
	return s.size() == 0 or s[0] == '#';
}

}  // namespace shark
