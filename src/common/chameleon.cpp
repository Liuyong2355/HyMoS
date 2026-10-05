/* 
   Chameleon.cpp

   Copyright (C) 2002-2004 René Nyffenegger

   This source code is provided 'as-is', without any express or implied
   warranty. In no event will the author be held liable for any damages
   arising from the use of this software.

   Permission is granted to anyone to use this software for any purpose,
   including commercial applications, and to alter it and redistribute it
   freely, subject to the following restrictions:

   1. The origin of this source code must not be misrepresented; you must not
      claim that you wrote the original source code. If you use this source code
      in a product, an acknowledgment in the product documentation would be
      appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
      misrepresented as being the original source code.

   3. This notice may not be removed or altered from any source distribution.

   René Nyffenegger rene.nyffenegger@adp-gmbh.ch
*/

// Modified for HyMoS: integer conversion support and interface adaptations.

#include "chameleon.h"

#include <cstdlib>
#include <sstream>

Chameleon::Chameleon(const std::string& value) : value_(value) {}

Chameleon::Chameleon(const char* value) : value_(value) {}

Chameleon::Chameleon(double value) {
  std::stringstream stream;
  stream << value;
  value_ = stream.str();
}

Chameleon::Chameleon(int value) {
  std::stringstream stream;
  stream << value;
  value_ = stream.str();
}

Chameleon::Chameleon(const Chameleon& other) : value_(other.value_) {}

Chameleon& Chameleon::operator=(const Chameleon& other) {
  value_ = other.value_;
  return *this;
}

Chameleon& Chameleon::operator=(double value) {
  std::stringstream stream;
  stream << value;
  value_ = stream.str();
  return *this;
}

Chameleon& Chameleon::operator=(int value) {
  std::stringstream stream;
  stream << value;
  value_ = stream.str();
  return *this;
}

Chameleon& Chameleon::operator=(const std::string& value) {
  value_ = value;
  return *this;
}

Chameleon::operator std::string() const {
  return value_;
}

Chameleon::operator double() const {
  return atof(value_.c_str());
}

Chameleon::operator int() const {
  return atoi(value_.c_str());
}
