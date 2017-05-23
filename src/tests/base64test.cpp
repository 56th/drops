#include "misc/base64.h"

#include <string>
#include<algorithm>
#include <vector>
#include <iterator>
#include <sstream>

int main ()
{
    std::cout << "string? ";
    std::string s;
    std::getline( std::cin, s);
    std::cout << "base64:\n";
    unsigned char* p= (unsigned char*) &s[0];
    DROPS::Base64Encoding::encode( p, p + s.size(), std::cout, /*wrap_lines*/ true);
    std::cout << std::endl;

    std::cout << "base64? ";
    std::getline( std::cin, s);
    std::istringstream iss( s);
    std::vector<unsigned char> v;
    DROPS::Base64Encoding::decode( iss, std::back_inserter<std::vector<unsigned char> >( v));
    for (size_t i= 0; i < v.size(); ++i)
        std::cerr << v[i] << std::flush;
    std::cout << std::endl;
    return 0;
}