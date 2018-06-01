
#include <BL_Arena.H>

const unsigned int bl::Arena::align_size;

bl::Arena::~Arena () {}

std::size_t
bl::Arena::align (std::size_t s)
{
    std::size_t x = s + (align_size-1);
    x -= x & (align_size-1);
    return x;
}
