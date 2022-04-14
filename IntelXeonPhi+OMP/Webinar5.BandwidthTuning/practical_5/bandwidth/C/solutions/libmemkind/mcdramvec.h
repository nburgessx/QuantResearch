#ifndef MCDRAMVEC_H_
#define MCDRAMVEC_H_
#include <vector>
#include <hbw_allocator.h>

namespace mcdramvec {


template <typename real_t> using vector = std::vector<real_t,hbw::allocator<real_t> >;


};

#endif
