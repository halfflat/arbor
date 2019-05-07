#include "../gtest.h"

#include <vector>

#include <backends/gpu/reduce_by_key.hpp>
#include <memory/memory.hpp>
#include <util/span.hpp>
#include <util/rangeutil.hpp>

using namespace arb;

template <typename T, typename I>
__global__
void reduce_kernel(const T* src, T* dst, const I* index, int n) {
    unsigned tid = threadIdx.x + blockIdx.x*blockDim.x;

    unsigned mask = __ballot_sync(0xffffffff, tid<n);
    if (tid<n) {
        gpu::reduce_by_key(src[tid], dst, index[tid], mask);
    }
}

template <typename T>
std::vector<T> reduce(const std::vector<T>& in, size_t n_out, const std::vector<int>& index, unsigned block_dim=128) {
    EXPECT_EQ(in.size(), index.size());
    EXPECT_TRUE(std::is_sorted(index.begin(), index.end()));

    using array = memory::device_vector<T>;
    using iarray = memory::device_vector<int>;

    int n = in.size();

    array  src = memory::make_const_view(in);
    iarray idx = memory::make_const_view(index);
    array  dst(n_out, 0);

    unsigned grid_dim = (n-1)/block_dim + 1;
    reduce_kernel<<<grid_dim, block_dim>>>(src.data(), dst.data(), idx.data(), n);

    std::vector<T> out(n_out);
    memory::copy(dst, memory::make_view(out));

    return out;
}

TEST(reduce_by_key, no_repetitions)
{
    int n = 64;
    {
        std::vector<float> in(n, 1);
        std::vector<int> index = util::assign_from(util::make_span(0, n));

        auto out = reduce(in, n, index);
        for (auto o: out) EXPECT_EQ(o, 1.0f);
    }
    {
        std::vector<double> in(n, 1);
        std::vector<int> index = util::assign_from(util::make_span(0, n));

        auto out = reduce(in, n, index);
        for (auto o: out) EXPECT_EQ(o, 1.0);
    }
}

TEST(reduce_by_key, single_repeated_index)
{
    // Perform reduction of a sequence of 1s of length n
    // The expected result is n
    for (auto n: {1, 2, 7, 31, 32, 33, 63, 64, 65, 128}) {
        std::vector<double> in(n, 1);
        std::vector<int> index(n, 0);

        auto out = reduce(in, 1, index, 32);
        EXPECT_EQ(double(n), out[0]);
    }
    // Perform reduction of an ascending sequence of {1,2,3,...,n}
    // The expected result is n*(n+1)/2
    for (auto n: {1, 2, 7, 31, 32, 33, 63, 64, 65, 128}) {
        std::vector<double> in = util::assign_from(util::make_span(1, n+1));
        std::vector<int> index(n, 0);

        auto out = reduce(in, 1, index);
        EXPECT_EQ(out[0], double((n+1)*n/2));
    }
}

TEST(reduce_by_key, scatter)
{
    std::vector<int> index = {0,0,0,1,2,2,2,2,3,3,7,7,7,7,7,11};
    unsigned n = util::max_value(index)+1;
    std::vector<double> in(index.size(), 1);
    std::vector<double> expected = {3., 1., 4., 2., 0., 0., 0., 5., 0., 0., 0., 1.};

    unsigned m = index.size();

    EXPECT_EQ(n, expected.size());

    auto out = reduce(in, n, index);
    EXPECT_EQ(expected, out);

    // rerun with 7 threads per thread block, to test
    //  * using more than one thread block
    //  * thread blocks that are not a multiple of 32
    //  * thread blocks that are less than 32

    out = reduce(in, n, index, 7);
    EXPECT_EQ(expected, out);
}

template <typename I>
__global__
void key_set_pos_kernel(const I* index, int n, key_set_pos* out) {
    unsigned tid = threadIdx.x + blockIdx.x*blockDim.x;

    unsigned mask = __ballot_sync(0xffffffff, tid<n);
    key_set_pos k(index[tid], mask);
    out[tid] = k;
}

std::vector<key_set_pos> run_key_set_pos(const std::vector<int>& index, unsigned block_dim=128) {
    using memory::device_vector;
    using memory::make_view;

    unsigned n = index.size();
    device_vector<int> index_gpu(memory::make_const_view(index));
    device_vector<key_set_pos> ks_gpu(n);

    unsigned grid_dim = (n-1)/block_dim + 1;

    std::vector<key_set_pos> ks(n);
    memory::copy(ks_gpu, memory::make_view(ks));

    return ks;
}

TEST(reduce_by_key, key_set_pos) {
    using ivector = std::vector<int>;
    using uvector = std::vector<unsigned>;

    {
	ivector index = { 2, 2, 2, 2, 3, 4, 4 };
	auto ks = run_key_set_pos(index);

	uvector expected_width = { 4, 3, 2, 1, 1, 2, 1 };
	uvector expected_is_root = { 1, 0, 0, 0, 1, 1, 0 };
	uvector expected_lane_id = { 0, 1, 2, 3, 4, 5, 6 };
	uvector expected_key_mask(7, 0x7fu);

	for (unsigned i = 0; i<7; ++i) {
	    EXPECT_EQ(expected_width[i], ks[i].width);
	    EXPECT_EQ(expected_lane_id[i], ks[i].lane_id);
	    EXPECT_EQ(expected_key_mask[i], ks[i].key_mask);
	    EXPECT_EQ(expected_is_root[i], ks[i].is_root);
	}
    }
}

// Test kernels that perform more than one reduction in a single invokation.
// Used to reproduce and test for synchronization issues on V100 GPUs.

template <typename T, typename I>
__global__
void reduce_twice_kernel(const T* src, T* dst, const I* index, int n) {
    unsigned tid = threadIdx.x + blockIdx.x*blockDim.x;

    unsigned mask = __ballot_sync(0xffffffff, tid<n);
    if (tid<n) {
        gpu::reduce_by_key(src[tid], dst, index[tid], mask);
        gpu::reduce_by_key(src[tid], dst, index[tid], mask);
    }
}

template <typename T>
std::vector<T> reduce_twice(const std::vector<T>& in, size_t n_out, const std::vector<int>& index, unsigned block_dim=128) {
    EXPECT_EQ(in.size(), index.size());
    EXPECT_TRUE(std::is_sorted(index.begin(), index.end()));

    using array = memory::device_vector<T>;
    using iarray = memory::device_vector<int>;

    int n = in.size();

    array  src = memory::make_const_view(in);
    iarray idx = memory::make_const_view(index);
    array  dst(n_out, 0);

    unsigned grid_dim = (n-1)/block_dim + 1;
    reduce_twice_kernel<<<grid_dim, block_dim>>>(src.data(), dst.data(), idx.data(), n);

    std::vector<T> out(n_out);
    memory::copy(dst, memory::make_view(out));

    return out;
}

TEST(reduce_by_key, scatter_twice)
{
    std::vector<int> index = {0,0,0,1,2,2,3,7,7,7,11};
    unsigned n = util::max_value(index)+1;
    std::vector<double> in(index.size(), 1);
    std::vector<double> expected = {6., 2., 4., 2., 0., 0., 0., 6., 0., 0., 0., 2.};

    unsigned m = index.size();

    EXPECT_EQ(n, expected.size());

    auto out = reduce_twice(in, n, index);
    EXPECT_EQ(expected, out);

    // rerun with 7 threads per thread block, to test
    //  * using more than one thread block
    //  * thread blocks that are not a multiple of 32
    //  * thread blocks that are less than 32

    out = reduce_twice(in, n, index, 7);
    EXPECT_EQ(expected, out);
}
