// https://gist.github.com/miki16g/8319640

#include <array>
#include <boost/generator_iterator.hpp>
#include <boost/random.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <cstdlib>
#include <iostream>
#include <string>

template <class T, std::size_t N>
inline T sum(std::array<T, N> &ar)
{
    T s = 0;
    for (T &e : ar)
        s += e;
    return s;
}

template <std::size_t N>
std::array<int, N> multinomial(int n, std::array<float, N> &pvals)
{
    static unsigned int seed = 0;
    boost::mt19937 rng;
    rng.seed(++seed);
    std::array<int, N> sample;
    sample.fill(0);
    try
    {
        if (sum(pvals) - pvals.at(N - 1) > 1.0)
        {
            throw "summation of pvals[:-1] >= 1.0";
        }
        auto remp = 1.0f;
        int rem = n, idx = 0;
        for (auto &pval : pvals)
        {
            auto p = pval / remp;
            boost::binomial_distribution<> binom(rem, p);
            boost::variate_generator<boost::mt19937 &, boost::binomial_distribution<>> bpdf(rng, binom);
            int hit = bpdf();
            if (hit > 0)
            {
                sample.at(idx) = hit;
            }
            rem -= hit;
            remp -= pval;
            idx += 1;
        }
    }
    catch (std::string errmsg)
    {
        std::cerr << errmsg << std::endl;
        exit(1);
    }
    return sample;
}
