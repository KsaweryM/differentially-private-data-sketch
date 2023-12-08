#include <iostream>
#include <iomanip>
#include <sstream>
#include <openssl/sha.h>
#include <cstring>
#include <assert.h>
#include <limits>
#include <cmath>
#include <random>

uint32_t get_random_number(uint32_t min, uint32_t max)
{
	static std::default_random_engine engine;
	std::uniform_int_distribution<uint32_t> distribution(min, max);

	return distribution(engine);
}

class Exp_sketch {
public:
    double estimate_cardinality(uint32_t sketch_size, uint32_t* multiset, uint32_t multiset_size)
    {
        double* sketch = new double[sketch_size];
        for (uint32_t i = 0; i < sketch_size; i++)
        {
            sketch[i] = std::numeric_limits<double>::infinity();
        }

        for (uint32_t i = 0; i < multiset_size; i++)
        {
            for (uint32_t j = 0; j < sketch_size; j++)
            {
                double m = static_cast<double>(multiset[i]);

                double U = get_uniform_sample(i, j);
                double E = -log(U) / m;

                assert(E >= 0.0);
                assert(sketch[j] >= 0.0);

                if (E < sketch[j])
                {
                    sketch[j] = E;
                }
            }
        }

        double estimation = estimator(sketch, sketch_size);

        delete[] sketch;

        return estimation;
    }

private:
    double estimator(double* sketch, uint32_t sketch_size)
    {
        double G_m = 0;

        for (uint32_t i = 0; i < sketch_size; i++)
        {
            G_m += sketch[i];
        }

        return static_cast<double>(sketch_size - 1) / G_m;
    }

    double get_uniform_sample(uint32_t i, uint32_t k)
    {
        const uint64_t source = (static_cast<uint64_t>(i) << 32) | static_cast<uint64_t>(k);

        unsigned char hash[SHA256_DIGEST_LENGTH];
        SHA256(reinterpret_cast<const unsigned char*>(&source), sizeof(source), hash);

        uint64_t uniform_value = 0;

        assert(sizeof(uniform_value) <= SHA256_DIGEST_LENGTH);
        std::memcpy(&uniform_value, hash, sizeof(uniform_value));

        const uint64_t max_value = -1;

        return static_cast<double>(uniform_value) / static_cast<double>(max_value);
    }
};


uint64_t sum_of_unique_elements(uint32_t* multiset, uint32_t size)
{
    uint64_t sum = 0;
    for (uint32_t i = 0; i < size; i++)
    {
        for (uint32_t j = 0; j < size; j++)
        {
            if (i == j)
            {
                continue;
            }

            if (multiset[i] == multiset[j])
            {
                break;
            }
        }

        sum += multiset[i];
    }

    return sum;
}

int main() {
    uint32_t multiset_size = 100000, min_value = 1, max_value = 100000, sketch_size = 100;
    // std::cin >> multiset_size >> min_value >> max_value >> sketch_size;

    uint32_t* multiset = new uint32_t[multiset_size];
    for (uint32_t i = 0; i < multiset_size; i++)
    {
        multiset[i] = get_random_number(min_value, max_value);
    }

    Exp_sketch exp_sketch;
    double estimation = exp_sketch.estimate_cardinality(sketch_size, multiset, multiset_size);
    std::cout << estimation << std::endl;

    double correct_sum = sum_of_unique_elements(multiset, multiset_size);
    std::cout << "correct_sum: " << correct_sum << ", div:" << abs(correct_sum - estimation)/correct_sum * 100.0 << "%" << std::endl;


    delete[] multiset;



    return 0;
}