#include <iostream>
#include <iomanip>
#include <sstream>
#include <openssl/sha.h>
#include <cstring>
#include <assert.h>
#include <limits>
#include <cmath>
#include <random>
#include <utility>

uint32_t get_random_number(uint32_t min, uint32_t max)
{
	static std::default_random_engine engine;
	std::uniform_int_distribution<uint32_t> distribution(min, max);

	return distribution(engine);
}

class Exp_sketch {
public:
    Exp_sketch(uint32_t sketch_size) : sketch(std::vector<double>(sketch_size, std::numeric_limits<double>::infinity()))
    {
        std::cout << "Basic ExpSketch" << std::endl;
    }

    void update_sketch(const std::vector<std::pair<uint32_t, uint32_t>>& multiset)
    {
        for (const auto& element: multiset)
        {
            uint32_t identifier = element.first;
            uint32_t attribute = element.second;

            for (uint32_t j = 0; j < sketch.size(); j++)
            {
                double m = static_cast<double>(attribute);

                double U = get_uniform_sample(identifier, j);
                double E = -log(U) / m;

                assert(E >= 0.0);
                assert(sketch[j] >= 0.0);

                if (E < sketch[j])
                {
                    sketch[j] = E;
                }
            }
        }
    }

    virtual double estimate()
    {
        double G_m = 0;

        for (uint32_t i = 0; i < sketch.size(); i++)
        {
            G_m += sketch[i];
        }

        return static_cast<double>(sketch.size() - 1) / G_m;
    }

protected:
    std::vector<double> sketch;

private:
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

class DF_exp_sketch : public Exp_sketch {
public:
    DF_exp_sketch(double epsilon, uint32_t max_attribute, uint32_t sketch_size) :
        Exp_sketch(sketch_size),
        _large_attribute(get_large_attribute(epsilon, max_attribute)),
        _max_exp_sketch_value(get_max_exp_sketch_value(epsilon, max_attribute))
    {
        std::cout << "DF ExpSketch with m * epsilon: " << epsilon * static_cast<double>(sketch_size) << std::endl;

        const uint32_t unique_identifier = std::numeric_limits<uint32_t>::max();
        std::vector<std::pair<uint32_t, uint32_t>> multiset = {std::make_pair(unique_identifier, _large_attribute)};

        update_sketch(multiset);
    }

    double estimate() override
    {
        double G_m = 0;

        for (uint32_t i = 0; i < sketch.size(); i++)
        {
            G_m += std::min(sketch[i], _max_exp_sketch_value);
        }

        return (static_cast<double>(sketch.size() - 1) / G_m) - static_cast<double>(_large_attribute);
    }

private:
    const uint32_t _large_attribute;
    const double _max_exp_sketch_value;

    uint32_t get_large_attribute(double epsilon, uint32_t max_attribute)
    {
        double e_power_y_minus_one = exp(epsilon) - 1.0;

        if (e_power_y_minus_one < 0.0)
        {
             throw std::invalid_argument("Epsilon is too small.");
        }

        double attribute = static_cast<double>(max_attribute) / e_power_y_minus_one;

        if (std::isinf(attribute))
        {
            throw std::invalid_argument("Epsilon is too small or Max attribute is too large.");
        }

        if (attribute > std::numeric_limits<uint32_t>::max())
        {
            throw std::invalid_argument("Generated large attribute is too large for uint32_t.");
        }

        return static_cast<uint32_t>(attribute);
    }

    double get_max_exp_sketch_value(double epsilon, uint32_t max_attribute)
    {
        return epsilon / static_cast<double>(max_attribute);
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

    std::vector<std::pair<uint32_t, uint32_t>> multiset;

    const double epsilon = 0.04;
    const uint32_t sketch_size = 10;
    const uint32_t nr_unique_elements = 10000;
    const uint32_t max_attribute = 100;
    const uint32_t max_nr_repetitions = 1;

    uint64_t sum_of_unique_elements = 0;

    std::cout << "The data has been initialized." << std::endl;
    for (uint32_t i = 0; i < nr_unique_elements; i++)
    {
        uint32_t argument = get_random_number(1, max_attribute);

        for (uint32_t j = 0; j < max_nr_repetitions; j++)
        {
            multiset.push_back(std::make_pair(i, argument));
        }

        sum_of_unique_elements += argument;
    }

    std::cout << "The multiset has been created." << std::endl;

    std::vector<Exp_sketch> sketches = {DF_exp_sketch(epsilon, max_attribute, sketch_size), Exp_sketch(sketch_size)};

    for (auto& sketch : sketches)
    {
        sketch.update_sketch(multiset);
        double estimated_sum = sketch.estimate();
        double exact_sum = static_cast<double>(sum_of_unique_elements);


        std::cout << "exact_sum: " << exact_sum << std::endl;
        std::cout << "Sum estimated by sketch: " << estimated_sum<< ", div:" << fabs(estimated_sum - exact_sum)/exact_sum * 100.0 << "%" << std::endl;

    }
    return 0;
}