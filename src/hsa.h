#ifndef MCMC_RANKING_HSA_H
#define MCMC_RANKING_HSA_H

#include <vector>

namespace mcmc {
    using namespace std;

    class HSA {
        vector<size_t> map;
        vector<size_t> elements;
        vector<bool> contain;
    public:
        HSA(size_t size);

        void insert(size_t el);

        void erase(size_t el);

        bool contains(size_t el);

        void swap(size_t el, size_t new_el);

        void clear();

        size_t size();

        size_t get(size_t ind);

        size_t get_index(size_t el);

        vector<size_t> get_all();
    };
}

#endif //MCMC_RANKING_HSA_H
