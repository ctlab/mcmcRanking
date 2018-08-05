#ifndef MCMC_RANKING_HSA_H
#define MCMC_RANKING_HSA_H

#include <vector>
#include <unordered_map>

namespace mcmc {
    using namespace std;

    class HSA {
        unordered_map<unsigned, size_t> map;
        vector<unsigned> elements;
        vector<bool> contain;
    public:
        HSA(size_t size);

        void clear();

        size_t size();

        unsigned get(size_t ind);

        size_t get_index(unsigned v);

        vector<unsigned> get_all();

        bool contains(unsigned el);

        void insert(unsigned el);

        void erase(unsigned el);

        void swap(unsigned el, unsigned new_el);
    };
}

#endif //MCMC_RANKING_HSA_H
