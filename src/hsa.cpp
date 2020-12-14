#include "mcmc.h"

static const size_t NONE = -(size_t)1;

namespace mcmc {
    HSA::HSA(size_t size) : map(size, NONE), elements() {
    }

    void HSA::insert(size_t el) {
        map[el] = elements.size();
        elements.push_back(el);
    }

    void HSA::erase(size_t el) {
        size_t ind = map[el];
        if (ind == NONE) {
            throw std::invalid_argument("erasing non-existing element of a HSA");
        }
        elements[ind] = elements.back();
        map[elements[ind]] = ind;
        map[el] = NONE;
        elements.pop_back();
    }

    bool HSA::contains(size_t el) {
        return map[el] != NONE;
    }

    void HSA::swap(size_t el, size_t new_el) {
        map[new_el] = map[el];
        map[el] = NONE;
        elements[map[new_el]] = new_el;
    }

    void HSA::clear() {
        fill(map.begin(), map.end(), NONE);
        elements.clear();
    }

    size_t HSA::size() {
        return elements.size();
    }

    size_t HSA::get(size_t ind) {
        return elements[ind];
    }

    size_t HSA::get_index(size_t el) {
        return map[el];
    }

    vector<size_t> HSA::get_all() {
        return elements;
    }
};
