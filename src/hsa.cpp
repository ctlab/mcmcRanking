#include "mcmc.h"

namespace mcmc {
    HSA::HSA(size_t size) : map(size, -1), elements(), contain(size, false) {
    }

    void HSA::insert(size_t el) {
        map[el] = elements.size();
        elements.push_back(el);
        contain[el] = true;
    }

    void HSA::erase(size_t el) {
        if (!contain[el]) {
            throw std::invalid_argument("erasing non-existing element of a HSA");
        }
        size_t ind = map[el];
        elements[ind] = elements[elements.size() - 1];
        map[elements[ind]] = ind;
        map[el] = -1;
        contain[el] = false;
        elements.pop_back();
    }

    bool HSA::contains(size_t el) {
        return contain[el];
    }

    void HSA::swap(size_t el, size_t new_el) {
        map[new_el] = map[el];
        map[el] = -1;
        elements[map[new_el]] = new_el;
        contain[el] = false;
        contain[new_el] = true;
    }

    void HSA::clear() {
        fill(map.begin(), map.end(), -1);
        elements.clear();
        fill(contain.begin(), contain.end(), false);
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
