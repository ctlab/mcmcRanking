#include "mcmc.h"

namespace mcmc {
    HSA::HSA(size_t size) : map(), elements(), contain(size, false) {
    }

    void HSA::clear() {
        map.clear();
        elements.clear();
        for (int i = 0; i < contain.size(); ++i) {
            contain[i] = false;
        }
    }

    size_t HSA::size() {
        return elements.size();
    }

    unsigned HSA::get(size_t ind) {
        return elements[ind];
    }

    size_t HSA::get_index(unsigned v) {
        return map.find(v)->second;
    }

    vector<unsigned> HSA::get_all() {
        return elements;
    }

    bool HSA::contains(unsigned el) {
        return contain[el];
    }

    void HSA::insert(unsigned el) {
        elements.push_back(el);
        map.insert(pair<unsigned, size_t>(el, elements.size() - 1));
        contain[el] = true;
    }

    void HSA::erase(unsigned el) {
        unordered_map<unsigned, size_t>::iterator it = map.find(el);
        size_t ind = it->second;
        map.erase(it);
        contain[el] = false;
        elements[ind] = elements[elements.size() - 1];
        elements.pop_back();
        if (ind != elements.size()) {
            it = map.find(elements[ind]);
            it->second = ind;
        }
    }

    void HSA::swap(unsigned el, unsigned new_el) {
        unordered_map<unsigned, size_t>::iterator it = map.find(el);
        size_t ind = it->second;
        map.erase(it);
        contain[el] = false;
        elements[ind] = new_el;
        map.insert(pair<unsigned, size_t>(new_el, ind));
        contain[new_el] = true;
    }
};
