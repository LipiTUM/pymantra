#include "LocalSearch.hpp"


template<typename T, typename F>
T opt(T a, T b, F comp) {
    if (comp(a, b)) return a;
    return b;
}


template<typename T, typename F>
T opt(T a, T b, T c, F comp) {
    if (comp(a, b)) {
        if (comp(a, c)) return a;
        return c;
    } else if (comp(b, c)) return b;
    return c;
}


template<typename T, typename F>
vector<T> LocalSearch::merge_solutions(
    vector<T>& ins, vector<T>& del, vector<T>& subs, F comp_fun
) {
    vector<T> merged(ins.size() + del.size() + subs.size());
    unsigned
    m = 0, i = 0, d = 0, s = 0;
    // do as long as no end is reached
    while ((i < ins.size()) & (d < del.size()) & (s < subs.size())) {
       merged[m] = opt(ins[i], del[d], subs[s], comp_fun);
       if (merged[m] == ins[i]) i++;
       else if (merged[m] == del[d]) d++;
       else s++;
       m++;
    }
    // all substitutions inserted
    while ((i < ins.size()) & (d < del.size())) {
        merged[m] = opt(ins[i], del[d], comp_fun);
        if (merged[m] == ins[i]) i++;
        else d++;
        m++;
    }
    // all deletions inserted
    while ((i < ins.size()) & (s < subs.size())) {
        merged[m] = opt(ins[i], subs[s], comp_fun);
        if (merged[m] == ins[i]) i++;
        else s++;
        m++;
    }
    // all insertions inserted
    while ((d < del.size()) & (s < subs.size())) {
        merged[m] = opt(del[d], subs[s], comp_fun);
        if (merged[m] == del[d]) d++;
        else s++;
        m++;
    }
    // only deletions left
    while (d < del.size()) {
        merged[m] = del[d];
        d++;
        m++;
    }
    // only substitutions left
    while (s < subs.size()) {
        merged[m] = subs[s];
        s++;
        m++;
    }
    // only insertions left
    while (i < ins.size()) {
        merged[m] = ins[i];
        i++;
        m++;
    }

    return merged;
}
