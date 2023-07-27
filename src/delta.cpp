#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

double cardinality(const std::string& seq, int k) {
    // Implementation of the cardinality function
    std::unordered_set<std::string> substrings;
    for (int i = 0; i <= seq.length() - k; ++i) {
        substrings.insert(seq.substr(i, k));
    }
    return substrings.size();
}

double delta(const std::string& seq) {
    // Implementation of the delta function
    if (seq.empty()) {
        return 0;
    }
    
    double maxRatio = 0.0;
    for (int k = 1; k <= seq.length(); ++k) {
        double ratio = cardinality(seq, k) / k;
        maxRatio = std::max(maxRatio, ratio);
    }
    return maxRatio;
}

double delta_fast(const std::string& seq) {
    if (seq.empty()) {
        return 0;
    }

    double maxRatio = 0.0;
    std::unordered_set<std::string> substrings;

    for (int k = 1; k <= seq.length(); ++k) {
        substrings.clear();
        for (int i = 0; i <= seq.length() - k; ++i) {
            substrings.insert(seq.substr(i, k));
        }
        double ratio = static_cast<double>(substrings.size()) / k;
        maxRatio = std::max(maxRatio, ratio);
    }

    return maxRatio;
}

PYBIND11_MODULE(delta, m) {
    m.doc() = "Python binding for delta function";

    m.def("delta", &delta, "Calculate the delta value of a sequence");
    m.def("delta_fast", &delta_fast, "Calculate the delta value of a sequence");
}
