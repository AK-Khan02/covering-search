#define main cover_13_6_3_verifier_main
#include "cover_13_6_3_final_verifier.cpp"
#undef main

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: export_representatives CERT OUT\n";
        return 2;
    }
    std::string cert_path = argv[1];
    std::string out_path = argv[2];

    LocalGenerator gen9(9);
    auto rows = gen9.generate();
    if (rows.size() != 3385) {
        std::cerr << "unexpected row count " << rows.size() << "\n";
        return 1;
    }

    std::ifstream in(cert_path);
    if (!in) {
        std::cerr << "could not open certificate\n";
        return 1;
    }
    int n_classes;
    in >> n_classes;
    std::vector<int> reps(n_classes);
    for (int i = 0; i < n_classes; ++i) in >> reps[i];

    std::ofstream out(out_path);
    out << n_classes << "\n";
    for (int cls = 0; cls < n_classes; ++cls) {
        int rep = reps[cls];
        if (rep < 0 || rep >= (int)rows.size()) {
            std::cerr << "bad representative index\n";
            return 1;
        }
        out << cls << " " << rep;
        for (int mask : rows[rep]) out << " " << mask;
        out << "\n";
    }
    return 0;
}
