#include <iostream>
#include <xpas/phylo_kmer_db.h>

xpas::phylo_kmer_db create_db()
{
    const size_t kmer_size = 3;
    const xpas::phylo_kmer::score_type omega = 1.0;
    const std::string tree;

    xpas::phylo_kmer_db db { kmer_size, omega, tree };

    /// branch 0
    db.unsafe_insert(0, { 0, 0.00f });
    db.unsafe_insert(1, { 0, 0.10f });
    db.unsafe_insert(2, { 0, 0.20f });

    /// branch 1
    db.unsafe_insert(1, { 1, 0.11f });
    db.unsafe_insert(2, { 1, 0.21f });
    db.unsafe_insert(3, { 1, 0.31f });

    /// branch 2
    db.unsafe_insert(2, { 2, 0.22f });
    db.unsafe_insert(3, { 2, 0.32f });
    db.unsafe_insert(4, { 2, 0.42f });

    return db;
}

void search(const xpas::phylo_kmer_db& db, xpas::phylo_kmer_db::key_type key)
{
    if (auto entries = db.search(key); entries)
    {
        std::cout << "Found " << key << ":\n";
        for (const auto& [branch, score] : *entries)
        {
            std::cout << "\tbranch " << branch << ": " << score << '\n';
        }
    }
    else
    {
        std::cout << "Key " << key << " not found.\n";
    }
}

int main()
{
    const auto db = create_db();
    search(db, 0);
    search(db, 2);
    search(db, 42);
}