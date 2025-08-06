#ifndef _EXPR_HASHER_HPP
#define _EXPR_HASHER_HPP

namespace std {

    template<>
    struct hash<z3::expr> {
        std::size_t operator()(const z3::expr &k) const {
            return k.hash();
        }
    };

    // Do not use Z3's == operator in the hash table
    template<>
    struct equal_to<z3::expr> {
        bool operator()(const z3::expr &lhs, const z3::expr &rhs) const {
            return z3::eq(lhs, rhs);
        }
    };
}

#endif
