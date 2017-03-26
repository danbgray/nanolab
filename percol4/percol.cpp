#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <queue>
#include <map>
#include <random>
#include <set>
#include <stack>
#include <string>
#include <tuple>
#include <vector>
#include <utility>
#include <mkl_lapack.h>
#include <mkl_vsl.h>

#if 1
#define REAL float
#define XPPSV sppsv
#define vxRngUniform vsRngUniform
#else
#define REAL double
#define XPPSV dppsv
#define vxRngUniform vdRngUniform
#endif

#define Int MKL_INT
#define USE_MKL_UNIFORM_RNG 1
#define VLEFT 1.0
#define VRIGHT 0.0

typedef std::pair<Int,Int> RowCol;

std::string green(const std::string &s) {
    return std::string("\033[1;32m") + s + std::string("\033[0m");
}

std::string red(const std::string &s) {
    return std::string("\033[1;31m") + s + std::string("\033[0m");
}

std::string human_readable_time(double s) {
    char buf[100];
    if (std::abs(s) > 1) sprintf(buf, "%.4lgs", s);
    else if (std::abs(s) > 1e-3) sprintf(buf, "%.4lgms", s*1e3);
    else if (std::abs(s) > 1e-6) sprintf(buf, "%.4lgus", s*1e6);
    else if (std::abs(s) > 1e-9) sprintf(buf, "%.4lgns", s*1e9);
    else sprintf(buf,"%.3lgs", s);
    return std::string(buf);
}

class Grid {
    Int M, N;
    std::vector<Int> nodes;
    std::vector<REAL> right_edges;
    std::vector<REAL> down_edges;
    std::map<RowCol, REAL> vx; // voltage at node[at(r, c)]
    Int at(Int m, Int n) const { return m + M * n; }
 public:
    Grid(Int _M, Int _N) : M(_M), N(_N)
    {
        resize(M, N);
    }
    void resize(Int _M, Int _N) {
        M = _M;
        N = _N;
        nodes.resize(M*N);
        right_edges.resize(M*N); // edges for column N-1 are invalid
        down_edges.resize(M*N); // edges for column 0, N-1 and row M-1 are invalid
    }
    Int num_nodes() const { return M*N; }
    Int num_edges() const { return M*(N-1) + (M-1)*(N-2); }
    void print() {
        printf("Have %i nodes, %i edges.\n", nodes.size(), M*(N-1)+(M-1)*(N-2));
        for (Int r = 0; r < M; ++r) {
            for (Int c = 0; c < N; ++c) {
                Int n = at(r, c);
                double re = c < N - 1 ? right_edges[n] : 0;
                printf("%c", nodes[n] ? 'o' : ' ');
                printf("%s", re ? "---" : "   ");
            }
            printf("\n ");
            for (Int c = 1; r < M-1 && c < N-1; ++c) {
                Int n = at(r, c);
                double de = down_edges[n];
                printf("   %c", de ? '|' : ' ');
            }
            printf("\n");
        }
    }
    void random_10(REAL prob1, Int seed) {
#if USE_MKL_UNIFORM_RNG
        VSLStreamStatePtr stream;
        vslNewStream(&stream, VSL_BRNG_MT19937, seed);
        vxRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, right_edges.size(), right_edges.data(), 0, 1);
        vxRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, down_edges.size(), down_edges.data(), 0, 1);
        vslDeleteStream( &stream );
        for (auto& e : right_edges) e = e < prob1 ? 1 : 0;
        for (auto& e : down_edges) e = e < prob1 ? 1 : 0;
#else
        //std::default_random_engine generator;
        std::mt19937 generator;
        generator.seed(std::seed_seq(&seed, &seed + 1));
        std::uniform_real_distribution<REAL> distribution(0, 1);
        std::generate_n(right_edges.begin(), right_edges.size(), [&](){
            return distribution(generator) < prob1 ? 1 : 0; });
        std::generate_n(down_edges.begin(), down_edges.size(), [&](){
            return distribution(generator) < prob1 ? 1 : 0; });
#endif
        for (Int r = 0; r < M; ++r) right_edges[at(r, N-1)] = 0;
        for (Int r = 0; r < M; ++r) down_edges[at(r, 0)] = 0;
        for (Int r = 0; r < M; ++r) down_edges[at(r, N-1)] = 0;
        for (Int c = 0; c < N; ++c) down_edges[at(M-1, c)] = 0;
    }
    Int num_nonzero_edges() const {
        auto nonzero = [](Int e){ return !!e; };
        Int n
            = std::count_if( &right_edges[0], &right_edges[M*N], nonzero)
            + std::count_if( &down_edges[0], &down_edges[M*N], nonzero);
        return n;
    }
    Int num_nonzero_nodes() const {
        Int n = std::count_if(&nodes[0], &nodes[M*N], [](bool e){ return e; });
        return n;
    }
    double fraction_1() {
        return double(num_nonzero_edges()) / double(num_edges());
    }
    bool connected() {
        std::stack<RowCol> snodes;
        _zero_nodes();
        for (Int m = 0; m < M; ++m) {
            if (right_edges[m]) {
                snodes.push(RowCol(m,1));
            }
        }
        bool reached_right = _flood_fill(snodes, true);
        if (M * N <= 100) print();
        return reached_right;
    }
    Int all_connections() {
        Int nconnected = 0;
        std::stack<RowCol> snodes;
        _zero_nodes();
        auto conn = nodes;
        std::vector<char> skip(M);
        for (Int m = 0; m < M; ++m) {
            if (right_edges[m] && !skip[m]) {
                _zero_nodes();
                snodes.push(RowCol(m,1));
                bool reached_right = _flood_fill(snodes, false);
                if (reached_right) {
                    _update_from_nodes(conn);
                    nconnected += 1;
                }
                for (Int r = m; r < M; ++r) {
                    if (nodes[r + M]) skip[r] = 1;
                }
            }
        }
        nodes.swap(conn);
        if (M * N <= 100) print();
        return nconnected;
    }
    void trim_whiskers() {
        std::stack<RowCol> snodes;
        for (Int c = 1; c < N - 1; ++c) {
            for (Int r = 0; r < M; ++r) {
                Int v = at(r, c);
                if (nodes[v] && _nlinks(r, c) < 2) {
                    snodes.push(RowCol(r,c));
                }
            }
        }
        while (!snodes.empty()) {
            Int r = snodes.top().first, c = snodes.top().second;
            snodes.pop();
            Int v = at(r, c);
            nodes[v] = 0;
            if (r > 0 && nodes[v - 1] && _nlinks(r - 1, c) < 2) snodes.push(RowCol(r - 1, c));
            if (r < M - 1 && nodes[v + 1] && _nlinks(r + 1, c) < 2) snodes.push(RowCol(r + 1, c));
            if (c > 0 && nodes[v - M] && _nlinks(r, c - 1) < 2) snodes.push(RowCol(r, c - 1));
            if (c < N - 1 && nodes[v + M] && _nlinks(r, c + 1) < 2) snodes.push(RowCol(r, c + 1));
        }
    }
    Int solve_currents() {
        const Int BITVARNO = 0x80000000;
        std::vector<REAL> A, B;
        std::map<Int, RowCol> rc_xno; // variable number -> row, col
        // Assign variable numbers to nodes
        Int nx = 0, nd = 0;
        for (Int c = 0; c < N; ++c) {
            for (Int r = 0; r < M; ++r) {
                if (nodes[at(r, c)]) {
                    if (c == 0 || c == N - 1) {
                        nodes[at(r, c)] = nd++ | BITVARNO;
                    } else {
                        rc_xno.emplace(nx, RowCol(r, c));
                        nodes[at(r, c)] = nx++ | BITVARNO;
                    }
                }
            }
        }
        printf("solve_currents: nx=%i nd=%i\n", nx, nd);
        // Init matrix A and rhs B, for Ax=B
        A.resize(nx * (nx + 1) / 2);
        B.resize(nx);
        auto LPACK = [nx](Int r, Int c) { return r + c*nx - c*(c+1)/2; };
        for (Int c = 0; c < N - 1; ++c) {
            for (Int r = 0; r < M; ++r) {
                if (c + 1 < N
                    && right_edges[at(r, c)]
                    && nodes[at(r, c)] && nodes[at(r, c + 1)])
                {
                    REAL sigma = 1.0 / right_edges[at(r, c)];
                    if (c == 0) {
                        Int a = nodes[at(r, c + 1)] & ~BITVARNO;
                        REAL voltb = VLEFT;
                        A[LPACK(a, a)] += right_edges[at(r, c)] * sigma;
                        B[a] += right_edges[at(r, c)] * voltb * sigma;
                    } else if (c + 1 == N - 1) {
                        Int a = nodes[at(r, c)] & ~BITVARNO;
                        REAL voltb = VRIGHT;
                        A[LPACK(a, a)] += right_edges[at(r, c)] * sigma;
                        B[a] += right_edges[at(r, c)] * voltb * sigma;
                    } else {
                        Int a = nodes[at(r, c)] & ~BITVARNO;
                        Int b = nodes[at(r, c + 1)] & ~BITVARNO;
                        if (a < b) std::swap(a, b);
                        A[LPACK(a, a)] += right_edges[at(r, c)] * sigma;
                        A[LPACK(b, b)] += right_edges[at(r, c)] * sigma;
                        A[LPACK(a, b)] -= right_edges[at(r, c)] * sigma;
                    }
                }
                if (r < M - 1
                    && down_edges[at(r, c)]
                    && nodes[at(r, c)] && nodes[at(r + 1, c)])
                {
                    REAL sigma = 1.0 / down_edges[at(r, c)];
                    Int a = nodes[at(r, c)] & ~BITVARNO;
                    Int b = nodes[at(r + 1, c)] & ~BITVARNO;
                    if (a < b) std::swap(a, b);
                    A[LPACK(a, a)] += down_edges[at(r, c)] * sigma;
                    A[LPACK(b, b)] += down_edges[at(r, c)] * sigma;
                    A[LPACK(a, b)] -= down_edges[at(r, c)] * sigma;
                }
            }
        }
        // Solve the system
        if (1) {
            Int ONERHS = 1;
            Int INFO = 0;
            XPPSV("L", &nx, &ONERHS, A.data(), B.data(), &nx, &INFO );
            if (INFO) {
                printf("INFO=%i\n", INFO);
                return INFO;
            }
            vx.clear();
            for (Int i = 0; i < nx; ++i) {
                vx.emplace(rc_xno[i], B[i]);
            }
        }
        return 0;
    }
    void write_currents(const std::string& filename) {
        std::ofstream out(filename);
        for (Int c = 0; c < N - 1; ++c) {
            for (Int r = 0; r < M; ++r) {
                if (nodes[at(r, c)] && nodes[at(r, c + 1)] && right_edges[at(r, c)]) {
                    REAL sigma = 1.0 / right_edges[at(r, c)];
                    REAL v0 = c == 0 ? VLEFT : vx.at({r, c});
                    REAL v1 = c + 1 == N - 1 ? VRIGHT : vx.at({r, c + 1});
                    REAL deltav = v0 - v1;
                    REAL current = std::abs(deltav * sigma);
                    out << c << " " << M - 1 - r << " " << c + 1 << " " << M - 1 - r << " " << current << std::endl;
                }
                if (c > 0 && nodes[at(r, c)] && nodes[at(r + 1, c)] && down_edges[at(r, c)]) {
                    REAL sigma = 1.0 / down_edges[at(r, c)];
                    REAL v0 = vx.at({r, c});
                    REAL v1 = vx.at({r + 1, c});
                    REAL deltav = v0 - v1;
                    REAL current = std::abs(deltav * sigma);
                    out << c << " " << M - 1 - r << " " << c << " " << M - 1 - (r + 1) << " " << current << std::endl;
                }
            }
        }
    }
    REAL left_current() {
        REAL res = 0;
        for (Int r = 0; r < M; ++r) {
            if (nodes[at(r, 0)] && nodes[at(r, 1)] && right_edges[at(r, 0)]) {
                REAL sigma = 1.0 / right_edges[at(r, 0)];
                REAL deltav = VLEFT - vx.at({r, 1});
                res += deltav * sigma;
            }
        }
        return res;
    }
    REAL right_current() {
        REAL res = 0;
        for (Int r = 0; r < M; ++r) {
            if (nodes[at(r, N-2)] && nodes[at(r, N-1)] && right_edges[at(r, N-2)]) {
                REAL sigma = 1.0 / right_edges[at(r, N - 2)];
                REAL deltav = vx.at({r, N-2}) - VRIGHT;
                res += deltav * sigma;
            }
        }
        return res;
    }
private:
    Int _nlinks(Int r, Int c) const {
        Int v = at(r, c);
        Int nu = r > 0 && down_edges[v - 1] && nodes[v - 1];
        Int nl = c > 0 && right_edges[v - M] && nodes[v - M];
        Int nr = c < N - 1 && right_edges[v] && nodes[v + M];
        Int nd = r < M - 1 && down_edges[v] && nodes[v + 1];
        return nu + nl + nr + nd;
    }
    void _update_from_nodes(decltype(nodes)& conn) {
        for (Int c = 0; c < N; ++c) {
            for (Int r = 0; r < M; ++r) {
                if (nodes[at(r, c)]) {
                    conn[at(r, c)] = 1;
                }
            }
        }
    }
    void _zero_nodes() {
        std::fill(nodes.begin(), nodes.end(), 0);
    }
    /// 'rrr' = return on reaching right side
    bool _flood_fill(std::stack<RowCol>& snodes, bool rrr)
    {
        bool rr = false;
        while (!snodes.empty()) {
            Int m = snodes.top().first, n = snodes.top().second;
            snodes.pop();
            if (n == N - 1) rr = true;
            if (rrr && rr) break;
            Int c = at(m, n);
            if (n > 0 && right_edges[c - M] && nodes[c - M] == 0) {
                nodes[c - M] = 1;
                snodes.push(RowCol(m, n - 1));
            }
            if (m > 0 && down_edges[c - 1] && nodes[c - 1] == 0) {
                nodes[c - 1] = 1;
                snodes.push(RowCol(m - 1, n));
            }
            if (m < M - 1 && down_edges[c] && nodes[c + 1] == 0) {
                nodes[c + 1] = 1;
                snodes.push(RowCol(m + 1, n));
            }
            if (n < N - 1 && right_edges[c] && nodes[c + M] == 0) {
                nodes[c + M] = 1;
                snodes.push(RowCol(m, n + 1));
            }
        }
        return rr;
    }
};

template<typename T>
double seconds(const T&& t) {
    namespace c = std::chrono;
    return 1e-9*c::duration_cast<c::nanoseconds>(t).count();
}

void dwim1(Int M, Int N, Int seed, REAL frac1) {

    using clock = std::chrono::high_resolution_clock;

    Grid grid(M, N);

    auto t0 = clock::now();
    grid.random_10(frac1, seed);
    auto t1 = clock::now();
    double dt_gen = seconds(t1 - t0);

    double fill = grid.fraction_1();

    t0 = clock::now();
    //bool con = grid.connected();
    Int con = grid.all_connections();
    t1 = clock::now();
    double dt_con = seconds(t1 - t0);

    std::string res = con ? green("connected") : red("broken");
    printf("%i-by-%i graph fill %lg%+.5lg (%i/%i) is %s (%i)\n",
        M, N, frac1, fill - frac1,
        grid.num_nonzero_edges(), grid.num_edges(),
        res.c_str(), con);
    printf("rng_time:%s/item(%s), connected_time:%s/item(%s)\n",
        human_readable_time(dt_gen/(M*N)).c_str(), human_readable_time(dt_gen).c_str(),
        human_readable_time(dt_con/(M*N)).c_str(), human_readable_time(dt_con).c_str());

    t0 = clock::now();
    Int n0 = grid.num_nonzero_nodes();
    t1 = clock::now();
    printf("nznodes %i/%i = %lg (%s) in connecting clusters\n",
        n0, M*N, double(n0)/(M*N),
        human_readable_time(seconds(t1 - t0)).c_str());

    t0 = clock::now();
    grid.trim_whiskers();
    t1 = clock::now();
    Int n1 = grid.num_nonzero_nodes();

    printf("nznodes %i/%i = %lg (%s) after trim_whiskers\n",
        n1, M*N, double(n1)/(M*N),
        human_readable_time(seconds(t1 - t0)).c_str());
    if (M*N <= 100) grid.print();

    if (con) {
        t0 = clock::now();
        Int failed = grid.solve_currents();
        t1 = clock::now();
        REAL lc = grid.left_current();
        REAL rc = grid.right_current();
        printf("currents = %lg|%lg (rel.defect=%lg) (%s)\n",
            lc, rc, (lc - rc) / std::max(std::abs(lc), std::abs(rc)),
            human_readable_time(seconds(t1 - t0)).c_str());
        if (!failed) grid.write_currents("cur.dat");
    }
}


int main(int argc, char *argv[]) {

    Int M = 50; //3000; // rows
    Int N = 50; //3000; // cols

    Int seed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    //seed = 3857272299;
    printf("Seed = %u\n", seed);
    dwim1(M, N, seed, 0.5);
    printf("Seed = %u\n", seed);

    return 0;
}
