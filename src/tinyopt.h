#pragma once

#include <cstring>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <type_traits>
#include <vector>

namespace to {

constexpr struct nothing_t {} nothing;

template <typename T>
struct maybe {
    bool ok = false;
    alignas(T) char data[sizeof(T)];

    maybe(): ok(false) {}
    maybe(nothing_t): ok(false) {}
    maybe(const T& v): ok(true) { construct(v); }
    maybe(T&& v): ok(true) { construct(std::move(v)); }
    template <typename U>
    maybe(const maybe<U>& m): ok(m.ok) { if (ok) construct(*m); }
    template <typename U>
    maybe(maybe<U>&& m): ok(m.ok) { if (ok) construct(std::move(*m)); }

    maybe& operator=(nothing_t) { destroy(); return *this; }
    maybe& operator=(const T& v) { assign(v); return *this; }
    maybe& operator=(T&& v) { assign(std::move(v)); return *this; }
    maybe& operator=(const maybe& m) { m.ok? assign(*m): destroy(); return *this; }
    maybe& operator=(maybe&& m) { m.ok? assign(std::move(*m)): destroy(); return *this; }

    const T& value() const & { assert_ok(); return *vptr(); }
    T&& value() && { assert_ok(); return std::move(*vptr()); }

    const T& operator*() const & { return *vptr(); }
    T&& operator*() && { return std::move(*vptr()); }
    explicit operator bool() const { return ok; }

private:
    T* vptr() { return reinterpret_cast<T*>(data); }
    const T* vptr() const { return reinterpret_cast<const T*>(data); }

    void construct(const T& v) { new (data) T(v); ok = true; }
    void construct(T&& v) { new (data) T(std::move(v)); ok = true; }

    void assign(const T& v) { ok? *vptr()=v: construct(v); }
    void assign(T&& v) { ok? *vptr()=std::move(v): construct(std::move(v)); }

    void destroy() { if (ok) (**this).~T(); ok = false; }

    void assert_ok() const { if (!ok) throw std::invalid_argument("is nothing"); }
};

template <>
struct maybe<void> {
    bool ok = false;

    maybe(): ok(false) {}
    template <typename U>
    maybe(const maybe<U>& m): ok(m.ok) {}
    template <typename U>
    maybe(maybe<U>&& m): ok(m.ok) {}
    maybe(nothing_t): ok(false) {}
    template <typename X> maybe(X&&): ok(true) {}

    maybe& operator=(nothing_t) { ok = false; return *this; }
    maybe& operator=(const maybe& m) { ok = m.ok; return *this; }
    template <typename U>
    maybe& operator=(U&& v) { ok = true; return *this; }

    explicit operator bool() const { return ok; }
};

template <typename F, typename T>
auto operator<<(F&& f, const maybe<T>& m) -> maybe<std::decay_t<decltype(f(*m))>> {
    if (m) return f(*m); else return nothing;
}

template <typename F>
auto operator<<(F&& f, const maybe<void>& m) -> maybe<std::decay_t<decltype(f())>> {
    if (m) return f(), maybe<void>(1); else return nothing;
}

template <typename T, typename U>
auto operator<<(T& x, const maybe<U>& m) -> maybe<std::decay_t<decltype(x=*m)>> {
    if (m) return x=*m; else return nothing;
}

template <typename T>
std::ostream& operator<<(std::ostream& x, const maybe<T>& m) {
    if (m) return x << *m; else return x << "()";
}

struct parse_opt_error: public std::runtime_error {
    parse_opt_error(const std::string& s): std::runtime_error(s) {}
    parse_opt_error(const char *arg, const std::string& s):
        std::runtime_error(s+": "+arg) {}
};

void usage(const char* argv0, const std::string& usage_str) {
    const char* basename = std::strrchr(argv0, '/');
    basename = basename? basename+1: argv0;

    std::cout << "Usage: " << basename << " " << usage_str << "\n";
}

void usage(const char* argv0, const std::string& usage_str, const std::string& parse_err) {
    const char* basename = std::strrchr(argv0, '/');
    basename = basename? basename+1: argv0;

    std::cerr << basename << ": " << parse_err << "\n";
    std::cerr << "Usage: " << basename << " " << usage_str << "\n";
}

template <typename V>
struct default_parser {
    maybe<V> operator()(const std::string& text) const {
        V v;
        std::istringstream stream(text);
        stream >> v;
        return stream? maybe<V>(v): nothing;
    }
};

template <typename V>
class keyword_parser {
    std::vector<std::pair<std::string, V>> map_;

public:
    template <typename KeywordPairs>
    keyword_parser(const KeywordPairs& pairs): map_(std::begin(pairs), std::end(pairs)) {}

    maybe<V> operator()(const std::string& text) const {
        for (const auto& p: map_) {
            if (text==p.first) return p.second;
        }
        return nothing;
    }
};

template <typename V>
class delimited_parser {
    std::string delim_;

public:
    explicit delimited_parser(const std::string& delim): delim_(delim) {}
    maybe<std::vector<V>> operator()(const std::string& text) const {
        std::vector<V> values;
        std::size_t p = 0, n = text.size();

        while (p<n) {
            std::size_t q = text.find(delim_, p);
            if (q==std::string::npos) q = n;

            V v;
            std::istringstream stream(text.substr(p, q-p));
            stream >> v;
            if (!stream) return nothing;
            else values.push_back(std::move(v));

            p = q+delim_.size();
        }
        return values;
    }
};

template <typename KeywordPairs>
auto keywords(const KeywordPairs& pairs) -> keyword_parser<decltype(std::begin(pairs)->second)> {
    return keyword_parser<decltype(std::begin(pairs)->second)>(pairs);
}

template <typename V = std::string, typename P = default_parser<V>, typename = std::enable_if_t<!std::is_same<V, void>::value>>
maybe<V> parse_opt(char **& argp, char shortopt, const char* longopt=nullptr, const P& parse = P{}) {
    const char* arg = argp[0];

    if (!arg || arg[0]!='-') {
        return nothing;
    }

    std::string text;

    if (arg[1]=='-' && longopt) {
        const char* rest = arg+2;
        const char* eq = std::strchr(rest, '=');

        if (!std::strcmp(rest, longopt)) {
            if (!argp[1]) throw parse_opt_error(arg, "missing argument");
            text = argp[1];
            argp += 2;
        }
        else if (eq && !std::strncmp(rest, longopt,  eq-rest)) {
            text = eq+1;
            argp += 1;
        }
        else {
            return nothing;
        }
    }
    else if (shortopt && arg[1]==shortopt && arg[2]==0) {
        if (!argp[1]) throw parse_opt_error(arg, "missing argument");
        text = argp[1];
        argp += 2;
    }
    else {
        return nothing;
    }

    auto v = parse(text);

    if (!v) throw parse_opt_error(arg, "failed to parse option argument");
    return v;
}

maybe<void> parse_opt(char **& argp, char shortopt, const char* longopt) {
    if (!*argp || *argp[0]!='-') {
        return nothing;
    }
    else if (argp[0][1]=='-' && longopt && !std::strcmp(argp[0]+2, longopt)) {
        ++argp;
        return true;
    }
    else if (shortopt && argp[0][1]==shortopt && argp[0][2]==0) {
        ++argp;
        return true;
    }
    else {
        return nothing;
    }
}

} // namespace to
