#include <ratio>
#include <string>
#include <chrono>
#include <system_error>

#include <time.h>

template <clockid_t clock>
struct posix_clock {
    using period = std::chrono::nanoseconds::period;
    using rep = std::chrono::nanoseconds::rep;
    using duration = std::chrono::duration<rep, period>;
    using time_point = std::chrono::time_point<posix_clock>;

    static time_point now() {
        timespec ts;
        int r = clock_gettime(clock, &ts);
        return r? throw std::system_error(errno, std::generic_category()): time_point(duration(ts.tv_sec*1000000000ll+ts.tv_nsec));
    }
};

using thread_cputime_clock = posix_clock<CLOCK_THREAD_CPUTIME_ID>;

void busy_wait_ns(double time_ns) {
    auto until = thread_cputime_clock::now()+std::chrono::nanoseconds((int)(time_ns));
    while (thread_cputime_clock::now()<until) ;
}

