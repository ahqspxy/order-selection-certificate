#include <cfloat>
#include <cfenv>
#include <cstdint>
#include <cstring>
#include <limits>
static_assert(FLT_EVAL_METHOD == 0, "Binary64 operations must not use excess precision");
static_assert(std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits == 53);
extern "C" int certification_environment() {
    if (std::fegetround() != FE_TONEAREST) return 1;
    volatile double x = std::numeric_limits<double>::denorm_min() * 1024;
    volatile double two = 2.0;
    double y = x * two;
    std::uint64_t bits;
    std::memcpy(&bits, &y, sizeof(bits));
    if (bits != 2048) return 2;
    volatile double one = 1.0;
    volatile double half_ulp = 0x1p-53;
    double tie = one + half_ulp;
    if (tie != 1.0) return 3;
    return 0;
}
