
#define FV3  1
#define MPAS 2

#if !defined(FV3_STOCH) && !defined(MPAS_STOCH)
#error "One of the macros FV3_STOCH and MPAS_STOCH (but not both) must be defined."
#elif defined(FV3_STOCH) && defined(MPAS_STOCH)
#error "The macros FV3_STOCH and MPAS_STOCH cannot both be defined."
#elif defined(FV3_STOCH)
#define DYCORE FV3
#elif defined(MPAS_STOCH)
#define DYCORE MPAS
#endif

