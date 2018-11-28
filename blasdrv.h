#include <string>
#include <dlfcn.h>

#define CHECK_HANDLE(h)                                     \
    if ((h) == NULL) {                                      \
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__,  \
                dlerror());                                 \
    }

#define SYMLINK_SO(f)                  \
    if (access(f "_so", F_OK) != -1 && \
        access(f ".so", F_OK) == -1) { \
        symlink(f "_so", f ".so");     \
    }

namespace {

    class cblas_mkl_t {
    public:
        typedef int MKL_INT;
        enum CBLAS_LAYOUT {
            CblasRowMajor  = 101,
            CblasColMajor  = 102
        };
        enum CBLAS_TRANSPOSE {
            CblasNoTrans   = 111,
            CblasTrans     = 112,
            CblasConjTrans = 113
        };
        enum CBLAS_UPLO {
            CblasUpper     = 121,
            CblasLower     = 122
        };
        enum CBLAS_SIDE {
            CblasLeft      = 141,
            CblasRight     = 142
        };
        void (*_get_version_string)(char *, int);
        void (*_set_num_threads_local)(int);
        void (*_set_num_threads)(int);
        int (*_get_max_threads)(void);
        void (*_set_num_stripes)(int);
        int (*_get_num_stripes)(void);
        void (*_domain_set_num_threads)(int);
        int (*_domain_get_max_threads)(void);
        void (*_set_dynamic)(int);
        int (*_get_dynamic)(void);
        float (*_sdot)(const MKL_INT, const float *, const MKL_INT,
                       const float *, const MKL_INT);
        void (*_scopy)(const MKL_INT, const float *, const MKL_INT,
                       float *, const MKL_INT);
        void (*_sscal)(const MKL_INT, const float, float *,
                       const MKL_INT);
        void (*_sgemv)(const CBLAS_LAYOUT, const CBLAS_TRANSPOSE,
                       const MKL_INT, const MKL_INT, const float,
                       const float *, const MKL_INT, const float *,
                       const MKL_INT, const float, float *,
                       const MKL_INT);
        void (*_ssymv)(const CBLAS_LAYOUT, const CBLAS_UPLO,
                       const MKL_INT, const float, const float *,
                       const MKL_INT, const float *, const MKL_INT,
                       const float, float *, const MKL_INT);
        void (*_sger)(const CBLAS_LAYOUT, const MKL_INT,
                      const MKL_INT, const float, const float *,
                      const MKL_INT, const float *, const MKL_INT,
                      float *, const MKL_INT);
        void (*_sgemm)(const CBLAS_LAYOUT, const CBLAS_TRANSPOSE,
                       const CBLAS_TRANSPOSE, const MKL_INT,
                       const MKL_INT, const MKL_INT, const float,
                       const float *, const MKL_INT, const float *,
                       const MKL_INT, const float, float *,
                       const MKL_INT);
        void (*_ssymm)(const CBLAS_LAYOUT, const CBLAS_SIDE,
                       const CBLAS_UPLO, const MKL_INT,
                       const MKL_INT, const float, const float *,
                       const MKL_INT, const float *, const MKL_INT,
                       const float, float *, const MKL_INT);
        void (*_vsmul)(const MKL_INT, const float [], const float [],
                       float []);
        cblas_mkl_t(void)
            : _get_version_string(NULL),
              _set_num_threads_local(NULL), _set_num_threads(NULL),
              _get_max_threads(NULL), _set_num_stripes(NULL),
              _get_num_stripes(NULL), _domain_set_num_threads(NULL),
              _domain_get_max_threads(NULL),
              _set_dynamic(NULL), _get_dynamic(NULL),
              _sdot(NULL), _scopy(NULL), _sscal(NULL), _sgemv(NULL),
              _ssymv(NULL), _sger(NULL), _sgemm(NULL), _ssymm(NULL),
              _vsmul(NULL)
        {
            SYMLINK_SO("libmkl_avx");
            SYMLINK_SO("libmkl_avx2");
            SYMLINK_SO("libmkl_def");
            SYMLINK_SO("libmkl_vml_avx");
            SYMLINK_SO("libmkl_vml_avx2");
            SYMLINK_SO("libmkl_vml_def");

            // The (non-"single") dynamic loading method for Intel
            // Math Kernel Library
            SYMLINK_SO("libiomp5");
            CHECK_HANDLE(dlopen("./libiomp5_so",
                                RTLD_NOW | RTLD_GLOBAL));
            // libmkl_core.so has circular dependencies and has to be
            // loaded lazily
            SYMLINK_SO("libmkl_core");
            CHECK_HANDLE(dlopen("./libmkl_core_so",
                                RTLD_LAZY | RTLD_GLOBAL));
            SYMLINK_SO("libmkl_intel_thread");
            CHECK_HANDLE(dlopen("./libmkl_intel_thread_so",
                                RTLD_NOW | RTLD_GLOBAL));
            SYMLINK_SO("libmkl_core");
            CHECK_HANDLE(dlopen("./libmkl_core_so",
                                RTLD_NOW | RTLD_GLOBAL));

            SYMLINK_SO("libmkl_intel_lp64");
            void *mkl_intel_lp64 = dlopen("./libmkl_intel_lp64_so",
                                          RTLD_NOW | RTLD_GLOBAL);

            CHECK_HANDLE(mkl_intel_lp64);

            _get_version_string =
                reinterpret_cast<void (*)(char *, int)>
                (dlsym(mkl_intel_lp64, "MKL_Get_Version_String"));
            CHECK_HANDLE(_get_version_string);
            _set_num_threads_local =
                reinterpret_cast<void (*)(int)>
                (dlsym(mkl_intel_lp64, "MKL_Set_Num_Threads_Local"));
            CHECK_HANDLE(_set_num_threads_local);
            _set_num_threads =
                reinterpret_cast<void (*)(int)>
                (dlsym(mkl_intel_lp64, "MKL_Set_Num_Threads"));
            CHECK_HANDLE(_set_num_threads);
            _get_max_threads =
                reinterpret_cast<int (*)(void)>
                (dlsym(mkl_intel_lp64, "MKL_Get_Max_Threads"));
            CHECK_HANDLE(_get_max_threads);
            _set_num_stripes =
                reinterpret_cast<void (*)(int)>
                (dlsym(mkl_intel_lp64, "MKL_Set_Num_Stripes"));
            CHECK_HANDLE(_set_num_stripes);
            _get_num_stripes =
                reinterpret_cast<int (*)(void)>
                (dlsym(mkl_intel_lp64, "MKL_Get_Num_Stripes"));
            CHECK_HANDLE(_get_num_stripes);
            _domain_set_num_threads =
                reinterpret_cast<void (*)(int)>
                (dlsym(mkl_intel_lp64, "MKL_Domain_Set_Num_Threads"));
            CHECK_HANDLE(_set_num_threads);
            _domain_get_max_threads =
                reinterpret_cast<int (*)(void)>
                (dlsym(mkl_intel_lp64, "MKL_Domain_Get_Max_Threads"));
            CHECK_HANDLE(_domain_get_max_threads);
            _set_dynamic =
                reinterpret_cast<void (*)(int)>
                (dlsym(mkl_intel_lp64, "MKL_Set_Dynamic"));
            CHECK_HANDLE(_set_dynamic);
            _get_dynamic =
                reinterpret_cast<int (*)(void)>
                (dlsym(mkl_intel_lp64, "MKL_Get_Dynamic"));
            CHECK_HANDLE(_get_dynamic);

            _sdot =
                reinterpret_cast<
                float (*)(const MKL_INT, const float *,
                          const MKL_INT, const float *,
                          const MKL_INT)>
                (dlsym(mkl_intel_lp64, "cblas_sdot"));
            CHECK_HANDLE(_sdot);
            _scopy =
                reinterpret_cast<
                void (*)(const MKL_INT, const float *,
                         const MKL_INT, float *, const MKL_INT)>
                (dlsym(mkl_intel_lp64, "cblas_scopy"));
            CHECK_HANDLE(_scopy);
            _sscal =
                reinterpret_cast<
                void (*)(const MKL_INT, const float, float *,
                         const MKL_INT)>
                (dlsym(mkl_intel_lp64, "cblas_sscal"));
            CHECK_HANDLE(_sscal);
            _sgemv =
                reinterpret_cast<
                void (*)(const CBLAS_LAYOUT, const CBLAS_TRANSPOSE,
                         const MKL_INT, const MKL_INT, const float,
                         const float *, const MKL_INT, const float *,
                         const MKL_INT, const float, float *,
                         const MKL_INT)>
                (dlsym(mkl_intel_lp64, "cblas_sgemv"));
            CHECK_HANDLE(_sgemv);
            _ssymv =
                reinterpret_cast<
                void (*)(const CBLAS_LAYOUT, const CBLAS_UPLO,
                         const MKL_INT, const float, const float *,
                         const MKL_INT, const float *, const MKL_INT,
                         const float, float *, const MKL_INT)>
                (dlsym(mkl_intel_lp64, "cblas_ssymv"));
            CHECK_HANDLE(_ssymv);
            _sger =
                reinterpret_cast<
                void (*)(const CBLAS_LAYOUT, const MKL_INT,
                         const MKL_INT, const float, const float *,
                         const MKL_INT, const float *, const MKL_INT,
                         float *, const MKL_INT)>
                (dlsym(mkl_intel_lp64, "cblas_sger"));
            CHECK_HANDLE(_sger);
            _sgemm =
                reinterpret_cast<
                void (*)(const CBLAS_LAYOUT, const CBLAS_TRANSPOSE,
                         const CBLAS_TRANSPOSE, const MKL_INT,
                         const MKL_INT, const MKL_INT, const float,
                         const float *, const MKL_INT, const float *,
                         const MKL_INT, const float, float *,
                         const MKL_INT)>
                (dlsym(mkl_intel_lp64, "cblas_sgemm"));
            CHECK_HANDLE(_sgemm);
            _ssymm =
                reinterpret_cast<
                void (*)(const CBLAS_LAYOUT, const CBLAS_SIDE,
                         const CBLAS_UPLO, const MKL_INT,
                         const MKL_INT, const float, const float *,
                         const MKL_INT, const float *, const MKL_INT,
                         const float, float *, const MKL_INT)>
                (dlsym(mkl_intel_lp64, "cblas_ssymm"));
            CHECK_HANDLE(_ssymm);
            _vsmul =
                reinterpret_cast<
                void (*)(const MKL_INT, const float [],
                         const float [], float [])>
                (dlsym(mkl_intel_lp64, "vsMul"));
            CHECK_HANDLE(_vsmul);
        }
        std::string version_str(void)
        {
            // The maximum buffer size MKL_Get_Version_String can
            // handle appears to be 4k, larger buffer will cause
            // MKL_Get_Version_String to silently return (and leaving
            // an empty string)
            char buf[4096] = { '\0' };

            if (_get_version_string != NULL) {
                (*_get_version_string)(buf, 4096);
            }
            return std::string(buf);
        }
    };
}
