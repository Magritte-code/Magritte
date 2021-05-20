#pragma once


#include <limits>

#include "io/io.hpp"
#include "tools/setOnce.hpp"




///  Create for each parameter "x":
///    - a (private) local variable "x__",
///    - a (private) Meyers' singleton for isSet "x_isSet()",
///    - a (public) Meyers' singleton for value "x()",
///    - a (public) setter function "set_x".
///  Remember to ensure in Parameters copy constructor to
///  ensure that the local variables are up to date, i.e.
///  ensure that their values equal the global values.

#define CREATE_PARAMETER(type, x)                                                           \
    private:                                                                                \
        SetOnce<type> x##__;                     /* Local (SetOnce) variable           */   \
        accel inline bool& x##_isSet () const    /* Meyers' singleton for global isSet */   \
        {                                                                                   \
            static bool isSet = false;           /* Global isSet                       */   \
            return isSet;                        /* Return global isSet                */   \
        }                                                                                   \
    public:                                                                                 \
        accel inline type& x () const            /* Meyers' singleton for global value */   \
        {                                                                                   \
            static type x##_value;               /* Global value                       */   \
            if (!x##_isSet())                    /* If global value is not set         */   \
            {                                                                               \
                x##_value   = x##__.get();       /* Set global with local value        */   \
                x##_isSet() = true;              /* Mark global value as set           */   \
            }                                                                               \
            return x##_value;                    /* Return global value                */   \
        }                                                                                   \
        inline void set_##x (const type value)   /* Setter function                    */   \
        {                                                                                   \
            x##__.set (value);                   /* Set local value                    */   \
            x() = value;                         /* Set global value                   */   \
        }                                                                                   \
        inline type get_##x () const             /* Getter function                    */   \
        {                                                                                   \
            return x();                          /* Return copy of global value        */   \
        }

#define CONSTRUCT_PARAMETER(type, x)                                                        \
    x##__ = SetOnce <type> (x());

#define COPY_PARAMETER(x)                                                                   \
    if (x##_isSet()) x##__.set(x());



///  Parameters: secure structure for the model parameters
//////////////////////////////////////////////////////////
struct Parameters
{
//    private:
//        SetOnce <size_t> npoints_;                ///< number of points
//        SetOnce <size_t> totnnbs_;                ///< total number of neighbors
//        SetOnce <size_t> nrays_;                  ///< number of rays (originating from each cell)
//        SetOnce <size_t> hnrays_;                 ///< half the number of rays (originating from each cell)
//        SetOnce <size_t> nrays_red_;              ///< number of rays reduced
//        SetOnce <size_t> order_min_;              ///< lowest (HEALPix) order of adaptive rays
//        SetOnce <size_t> order_max_;              ///< highest (HEALPix) order of adaptive rays
//        SetOnce <size_t> nboundary_;              ///< number of points on the boundary
//        SetOnce <size_t> nfreqs_;                 ///< number of frequency bins
//        SetOnce <size_t> nfreqs_red_;             ///< number of frequency bins reduced
//        SetOnce <size_t> nspecs_;                 ///< number of chemical species
//        SetOnce <size_t> nlspecs_;                ///< number of line producing species
//        SetOnce <size_t> nlines_;                 ///< number of line transitions
//        SetOnce <size_t> nquads_;                 ///< number of frequency quadrature points

//        SetOnce <double> pop_prec_;               ///< required precision for populations

//        SetOnce <bool>   use_scattering_;         ///< true if scattering is used
//        SetOnce <bool>   use_Ng_acceleration_;    ///< true if Ng acceleration is used
//        SetOnce <bool>   spherical_symmetry_;     ///< true if spherical symmetry is used
//        SetOnce <bool>   adaptive_ray_tracing_;   ///< true if adaptive ray-tracing is used

//    public:
//        inline void set_npoints   (const size_t value) {   npoints_.set(value);}
//        inline void set_totnnbs   (const size_t value) {   totnnbs_.set(value);}
//        inline void set_nrays     (const size_t value) {     nrays_.set(value);}
//        inline void set_hnrays    (const size_t value) {    hnrays_.set(value);}
//        inline void set_nrays_red (const size_t value) { nrays_red_.set(value);}
//        inline void set_order_min (const size_t value) { order_min_.set(value);}
//        inline void set_order_max (const size_t value) { order_max_.set(value);}
//        inline void set_nboundary (const size_t value) { nboundary_.set(value);}
//        inline void set_nfreqs    (const size_t value) {    nfreqs_.set(value);}
//        inline void set_nfreqs_red(const size_t value) {nfreqs_red_.set(value);}
//        inline void set_nspecs    (const size_t value) {    nspecs_.set(value);}
//        inline void set_nlspecs   (const size_t value) {   nlspecs_.set(value);}
//        inline void set_nlines    (const size_t value) {    nlines_.set(value);}
//        inline void set_nquads    (const size_t value) {    nquads_.set(value);}

//        inline void set_pop_prec(const double value) {pop_prec_.set(value);}

//        inline void set_use_scattering      (const bool value) {use_scattering_      .set(value);}
//        inline void set_use_Ng_acceleration (const bool value) {use_Ng_acceleration_ .set(value);}
//        inline void set_spherical_symmetry  (const bool value) {spherical_symmetry_  .set(value);}
//        inline void set_adaptive_ray_tracing(const bool value) {adaptive_ray_tracing_.set(value);}

//        accel inline size_t npoints   () const {return    npoints_.get();}
//        accel inline size_t totnnbs   () const {return    totnnbs_.get();}
//        accel inline size_t nrays     () const {return      nrays_.get();}
//        accel inline size_t hnrays    () const {return     hnrays_.get();}
//        accel inline size_t nrays_red () const {return  nrays_red_.get();}
//        accel inline size_t order_min () const {return  order_min_.get();}
//        accel inline size_t order_max () const {return  order_max_.get();}
//        accel inline size_t nboundary () const {return  nboundary_.get();}
//        accel inline size_t nfreqs    () const {return     nfreqs_.get();}
//        accel inline size_t nfreqs_red() const {return nfreqs_red_.get();}
//        accel inline size_t nspecs    () const {return     nspecs_.get();}
//        accel inline size_t nlspecs   () const {return    nlspecs_.get();}
//        accel inline size_t nlines    () const {return     nlines_.get();}
//        accel inline size_t nquads    () const {return     nquads_.get();}

//        accel inline double pop_prec() const {return pop_prec_.get();}

//        accel inline bool use_scattering      () const {return use_scattering_      .get();}
//        accel inline bool use_Ng_acceleration () const {return use_Ng_acceleration_ .get();}
//        accel inline bool spherical_symmetry  () const {return spherical_symmetry_  .get();}
//        accel inline bool adaptive_ray_tracing() const {return adaptive_ray_tracing_.get();}


    long n_off_diag = 0;

    bool writing_populations_to_disk=false;   ///< Toggle for writing the level populations after each iteration

    double max_width_fraction = 0.5;

    void read (const Io &io);
    void write(const Io &io) const;

    CREATE_PARAMETER (string, model_name);

    CREATE_PARAMETER (Size, dimension );
    CREATE_PARAMETER (Size, npoints   );
    CREATE_PARAMETER (Size, totnnbs   );
    CREATE_PARAMETER (Size, nrays     );
    CREATE_PARAMETER (Size, hnrays    );
    CREATE_PARAMETER (Size, nrays_red );
    CREATE_PARAMETER (Size, order_min );
    CREATE_PARAMETER (Size, order_max );
    CREATE_PARAMETER (Size, nboundary );
    CREATE_PARAMETER (Size, nfreqs    );
    CREATE_PARAMETER (Size, nspecs    );
    CREATE_PARAMETER (Size, nlspecs   );
    CREATE_PARAMETER (Size, nlines    );
    CREATE_PARAMETER (Size, nquads    );

    CREATE_PARAMETER (Real, pop_prec);

    CREATE_PARAMETER (bool, use_scattering      );
    CREATE_PARAMETER (bool, use_Ng_acceleration );
    CREATE_PARAMETER (bool, spherical_symmetry  );
    CREATE_PARAMETER (bool, adaptive_ray_tracing);

    Parameters ()
    {
        CONSTRUCT_PARAMETER (string, model_name);

        CONSTRUCT_PARAMETER (Size, dimension );
        CONSTRUCT_PARAMETER (Size, npoints   );
        CONSTRUCT_PARAMETER (Size, totnnbs   );
        CONSTRUCT_PARAMETER (Size, nrays     );
        CONSTRUCT_PARAMETER (Size, hnrays    );
        CONSTRUCT_PARAMETER (Size, nrays_red );
        CONSTRUCT_PARAMETER (Size, order_min );
        CONSTRUCT_PARAMETER (Size, order_max );
        CONSTRUCT_PARAMETER (Size, nboundary );
        CONSTRUCT_PARAMETER (Size, nfreqs    );
        CONSTRUCT_PARAMETER (Size, nspecs    );
        CONSTRUCT_PARAMETER (Size, nlspecs   );
        CONSTRUCT_PARAMETER (Size, nlines    );
        CONSTRUCT_PARAMETER (Size, nquads    );

        CONSTRUCT_PARAMETER (Real, pop_prec);

        CONSTRUCT_PARAMETER (bool, use_scattering      );
        CONSTRUCT_PARAMETER (bool, use_Ng_acceleration );
        CONSTRUCT_PARAMETER (bool, spherical_symmetry  );
        CONSTRUCT_PARAMETER (bool, adaptive_ray_tracing);
    }

    Parameters (const Parameters& parameters)
    {
        COPY_PARAMETER (model_name);

        COPY_PARAMETER (dimension );
        COPY_PARAMETER (npoints   );
        COPY_PARAMETER (totnnbs   );
        COPY_PARAMETER (nrays     );
        COPY_PARAMETER (hnrays    );
        COPY_PARAMETER (nrays_red );
        COPY_PARAMETER (order_min );
        COPY_PARAMETER (order_max );
        COPY_PARAMETER (nboundary );
        COPY_PARAMETER (nfreqs    );
        COPY_PARAMETER (nspecs    );
        COPY_PARAMETER (nlspecs   );
        COPY_PARAMETER (nlines    );
        COPY_PARAMETER (nquads    );

        COPY_PARAMETER (pop_prec);

        COPY_PARAMETER (use_scattering      );
        COPY_PARAMETER (use_Ng_acceleration );
        COPY_PARAMETER (spherical_symmetry  );
        COPY_PARAMETER (adaptive_ray_tracing);
    }
};
