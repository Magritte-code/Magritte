#pragma once


#include "io/io.hpp"
#include "io/python/io_python.hpp"
#include "parameters/parameters.hpp"
#include "tools/types.hpp"
#include "tools/timer.hpp"
#include "geometry/geometry.hpp"
#include "chemistry/chemistry.hpp"
#include "thermodynamics/thermodynamics.hpp"
#include "lines/lines.hpp"
#include "radiation/radiation.hpp"
#include "image/image.hpp"


struct Model
{
    std::shared_ptr<Parameters> parameters;

    Geometry       geometry;
    Chemistry      chemistry;
    Thermodynamics thermodynamics;
    Lines          lines;
    Radiation      radiation;
    vector<Image>  images;

    enum SpectralDiscretisation {SD_None, SD_Lines, SD_Image}
         spectralDiscretisation = SD_None;


    Model ()
    : parameters     (new Parameters())
    , geometry       (    parameters  )
    , chemistry      (    parameters  )
    , thermodynamics (    parameters  )
    , lines          (    parameters  )
    , radiation      (    parameters  )
    {};

    Model (const string name)
    : parameters     (new Parameters())
    , geometry       (    parameters  )
    , chemistry      (    parameters  )
    , thermodynamics (    parameters  )
    , lines          (    parameters  )
    , radiation      (    parameters  )
    {
        parameters->set_model_name (name);
        read ();
    }

    void read  (const Io& io);
    void write (const Io& io) const;

    void read  ()       {read  (IoPython ("hdf5", parameters->model_name()));};
    void write () const {write (IoPython ("hdf5", parameters->model_name()));};

    void read  (const string model_name)       {read  (IoPython ("hdf5", model_name));};
    void write (const string model_name) const {write (IoPython ("hdf5", model_name));};

    int compute_inverse_line_widths               ();
    int compute_spectral_discretisation           ();
    int compute_spectral_discretisation           (
        const Real width );
    int compute_spectral_discretisation           (
        const long double nu_min,
        const long double nu_max );
    int compute_LTE_level_populations             ();
    int compute_radiation_field                   ();
    int compute_radiation_field_feautrier_order_2 ();
    int compute_radiation_field_shortchar_order_0 ();
    int compute_Jeff                              ();
    int compute_Jeff_sparse                       ();
    int compute_level_populations_from_stateq     ();
    int compute_level_populations                 (
        const bool  use_Ng_acceleration,
        const long  max_niterations     );
    int compute_level_populations_sparse          (
        const bool  use_Ng_acceleration,
        const long  max_niterations     );
    int compute_image                             (const Size ray_nr);
    int compute_image_optical_depth               (const Size ray_nr);

    Double1 error_max;
    Double1 error_mean;

    Matrix<Real> eta;
    Matrix<Real> chi;

    Matrix<Real>  eta_ray;
    Matrix<Real>  chi_ray;
    Matrix<Real> dtau_ray;
    Matrix<Real>    u_ray;

    Matrix<Real> boundary_condition;

    int compute_image_for_point (const Size ray_nr, const Size p);

    int compute_radiation_field_feautrier_order_2_uv     ();
    int compute_radiation_field_feautrier_order_2_anis   ();
    int compute_radiation_field_feautrier_order_2_sparse ();

    int set_eta_and_chi       (const Size rr);
    int set_boundary_condition();

    Vector<Real> density;
    Matrix<Real> column;

    int set_column ();


    // PORTAL extension
    ////////////////////////////////
    Size  nalign;
    Size1 lau;
    Size1 lal;
    Size1 lcu;
    Size1 lcl;
    Size1 j_lev;
    Size2 a_lev;
    Double2 J0;
    Double2 J2;
    Double5 rp;
    Double5 rm;
    Double2 tp;
    Double1 t_A;
    Double1 up_loc;
    Double1 down_loc;

    MatrixXd M;

    Vector <Vector3D> b;

    Double2 population_align;

    Real1 k_abs_0;
    Real1 k_stm_0;
    Real1 k_abs_2;
    Real1 k_stm_2;


    int PORTAL_solve_statistical_equilibrium ();
    int PORTAL_image                         (const Size ray_nr, const Size l);
    ////////////////////////////////
};
