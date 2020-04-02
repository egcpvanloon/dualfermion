# Generated automatically using the command :
# c++2py ../../c++/app4triqs/solver_core.hpp -p --members_read_only -N app4triqs -a app4triqs -m solver_core -o solver_core --moduledoc="The app4triqs solve_core module" -C pytriqs --cxxflags="-std=c++17" --target_file_only
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "solver_core", doc = r"The app4triqs solve_core module", app_name = "app4triqs")

# Imports
module.add_imports(*['pytriqs.gf', 'pytriqs.operators'])

# Add here all includes
module.add_include("app4triqs/solver_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/optional.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/variant.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/operators_real_complex.hpp>
#include <triqs/cpp2py_converters/real_or_complex.hpp>
#include <triqs/cpp2py_converters/h5.hpp>

using namespace app4triqs;
""")


# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "app4triqs::solver_core",   # name of the C++ class
        doc = r"""The Solver class""",   # doc of the C++ class
        hdf5 = True,
)

c.add_member(c_name = "G_tau",
             c_type = "app4triqs::g_tau_t",
             read_only= True,
             doc = r"""Greens function in imaginary time""")

c.add_member(c_name = "G_iw",
             c_type = "app4triqs::g_iw_t",
             read_only= True,
             doc = r"""Greens function in Matsubara frequencies""")

c.add_member(c_name = "Sigma_iw",
             c_type = "app4triqs::g_iw_t",
             read_only= True,
             doc = r"""Self-energy in Matsubara frequencies""")

c.add_member(c_name = "constr_params",
             c_type = "app4triqs::constr_params_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "last_solve_params",
             c_type = "std::optional<solve_params_t>",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "G0_iw",
             c_type = "app4triqs::g_iw_t",
             read_only= True,
             doc = r"""Noninteracting Green Function in Matsubara frequencies""")

c.add_constructor("""(**app4triqs::constr_params_t)""", doc = r"""Construct a APP4TRIQS solver



+----------------+-----------------------------------+---------+----------------------------------+
| Parameter Name | Type                              | Default | Documentation                    |
+================+===================================+=========+==================================+
| n_tau          | int                               | 5001    | Number of tau points             |
+----------------+-----------------------------------+---------+----------------------------------+
| n_iw           | int                               | 500     | Number of Matsubara frequencies  |
+----------------+-----------------------------------+---------+----------------------------------+
| beta           | double                            | --      | Inverse temperature              |
+----------------+-----------------------------------+---------+----------------------------------+
| gf_struct      | triqs::hilbert_space::gf_struct_t | --      | Block structure of the gf        |
+----------------+-----------------------------------+---------+----------------------------------+
""")

c.add_method("""void solve (**app4triqs::solve_params_t)""",
             doc = r"""Solve method that performs APP4TRIQS calculation



+----------------+--------------------------------------+-----------------------------------+--------------------------------------------------+
| Parameter Name | Type                                 | Default                           | Documentation                                    |
+================+======================================+===================================+==================================================+
| h_int          | triqs::operators::many_body_operator | --                                | Interaction Hamiltonian                          |
+----------------+--------------------------------------+-----------------------------------+--------------------------------------------------+
| max_time       | int                                  | -1                                | Maximum running time in seconds (-1 : no limit)  |
+----------------+--------------------------------------+-----------------------------------+--------------------------------------------------+
| verbosity      | int                                  | mpi::communicator().rank()==0?3:0 | Verbosity                                        |
+----------------+--------------------------------------+-----------------------------------+--------------------------------------------------+
| post_process   | bool                                 | true                              | Perform post processing                          |
+----------------+--------------------------------------+-----------------------------------+--------------------------------------------------+
""")

c.add_method("""std::string hdf5_scheme ()""",
             is_static = True,
             doc = r"""""")

c.add_property(name = "post_process",
               getter = cfunction("void post_process ()"),
               doc = r"""""")

module.add_class(c)


# Converter for solve_params_t
c = converter_(
        c_type = "app4triqs::solve_params_t",
        doc = r"""The parameters for the solve function""",
)
c.add_member(c_name = "h_int",
             c_type = "triqs::operators::many_body_operator",
             initializer = """  """,
             doc = r"""Interaction Hamiltonian""")

c.add_member(c_name = "max_time",
             c_type = "int",
             initializer = """ -1 """,
             doc = r"""Maximum running time in seconds (-1 : no limit)""")

c.add_member(c_name = "verbosity",
             c_type = "int",
             initializer = """ mpi::communicator().rank()==0?3:0 """,
             doc = r"""Verbosity""")

c.add_member(c_name = "post_process",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Perform post processing""")

module.add_converter(c)

# Converter for constr_params_t
c = converter_(
        c_type = "app4triqs::constr_params_t",
        doc = r"""The parameters for the solver construction""",
)
c.add_member(c_name = "n_tau",
             c_type = "int",
             initializer = """ 5001 """,
             doc = r"""Number of tau points""")

c.add_member(c_name = "n_iw",
             c_type = "int",
             initializer = """ 500 """,
             doc = r"""Number of Matsubara frequencies""")

c.add_member(c_name = "beta",
             c_type = "double",
             initializer = """  """,
             doc = r"""Inverse temperature""")

c.add_member(c_name = "gf_struct",
             c_type = "triqs::hilbert_space::gf_struct_t",
             initializer = """  """,
             doc = r"""Block structure of the gf""")

module.add_converter(c)


module.generate_code()