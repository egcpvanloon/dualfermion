from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "dpt_core", doc = "dual fermion runr", app_name = "triqs_dualfermion")

# Imports
#import pytriqs.atom_diag
import pytriqs.gf
import pytriqs.operators
import pytriqs.statistics.histograms

# Add here all includes
module.add_include("triqs_dualfermion/dpt_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <cpp2py/converters/map.hpp>
#include <cpp2py/converters/optional.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/set.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/variant.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/operators_real_complex.hpp>
#include <triqs/cpp2py_converters/h5.hpp>
#include <boost/mpi.hpp>

using namespace triqs_dualfermion;
""")

module.add_enum("block_order", ['block_order::AABB', 'block_order::ABBA'], "triqs_dualfermion", """Order of block indices for Block2Gf objects""")

# The class dpt_core
c = class_(
        py_type = "DptCore",  # name of the python class
        c_type = "triqs_dualfermion::dpt_core",   # name of the C++ class
        doc = """Core class of the dualfermion program""",   # doc of the C++ class
        hdf5 = True,
)

c.add_constructor("""(**triqs_dualfermion::constr_parameters_t)""", doc = """Initialize the dual perturbation theory""")

c.add_method("""void run (**triqs_dualfermion::run_parameters_t)""",
             doc = """Run the dual perturbation theory once.""")

c.add_method("""std::string hdf5_scheme ()""",
             is_static = True,
             doc = """""")

c.add_property(name = "last_constr_parameters",
               getter = cfunction("triqs_dualfermion::constr_parameters_t last_constr_parameters ()"),
               doc = """Set of parameters used in the construction of the ``dpt_core`` class.""")

c.add_property(name = "last_run_parameters",
               getter = cfunction("triqs_dualfermion::run_parameters_t last_run_parameters ()"),
               doc = """Set of parameters used in the last call to ``run()``.""")

c.add_property(name = "gimp",
               getter = cfunction("block_gf_view<triqs::gfs::imfreq> gimp ()"),
               doc = """:math:`g(i\\omega)` in imaginary frequencies.""")

c.add_property(name = "Delta",
               getter = cfunction("block_gf_view<triqs::gfs::imfreq> Delta ()"),
               doc = """:math:`\\Delta(i\\omega)` in imaginary frequencies.""")

c.add_property(name = "G2_iw",
               getter = cfunction("block2_gf_view<imfreq_cube_mesh_t, tensor_valued<4>> G2_iw ()" ),
               doc = """Two-particle Green\'s function :math:`G^{(2)}(i\\omega,i\\nu\',i\\nu\'\')` (1 bosonic, 2 fermionic frequencies)""")

c.add_property(name = "Hk",
               getter = cfunction("block_gf_view<triqs::gfs::brillouin_zone> Hk ()"),
               doc = """:math:`H(k)` on the Brillouin Zone""")

c.add_property(name = "status",
               getter = cfunction("int status ()"),
               doc = """Status of the ``run()`` on exit.""")

module.add_class(c)


# Converter for run_parameters_t
c = converter_(
        c_type = "triqs_dualfermion::run_parameters_t",
        doc = """""",
)
c.add_member(c_name = "calculate_sigma",
             c_type = "bool",
             initializer = """ false """,
             doc = """Will the dual self-energy be calculated?""")
c.add_member(c_name = "calculate_sigma1",
             c_type = "bool",
             initializer = """ true """,
             doc = """Will the first-order dual self-energy be calculated?""")
c.add_member(c_name = "calculate_sigma2",
             c_type = "bool",
             initializer = """ true """,
             doc = """Will the second-order dual self-energy be calculated?""")
c.add_member(c_name = "sigmad_subset",
             c_type = "triqs::hilbert_space::gf_struct_t",
             initializer = """ """,
             doc = """Subset of indices to use for the evaluation of sigmad""")
c.add_member(c_name = "verbosity",
             c_type = "int",
             initializer = """ ((boost::mpi::communicator().rank()==0)?3:0) """,
             doc = """Verbosity level\n     default: 3 on MPI rank 0, 0 otherwise.""")
c.add_member(c_name = "delta_initial",
             c_type = "bool",
             initializer = """ false """,
             doc = """Should we only calculate the initial guess for Delta?""")


module.add_converter(c)

# Converter for constr_parameters_t
c = converter_(
        c_type = "triqs_dualfermion::constr_parameters_t",
        doc = """""",
)
c.add_member(c_name = "beta",
             c_type = "double",
             initializer = """  """,
             doc = """Inverse temperature""")

c.add_member(c_name = "gf_struct",
             c_type = "triqs::hilbert_space::gf_struct_t",
             initializer = """  """,
             doc = """block structure of the gf""")

c.add_member(c_name = "Hk",
             c_type = "triqs::gfs::block_gf<brillouin_zone, matrix_valued>",
             initializer = """  """,
             doc = """Dispersion""")

c.add_member(c_name = "n_iw",
             c_type = "int",
             initializer = """ 128 """,
             doc = """Number of Matsubara frequencies for gf<imfreq, matrix_valued>""")

c.add_member(c_name = "n_iw2",
             c_type = "int",
             initializer = """ 32 """,
             doc = """Number of Fermionic Matsubara frequencies for G2""")

c.add_member(c_name = "n_iW",
             c_type = "int",
             initializer = """ 31 """,
             doc = """Number of Bosonic Matsubara frequencies for G2""")


module.add_converter(c)


module.generate_code()
