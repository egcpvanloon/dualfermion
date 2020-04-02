# Generated automatically using the command :
# c++2py ../../c++/app4triqs/post_process.hpp --members_read_only -N app4triqs -a app4triqs -m post_process -o post_process -C pytriqs --moduledoc="The app4triqs postprocess functionality" --cxxflags="-std=c++17" --target_file_only
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "post_process", doc = r"The app4triqs postprocess functionality", app_name = "app4triqs")

# Imports

# Add here all includes
module.add_include("app4triqs/post_process.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""

using namespace app4triqs;
""")




module.generate_code()