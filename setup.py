import os
import sys
import pathlib
import warnings
import sysconfig
from setuptools import setup, Extension
from sysconfig import get_config_vars
import subprocess
from Cython.Build import cythonize


base = pathlib.Path(__file__).parent.resolve()
try:
    compiler = get_config_vars()["CXX"].split(" ")[0]
except KeyError:
    # not sure why this happens
    compiler = "msvc"
BASE_FLAGS = ["-Os", "-std=c++17"]


if sys.platform == "darwin":
    osys = "mac"
    OMP_LINKER = ["-fopenmp", "-lomp"]
    if compiler in {"gcc", "g++"}:
        OMP_FLAGS = ["-Xpreprocessor", "-fopenmp"]
    else:
        OMP_FLAGS = ["-Xpreprocessor", "-fopenmp"]
elif "win" in sys.platform or sys.platform == "msys":
    osys = "windows"
    # TODO: change once openmp is supported for windows
    OMP_FLAGS = None
    OMP_LINKER = None
else:
    OMP_FLAGS = ["-fopenmp"]
    OMP_LINKER = ["-fopenmp"]
    osys = "linux"


class CompileError(Exception):
    message: str


def conda_env():
    return os.path.exists(os.path.join(sys.prefix, "conda-meta"))


def find_boost(ignore_conda=False):
    # first try to find the boost submodule
    base_path = pathlib.Path(__file__).parent / "pymantra"
    base_path = base_path / "network" / "enrichment" / "boost"
    if base_path.is_dir():
        return str(base_path.absolute()), None

    if conda_env() and not ignore_conda:
        library_dir = sysconfig.get_config_var("LIBDIR")
        basepath = pathlib.Path(library_dir).parent.resolve()
        if "Lib" in library_dir:
            include_dir = str(basepath / "Include")
        else:
            include_dir = str(basepath / "include")
        if os.path.isdir(os.path.join(include_dir, "boost")):
            print("boost located in conda environment")
            return include_dir, library_dir
    if osys != "windows":
        # NOTE: this should also cover the typical macOS path
        include_dir = "/usr/local/include"
        if os.path.isdir(include_dir + os.path.sep + "boost"):
            if os.path.exists("/usr/local/lib/x86_64-linux-gnu/libboost*"):
                return include_dir, "/usr/local/lib/x86_64-linux-gnu"
            elif os.path.exists("/usr/local/lib64/libboost*"):
                return include_dir, "/usr/local/lib64"
            elif os.path.exists("/usr/local/lib/libboost*"):
                return include_dir, "/usr/local/lib"
        include_dir = "/usr/include"
        if os.path.isdir(include_dir + os.path.sep + "boost"):
            if os.path.join("/usr/lib/x86_64-linux-gnu/libboost*"):
                return include_dir, "/usr/lib/x86_64-linux-gnu"
            elif os.path.join("/usr/lib64/libboost*"):
                return include_dir, "/usr/lib64"
            elif os.path.join("/usr/lib/libboost*"):
                return include_dir, "/usr/lib"
    print(
        "boost library could not be located, please install it either via "
        "conda (conda install -c conda-forge boost) or manually install it "
        "and specify the path using --include-dir and --library-dir if the "
        "installation fails."
    )
    return None, None


# ============================ #
# Compile requirements testing #
# ============================ #
def _run_compile_test(call, exec_name, err_message, raise_warning):
    process = subprocess.Popen(call, stdout=subprocess.PIPE)
    _ = process.communicate()[0]
    if process.returncode != 0:
        if raise_warning:
            warnings.warn(err_message)
        return False

    test_call = [f"./{exec_name}"]
    process = subprocess.Popen(test_call, stdout=subprocess.PIPE)
    _ = process.communicate()[0]

    os.remove(exec_name)
    return process.returncode == 0


def boost_test(include_dir, lib_dir):
    base_call = get_config_vars()["CXX"].split(" ") + BASE_FLAGS

    if include_dir is not None:
        if compiler == "msvc":
            # according to microsoft's documentation the dash should work
            # too, but didn't work in tests. mingw should work with dashes
            base_call = ["/I", include_dir]
        else:
            base_call += ["-I", include_dir]
    if lib_dir is not None:
        if compiler == "msvc":
            base_call = ["/L", lib_dir]
        else:
            base_call += ["-L", lib_dir]
    compile_call = base_call + [str(base / "pymantra" / "boost_available.cpp")]

    test_command = "boost_test"
    compile_call += ["-o", test_command]
    passed = _run_compile_test(compile_call, test_command, "", False)

    if not passed:
        fail_message = \
            "Unable to compile correctly with boost location! " \
            f"Searched location: {include_dir}. We recommend either " \
            "installing boost (>= 1.77) with conda via conda-forge or " \
            "to manually install it. If you cannot resolve this issue " \
            "please contact the developers to figure out the best option " \
            "for your case."
        base_call.append(str(base / "pymantra" / "boost_version.cpp"))
        base_call += ["-o", "boost_version"]

        process = subprocess.Popen(base_call, stdout=subprocess.PIPE)
        _ = process.communicate()[0]
        if process.returncode != 0:
            raise CompileError(fail_message)

        test_call = ["./boost_version"]
        process = subprocess.Popen(test_call, stdout=subprocess.PIPE)
        output = process.communicate()[0]

        if process.returncode != 0:
            raise CompileError(fail_message)
        os.remove("boost_version")

        version_str = output.decode("utf-8").strip()
        version = [int(x) for x in version_str.split("_")]
        if version[0] < 2:
            if version[0] < 1 or version[1] < 77:
                raise CompileError(
                    f"Unsupported boost version {version_str}. Minimum "
                    "version of 1_77 required. Please update boost and run "
                    "the installation again."
                )


def build_with_openmp():
    warning = "Failed to compile with OpenMP. Proceeding without " \
              "parallelized implementations."
    try:
        compile_command = get_config_vars()["CXX"]
    except KeyError:
        warnings.warn(warning)
        return False

    compile_call = compile_command.split(" ") + BASE_FLAGS + OMP_FLAGS
    compile_call.append(str(base / "pymantra" / "openmp_available.cpp"))
    test_command = "openmp_test"
    compile_call += OMP_LINKER + ["-o", test_command]

    return _run_compile_test(compile_call, test_command, warning, True)


def contains_index(flag):
    for i, part in enumerate(sys.argv):
        if part == flag:
            return i
    return None


# ======================= #
# Building c++ extensions #
# ======================= #
boost_include = None
boost_lib = None
include_idx = contains_index("--include-dir")
if include_idx is not None:
    boost_include = sys.argv[include_idx + 1]
    # NOTE: for newer setuptools versions it is required to remove
    #       all self-defined flags from sys.argv before running setup()
    del sys.argv[include_idx]
    del sys.argv[include_idx]

lib_idx = contains_index("--library-dir")
if lib_idx is not None:
    boost_lib = sys.argv[lib_idx + 1]
    del sys.argv[lib_idx]
    del sys.argv[lib_idx]

    if boost_lib == "None":
        boost_lib = None

if boost_include is None:
    conda_idx = contains_index("--ignore-conda")
    if conda_idx is not None:
        print(
            "Found --ignore-conda flag, trying to find system-location "
            "of boost"
        )
        del sys.argv[conda_idx]
    if boost_lib is None:
        boost_include, boost_lib = find_boost(
            ignore_conda=conda_idx is not None)
    else:
        boost_include, _ = find_boost(
            ignore_conda=conda_idx is not None)
if not (osys == "windows" and os.environ.get("GITHUB_ACTION") == "true"):
    # this avoids having to manually call the compiler in a github action
    boost_test(boost_include, boost_lib)


if osys == "windows":
    warnings.warn(
        "OpenMP currently not supported for Windows. Proceeding "
        "without parallelized implementations."
    )
    if compiler == "msvc":
        # TODO: reset once figure out how to use openmp with standards < 2
        # OPENMP_FLAG = ["/openmp", "/Os", "/std:c++17"]
        FLAGS = ["/std:c++17"]
        LINKER = None
    elif compiler.startswith("mingw"):
        # TODO: reset once figure out how to use openmp with standards < 2
        # OPENMP_FLAG = ["-fopenmp", "-Os", "-std=c++17"]
        # OPENMP_LINK = ["-fopenmp"]
        FLAGS = ["-std:c++17"]
        LINKER = None
    else:
        raise ValueError(
            f"Unsupported compiler {compiler}. "
            "Please use msvc or mingw"
        )
elif build_with_openmp():
    # check if compilation with openmp works
    # if yes add openmp flags + macro
    FLAGS = ["-DUSE_OPENMP"] + OMP_FLAGS + BASE_FLAGS
    LINKER = OMP_LINKER
else:
    FLAGS = BASE_FLAGS
    LINKER = None


enrichment_path = "./pymantra/network/enrichment"
# ============ #
# Local Search #
# ============ #
lso_files = [
    "Exceptions.cpp", "pyutils.cpp", "statsutils.cpp",
    os.path.join("LSO", "LocalSearch.cpp"),
    os.path.join("LSO", "lso_utils.cpp"),
    os.path.join("LSO", "objective_functions.cpp"),
    os.path.join("LSO", "reaction_graph.cpp")
]
lso_files = [os.path.join(enrichment_path, file) for file in lso_files]
if boost_include is not None:
    include_dirs = [boost_include, enrichment_path]
else:
    include_dirs = [enrichment_path]
lso_module = Extension(
    "pymantra.network.enrichment.LSO.lso",
    sources=[os.path.join(enrichment_path, "LSO", "lso.pyx")] + lso_files,
    include_dirs=include_dirs,
    library_dirs=[boost_lib] if boost_lib is not None else None,
    extra_compile_args=FLAGS,
    extra_link_args=LINKER,
)

# ======== #
# Spearman #
# ======== #
spearman_files = [
    os.path.join(enrichment_path, file) for file in
    ["Exceptions.cpp", "pyutils.cpp", "statsutils.cpp"]
]
spearman_module = Extension(
    "pymantra.network.enrichment.spearman",
    sources=[os.path.join(enrichment_path, "spearman.pyx")] + spearman_files,
    include_dirs=include_dirs,
    library_dirs=[boost_lib] if boost_lib is not None else None,
    extra_compile_args=FLAGS,
    extra_link_args=LINKER
)

# ======================= #
# Building actual package #
# ======================= #
lso_ext = cythonize(lso_module, compiler_directives={"language_level": "3"})
spearman_ext = cythonize(
    spearman_module, compiler_directives={"language_level": "3"})

try:
    # NOTE: for some reason each extension is a list
    setup(ext_modules=lso_ext + spearman_ext)
except RuntimeError:
    raise RuntimeError(
        "pymantra installation failed, most likely due to errors during "
        "the compilation process. Please make sure you have boost "
        "installed. If you are on windows, also make sure you have the "
        "Visual Studio C++ extension installed. If you know the location of "
        "your boost installation you can pass it manually using the "
        "--include-dir and --library-dir flags. In case you are using conda, "
        "boost can best installed using conda install -c conda-forge boost."
    )
