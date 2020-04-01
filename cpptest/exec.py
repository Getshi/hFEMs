import os
import sys
import subprocess
import argparse
import shutil


def build(sourcedir, builddir, target=[], debug=False, toolchain=None):
    cfg = 'Debug' if debug else 'Release'
    cmake_args = ['-DCMAKE_BUILD_TYPE=' + cfg]
    if not toolchain is None:
        cmake_args += ['-DCMAKE_TOOLCHAIN_FILE=' + toolchain]
    build_args = ['--config', cfg]
    build_args += ['--', '-j']

    os.makedirs(builddir, exist_ok=True)

    subprocess.check_call(['cmake', sourcedir] +
                          cmake_args, cwd=builddir)
    subprocess.check_call(['cmake', '--build', '.'] +
                          build_args + [target], cwd=builddir)

# get command line arguments
ap = argparse.ArgumentParser()
ap.add_argument("-d", "--debug", action='store_true',
                help="...")
ap.add_argument("-b", "--build", default="1",
                help="...")
ap.add_argument("-c", "--clean", action='store_true',
                help="...")

args = ap.parse_args()
args.build = args.build != "0"
target = "cpptest"

# BUILD
sourcedir = os.getcwd()  # where the main CMakeLists.txt is
builddir = os.path.join(
    os.getcwd(), "build-Debug" if args.debug else "build-Release")


if args.clean:
    shutil.rmtree(builddir)

if args.build:
    toolchain = "../vcpkg/scripts/buildsystems/vcpkg.cmake"
    build(sourcedir=sourcedir, builddir=builddir, target=target,
        debug=args.debug, toolchain=toolchain)

# RUN
workdir = os.getcwd() # run from parent
executable = os.path.join(builddir, target)

# TODO this doesnt give segfault
print("Executing:", executable)
try:
    subprocess.run(executable, cwd=workdir, shell=True)
except subprocess.CalledProcessError as exc:
    print("PY: CalledProcessError")
    print(exc.returncode, exc.output, exc.stderr)
except KeyboardInterrupt:
    print("PY: Aborting execution")