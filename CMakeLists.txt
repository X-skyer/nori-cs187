cmake_minimum_required (VERSION 2.8.3)
project(nori)

add_subdirectory(ext ext_build)

include_directories(
  # Nori include files
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  # tinyformat string formatting library
  ${TFM_INCLUDE_DIR}
  # Eigen linear algebra library
  ${EIGEN_INCLUDE_DIR}
  # OpenEXR high dynamic range bitmap library
  ${OPENEXR_INCLUDE_DIRS}
  # Intel Thread Building Blocks
  ${TBB_INCLUDE_DIR}
  # Pseudorandom number generator
  ${PCG32_INCLUDE_DIR}
  # PugiXML parser
  ${PUGIXML_INCLUDE_DIR}
  # Helper functions for statistical hypothesis tests
  ${HYPOTHESIS_INCLUDE_DIR}
  # GLFW library for OpenGL context creation
  ${GLFW_INCLUDE_DIR}
  # GLEW library for accessing OpenGL functions
  ${GLEW_INCLUDE_DIR}
  # NanoVG drawing library
  ${NANOVG_INCLUDE_DIR}
  # NanoGUI user interface library
  ${NANOGUI_INCLUDE_DIR}
  # Portable filesystem API
  ${FILESYSTEM_INCLUDE_DIR}
  ext
)

# The following lines build the main executable. If you add a source
# code file to Nori, be sure to include it in this list.
add_executable(nori

  # Header files
  include/nori/bbox.h
  include/nori/bitmap.h
  include/nori/block.h
  include/nori/bsdf.h
  include/nori/bvh.h
  include/nori/camera.h
  include/nori/color.h
  include/nori/common.h
  include/nori/complexior.h
  include/nori/distributions.h
  include/nori/dpdf.h
  include/nori/emitter.h
  include/nori/frame.h
  include/nori/integrator.h
  include/nori/medium.h
  include/nori/mesh.h
  include/nori/object.h
  include/nori/optionsparser.h
  include/nori/parser.h
  include/nori/phase.h
  include/nori/pointlight.h
  include/nori/proplist.h
  include/nori/ray.h
  include/nori/rfilter.h
  include/nori/sample.h
  include/nori/sampler.h
  include/nori/scene.h
  include/nori/timer.h
  include/nori/transform.h
  include/nori/vector.h
  include/nori/warp.h

  # Source code files
  src/arealight.cpp
  src/avintegrator.cpp
  src/bitmap.cpp
  src/block.cpp
  src/bvh.cpp
  src/chi2test.cpp
  src/common.cpp
  src/roughconductor.cpp
  src/smoothconductor.cpp
  src/dielectric.cpp
  src/diffuse.cpp  
  src/directlightingintegrator.cpp
  src/distributions.cpp
  src/emitters/distantdisk.h
  src/emitters/distantdisk.cpp
  src/environmentlight.cpp
  src/gui.cpp
  src/homogeneous.cpp
  src/independent.cpp
  src/integrators/direct_ems.cpp
  src/integrators/direct_mats.cpp
  src/integrators/direct_mis.cpp
  src/integrators/path_mats.cpp
  src/integrators/path_mis.cpp
  src/main.cpp
  src/medium.cpp
  src/mesh.cpp
  src/normals.cpp
  src/obj.cpp
  src/object.cpp
  src/parser.cpp
  src/perspective.cpp
  src/pointlight.cpp
  src/phase.cpp
  src/proplist.cpp
  src/render.cpp
  src/rfilter.cpp
  src/roughdielectric.cpp
  src/scene.cpp
  src/ttest.cpp
  src/warp.cpp
  src/microfacet.cpp
  src/photon.cpp
  src/mirror.cpp
  src/dielectric.cpp
  src/photonmapper.cpp
  src/arealight.cpp
  src/volpath_singlescattering.cpp
)


# The following lines build the warping test application
add_executable(warptest
  include/nori/warp.h
  src/warp.cpp
  src/warptest.cpp
  src/microfacet.cpp
  src/object.cpp
  src/proplist.cpp
  src/common.cpp
)

# The following lines build the tonemapper
add_executable(tonemapper
        include/nori/bitmap.h
        src/bitmap.cpp
        src/common.cpp
        src/hdrToLdr.cpp)

target_link_libraries(nori tbb_static pugixml IlmImf nanogui ${NANOGUI_EXTRA_LIBS})
target_link_libraries(warptest tbb_static nanogui ${NANOGUI_EXTRA_LIBS})
target_link_libraries(tonemapper IlmImf)

# vim: set et ts=2 sw=2 ft=cmake nospell:
