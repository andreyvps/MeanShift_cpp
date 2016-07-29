cmake_minimum_required (VERSION 2.8)

project (MeanShift)

add_library(meanshift STATIC MeanShift.cpp MeanShift.h)

install(FILES MeanShift.h DESTINATION include/MeanShift)

install(TARGETS meanshift EXPORT MeanShiftTargetsRelease
    CONFIGURATIONS Release
    ARCHIVE DESTINATION Release/lib
    INCLUDES DESTINATION include
)

install(TARGETS meanshift EXPORT MeanShiftTargetsDebug
    CONFIGURATIONS Debug
    ARCHIVE DESTINATION Debug/lib
    INCLUDES DESTINATION include
)

install(EXPORT MeanShiftTargetsRelease
    FILE MeanShiftTargets.cmake
    NAMESPACE MeanShift::
    DESTINATION cmake
    CONFIGURATIONS Release
)

install(EXPORT MeanShiftTargetsDebug
    FILE MeanShiftTargets.cmake
    NAMESPACE MeanShift::
    DESTINATION cmake
    CONFIGURATIONS Debug
)
