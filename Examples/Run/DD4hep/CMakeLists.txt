add_executable(
  ActsExampleDD4hepDetectorGdml
  DD4hepDetectorGdmlExample.cpp)
target_link_libraries(
  ActsDD4hepDetectorGdml
  PRIVATE
    ActsCore ActsExamplesCommon
    ActsExamplesFramework
    Boost::program_options)

install(
  TARGETS ActsExampleDD4hepDetectorGdml
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
