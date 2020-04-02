configure_file(dualfermion-config.cmake.in dualfermion-config.cmake @ONLY)
configure_file(dualfermion-config-version.cmake.in dualfermion-config-version.cmake @ONLY)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/dualfermion-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/dualfermion-config-version.cmake
  DESTINATION lib/cmake/dualfermion
)

install(EXPORT dualfermion-targets NAMESPACE dualfermion:: DESTINATION lib/cmake/dualfermion)
