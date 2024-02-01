cmake_minimum_required(VERSION 3.26)

project(genomic_tools)

find_package(Boost 1.82.0 COMPONENTS iostreams REQUIRED)
find_package(fmt REQUIRED)

add_executable(genomic_validator
	genomic_validator.cxx
)

target_compile_features(genomic_validator PRIVATE cxx_std_23)
target_link_libraries(genomic_validator PRIVATE ${Boost_LIBRARIES} fmt::fmt)

