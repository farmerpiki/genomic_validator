cmake_minimum_required(VERSION 3.26)

project(cpp_template)

#find_package(Boost 1.82.0 COMPONENTS program_options system json filesystem REQUIRED)
#find_package(PostgreSQL REQUIRED)
#find_package(OpenSSL REQUIRED)
find_package(fmt REQUIRED)

add_executable(${PROJECT_NAME}
	main.cxx
)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_23)
#target_link_libraries(${PROJECT_NAME} PRIVATE ${Boost_LIBRARIES} ${PostgreSQL_LIBRARIES} fmt::fmt)
target_link_libraries(${PROJECT_NAME} PRIVATE fmt::fmt)

