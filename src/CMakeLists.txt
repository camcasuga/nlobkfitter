configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gitsha1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp" @ONLY)

add_executable(
	plottool
	tool_plot.cpp
	massive_swarmscan.cpp
	dipole.cpp
	data.cpp
	solver.cpp
	mv.cpp
	ic.cpp
	ic_datafile.cpp
	nlobk_config.cpp
	nlodissigmar.cpp	
	nlodissigmar_massiveq.cpp	
	nlodis_config.cpp
	helper.cpp
	${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp
)

add_executable(
	nlofit	
	#main.cpp
	massive_swarmscan.cpp
	dipole.cpp
	data.cpp
	solver.cpp
	mv.cpp
	ic.cpp
	ic_datafile.cpp
	nlobk_config.cpp
	nlodissigmar.cpp	
	nlodissigmar_massiveq.cpp	
	nlodis_config.cpp
	helper.cpp
	${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp
)
	
include_directories(.)

find_package(GSL REQUIRED)    #
target_link_libraries(
	nlofit
	PRIVATE
	cuba
	Minuit2
	amplitude
	GSL::gsl
	GSL::gslcblas
)


