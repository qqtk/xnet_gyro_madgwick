cmake_minimum_required(VERSION 2.8.3)
project(xnet_gyro_madgwick)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  roscpp
  sensor_msgs
  std_msgs
  tf
  robbase_msg
)
#  robbase_msg # re' bitbake donot generate build/--/include/WheelSpeed.h 're'solved after changing gitCode'and'md5sum' y16m8d02'h1750pm'h1842pm

## Generate messages in the 'msg' folder
# add_message_files(
#  DIRECTORY
#  msg
#  FILES
#  RazorImu.msg
#)

## Generate added messages and services with any dependencies listed here
# generate_messages(
#  DEPENDENCIES
#  std_msgs
#  sensor_msgs
# )

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if you package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
catkin_package(
  INCLUDE_DIRS include
  CATKIN_DEPENDS
        roscpp message_runtime std_msgs geometry_msgs
)
 # CATKIN_DEPENDS robbase_msg roscpp sensor_msgs std_msgs tf

###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
# 
# include_directories(include) # a=JunZhen
include_directories(include
  ${catkin_INCLUDE_DIRS}
)

# link_directories(/home/WrokGroup/catkin_ws/src/encoder_test/lib)

# add_library(encoder_test
#   src/${PROJECT_NAME}/encoder_test.cpp
# )

## Add cmake target dependencies of the library
## as an example, code may need to be generated before libraries
## either from message generation or dynamic reconfigure
# add_dependencies(encoder_test ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Declare a C++ executable
add_executable(xnet_gyro_madgwick_node src/xnet_gyro_madgwick.cpp)

## Add cmake target dependencies of the executable
## same as for the library above
# add_dependencies(encoder_test_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Specify libraries to link a library or executable target against
 target_link_libraries(xnet_gyro_madgwick_node
   ${catkin_LIBRARIES}
	libxNetDriver.so
 )
# add_dependencies(xnet_gyro_madgwick_node xnet_gyro_madgwick_gencpp)

#############
## Install ##
#############

# install(TARGETS xnet_gyro_madgwick_node
# )

# all install targets should use catkin DESTINATION variables
# See http://ros.org/doc/api/catkin/html/adv_user_guide/variables.html

## Mark executable scripts (Python etc.) for installation
## in contrast to setup.py, you can choose the destination
# install(PROGRAMS
#   scripts/my_python_script
#   DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

## Mark executables and/or libraries for installation
# install(TARGETS encoder_test encoder_test_node
#   ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

## Mark cpp header files for installation
# install(DIRECTORY include/${PROJECT_NAME}/
#   DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
#   FILES_MATCHING PATTERN "*.h"
#   PATTERN ".svn" EXCLUDE
# )
