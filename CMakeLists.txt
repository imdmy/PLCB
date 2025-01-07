cmake_minimum_required(VERSION 3.15)

# 项目名称和语言
project(BoundingboxCalculation LANGUAGES CXX)

# 指定 C++ 标准
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 设置解决方案配置
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE)

# 添加可执行文件
add_executable(BoundingboxCalculation ${CMAKE_SOURCE_DIR}/src/BoundingboxCal_main.cpp ${CMAKE_SOURCE_DIR}/src/Boundingbox.cpp)

#add_executable(BoundingboxCalculation ${CMAKE_SOURCE_DIR}/src/ConsoleApplication.cpp)


# 为目标添加头文件目录
target_include_directories(BoundingboxCalculation PRIVATE ${CMAKE_SOURCE_DIR}/include)
