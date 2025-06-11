cmake_minimum_required(VERSION 3.1)
project(stap++)

set(CMAKE_CXX_STANDARD 11)

option(STAP++_DEBUG "Output detailed results for debugging." OFF)
if(STAP++_DEBUG)
    add_definitions(-D_DEBUG_)
endif()

# === 添加 Eigen 支持 ===
# 假设你系统中已安装 Eigen（例如 /usr/include/eigen3）
# 或你将 Eigen 解压到项目目录的 "external/eigen" 中

# 方法 1：系统安装（推荐）
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

# 方法 2：使用自带 Eigen（如果没有系统安装）
# include_directories(external/eigen) # 取消注释以使用项目内的 Eigen 头文件路径

# 项目源代码目录
include_directories(h)

# 源代码和头文件
aux_source_directory(cpp SRC)
file(GLOB_RECURSE HEAD h/*.h)

# 分组（用于 Visual Studio 等 IDE）
source_group("Source Files" FILES ${SRC})
source_group("Header Files" FILES ${HEAD})

# 生成可执行文件
add_executable(stap++ ${SRC} ${HEAD})
