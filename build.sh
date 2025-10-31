#!/bin/bash
# Tic-tac 构建脚本

set -e

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 显示帮助信息
show_help() {
    echo "Tic-tac 构建脚本"
    echo ""
    echo "用法: ./build.sh [选项]"
    echo ""
    echo "选项:"
    echo "  clean       - 清理构建目录"
    echo "  debug       - 调试模式构建"
    echo "  release     - 发布模式构建（默认）"
    echo "  rebuild     - 清理后重新构建"
    echo "  help        - 显示此帮助信息"
    echo ""
}

# 清理构建目录
clean_build() {
    echo -e "${YELLOW}清理构建目录...${NC}"
    rm -rf build
    rm -f tic-tac
    echo -e "${GREEN}清理完成${NC}"
}

# 构建项目
build_project() {
    local build_type=$1
    
    echo -e "${GREEN}开始构建 Tic-tac (${build_type} 模式)...${NC}"
    
    # 创建构建目录
    mkdir -p build
    cd build
    
    # 运行CMake配置
    echo -e "${YELLOW}运行 CMake 配置...${NC}"
    cmake -DCMAKE_BUILD_TYPE=${build_type} ..
    
    # 编译
    echo -e "${YELLOW}编译中...${NC}"
    make -j$(nproc)
    
    # 复制可执行文件到根目录
    if [ -f bin/tic-tac ]; then
        cp bin/tic-tac ..
        echo -e "${GREEN}可执行文件已复制到根目录${NC}"
    fi
    
    cd ..
    
    echo -e "${GREEN}构建成功！${NC}"
    echo -e "${GREEN}可执行文件: ./tic-tac${NC}"
}

# 主逻辑
case "${1:-release}" in
    clean)
        clean_build
        ;;
    debug)
        build_project "Debug"
        ;;
    release)
        build_project "Release"
        ;;
    rebuild)
        clean_build
        build_project "Release"
        ;;
    help|--help|-h)
        show_help
        ;;
    *)
        echo -e "${RED}未知选项: $1${NC}"
        show_help
        exit 1
        ;;
esac
