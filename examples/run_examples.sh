#!/bin/bash

# Tic-tac 示例运行脚本
# 该脚本演示了如何运行不同类型的计算

set -e  # 遇到错误立即退出

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 打印带颜色的消息
print_message() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_header() {
    echo -e "${BLUE}=== $1 ===${NC}"
}

# 检查程序是否已构建
check_program() {
    if [ ! -f "./tic-tac" ]; then
        print_error "程序未找到，请先运行 'make' 编译程序"
        exit 1
    fi
}

# 生成配置文件
generate_config() {
    print_message "生成默认配置文件..."
    python3 config.py save data/input.txt
}

# 示例1：基本的三核子散射计算
example_basic() {
    print_header "示例1: 基本三核子散射计算"
    
    # 生成配置
    generate_config
    
    # 修改配置为基本计算
    cat > data/input.txt << EOF
# 基本三核子散射计算配置
two_J_3N_max=1
Np_WP=30
Nq_WP=30
J_2N_max=1
tensor_force=true
solve_faddeev=true
solve_dense=false
potential_model=LO_internal
output_folder=output/basic_scattering
EOF
    
    print_message "开始基本散射计算..."
    ./tic-tac --input data/input.txt
    print_message "基本散射计算完成！"
}

# 示例2：高精度计算
example_high_precision() {
    print_header "示例2: 高精度计算"
    
    cat > data/input.txt << EOF
# 高精度计算配置
two_J_3N_max=3
Np_WP=100
Nq_WP=100
J_2N_max=3
tensor_force=true
solve_faddeev=true
solve_dense=false
potential_model=N3LO_Idaho
output_folder=output/high_precision
P123_omp_num_threads=8
EOF
    
    print_warning "高精度计算需要更长时间和更多内存..."
    ./tic-tac --input data/input.txt
    print_message "高精度计算完成！"
}

# 示例3：参数扫描
example_parameter_scan() {
    print_header "示例3: 能量参数扫描"
    
    # 创建能量文件
    cat > data/lab_energies_scan.txt << EOF
# 实验室系能量 (MeV)
1.0
2.0
3.0
5.0
10.0
15.0
20.0
EOF
    
    cat > data/input.txt << EOF
# 参数扫描配置
two_J_3N_max=1
Np_WP=50
Nq_WP=50
J_2N_max=2
tensor_force=true
solve_faddeev=true
potential_model=LO_internal
energy_input_file=data/lab_energies_scan.txt
output_folder=output/parameter_scan
EOF
    
    print_message "开始能量参数扫描..."
    ./tic-tac --input data/input.txt
    print_message "参数扫描完成！"
}

# 示例4：测试运行（快速验证）
example_test() {
    print_header "示例4: 快速测试运行"
    
    cat > data/input.txt << EOF
# 快速测试配置
two_J_3N_max=1
Np_WP=20
Nq_WP=20
J_2N_max=1
tensor_force=false
solve_faddeev=true
solve_dense=false
potential_model=LO_internal
output_folder=output/test_run
production_run=false
EOF
    
    print_message "开始快速测试..."
    ./tic-tac --input data/input.txt
    print_message "测试运行完成！"
}

# 示例5：氘核-质子Ay计算（190 MeV/u）
example_deuteron_ay() {
    print_header "示例5: 氘核-质子分析能力(Ay)计算"
    
    print_message "这是一个专业的极化散射计算示例"
    print_message "计算190 MeV/u氘核与质子散射的分析能力"
    
    if [ ! -f "examples/deuteron_proton_Ay.py" ]; then
        print_error "氘核Ay计算脚本不存在"
        return 1
    fi
    
    print_warning "注意：这是一个高精度计算，可能需要数小时完成"
    echo -e "继续吗？(y/N): \c"
    read response
    
    if [ "$response" = "y" ] || [ "$response" = "Y" ]; then
        print_message "启动氘核-质子Ay计算..."
        python3 examples/deuteron_proton_Ay.py
    else
        print_message "已取消计算"
        print_message "要手动运行，请使用: python3 examples/deuteron_proton_Ay.py"
    fi
}

# 清理输出文件
clean_output() {
    print_header "清理输出文件"
    if [ -d "output" ]; then
        print_warning "删除输出目录..."
        rm -rf output/*
        print_message "输出文件已清理"
    else
        print_message "无需清理，输出目录不存在"
    fi
}

# 显示帮助信息
show_help() {
    echo "Tic-tac 示例运行脚本"
    echo ""
    echo "用法: $0 [选项]"
    echo ""
    echo "选项:"
    echo "  basic      - 运行基本三核子散射计算"
    echo "  precision  - 运行高精度计算"
    echo "  scan       - 运行参数扫描计算"
    echo "  test       - 运行快速测试"
    echo "  ay         - 氘核-质子分析能力计算(190 MeV/u)"
    echo "  clean      - 清理输出文件" 
    echo "  help       - 显示此帮助信息"
    echo ""
    echo "示例:"
    echo "  $0 test        # 快速测试"
    echo "  $0 basic       # 基本计算"
    echo "  $0 precision   # 高精度计算"
    echo "  $0 ay          # 氘核分析能力计算"
}

# 主函数
main() {
    print_header "Tic-tac 三体核物理计算示例"
    
    # 检查程序
    check_program
    
    # 创建必要目录
    mkdir -p data output
    
    # 处理命令行参数
    case "${1:-help}" in
        "basic")
            example_basic
            ;;
        "precision")
            example_high_precision
            ;;
        "scan")
            example_parameter_scan
            ;;
        "test")
            example_test
            ;;
        "ay")
            example_deuteron_ay
            ;;
        "clean")
            clean_output
            ;;
        "help"|*)
            show_help
            ;;
    esac
    
    print_message "脚本执行完成！"
}

# 运行主函数
main "$@"