/*
 * Tic-tac 三体核物理计算程序 — 薄入口
 *
 * main 只负责把进程控制权交给应用编排层 tictac::app::run_solver，后者在
 * src/app/solver_pipeline.cpp 中按文档化 WPCD 工作流串接求解器内核。
 * 求解器算法本身位于 src/core，本文件不包含任何数值/编排逻辑。
 */

#include "app/solver_pipeline.h"

int main(int argc, char* argv[]){
	return tictac::app::run_solver(argc, argv);
}
