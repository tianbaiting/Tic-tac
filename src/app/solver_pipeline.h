#ifndef TICTAC_APP_SOLVER_PIPELINE_H
#define TICTAC_APP_SOLVER_PIPELINE_H

namespace tictac::app {

// [EN] Application-level entry point: orchestrates the full WPCD Faddeev solver
// pipeline (parameters -> PW basis -> fWP -> P123 -> V -> SWP -> on-shell ->
// resolvent -> Faddeev -> output) and returns the process exit code. Solver
// internals live under src/core; this layer only sequences them exactly as the
// documented workflow prescribes, without altering the kernel algorithms.
// / [CN] 应用层入口：编排完整的 WPCD Faddeev 求解流水线
// (参数 -> 分波基 -> fWP -> P123 -> V -> SWP -> on-shell -> resolvent
// -> Faddeev -> 输出) 并返回进程退出码。求解器内核位于 tictac::core；
// 本层仅按文档工作流顺序串接，不改写内核算法。
int run_solver(int argc, char* argv[]);

} // namespace tictac::app

#endif // TICTAC_APP_SOLVER_PIPELINE_H
