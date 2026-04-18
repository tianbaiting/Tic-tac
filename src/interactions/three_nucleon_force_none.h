#ifndef THREE_NUCLEON_FORCE_NONE_H
#define THREE_NUCLEON_FORCE_NONE_H

#include <string>

#include "three_nucleon_force_model.h"

// [EN] Null-object 3NF implementation. Always disabled; the kernel assembler's `if (tnf.enabled())` test
// evaluates to false, so the 3NF contribution path is never entered and the solver runs through the
// original 2NF-only hot path. / [CN] 3NF 的 null 对象实现；永远处于关闭状态，内核组装器的
// `if (tnf.enabled())` 判断返回 false，因此 3NF 贡献路径永远不会被触发，求解器走的是原来的 2NF-only 热路径。
class three_nucleon_force_none : public three_nucleon_force_model
{
public:
	three_nucleon_force_none() = default;
	~three_nucleon_force_none() override = default;

	bool enabled() const override { return false; }
	std::string name() const override { return "none"; }
};

#endif // THREE_NUCLEON_FORCE_NONE_H
