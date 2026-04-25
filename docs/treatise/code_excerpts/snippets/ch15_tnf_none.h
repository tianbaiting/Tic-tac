// ===============================================================
// 抽取自仓库 [current]: src/interactions/three_nucleon_force_none.h
// 行号区段：12..22
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
class three_nucleon_force_none : public three_nucleon_force_model
{
public:
	three_nucleon_force_none() = default;
	~three_nucleon_force_none() override = default;

	bool enabled() const override { return false; }
	std::string name() const override { return "none"; }
};

#endif // THREE_NUCLEON_FORCE_NONE_H
