#ifndef REG_PTV_H
#define REG_PTV_H

#include <tulip/TulipPlugin.h>
#include <utility>

class Reg_pTV:public tlp::Algorithm {
private:
	unsigned int iter_max;
	tlp::DoubleProperty *f0, *wl, *fn_result;
	bool is_isotropic, is_anisotropic;
	double p, q, epsilon;

	unsigned int export_interval;
	std::string export_directory;

public:
	Reg_pTV(const tlp::AlgorithmContext& context);
	~Reg_pTV() {}
	bool run();
	bool check(std::string &);

private:
	void exportGraph(const int i);
};

#endif
