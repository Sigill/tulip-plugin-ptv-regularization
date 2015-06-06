#ifndef REG_PTV_H
#define REG_PTV_H

#include <tulip/TulipPluginHeaders.h>
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
	PLUGININFORMATIONS("Regularization p-TV", "Cyrille FAUCHEUX", "2013/01/17", "Perform a p-TV regularization.", "1.0", "Graph")
	Reg_pTV(tlp::PluginContext* context);
	~Reg_pTV() {}
	bool run();
	bool check(std::string &);

private:
	void exportGraph(const int i);
};

PLUGIN(Reg_pTV)

#endif
