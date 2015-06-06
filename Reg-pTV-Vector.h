#ifndef REG_PTV_VECTOR_H
#define REG_PTV_VECTOR_H

#include <tulip/TulipPluginHeaders.h>
#include <utility>

class Reg_pTV_Vector:public tlp::Algorithm {
private:
	unsigned int iter_max;
	tlp::DoubleProperty *wl;
	tlp::DoubleVectorProperty *f0, *fn_result;
	bool is_isotropic, is_anisotropic;
	double p, q, epsilon;
	int f0_size;

	unsigned int export_interval;
	std::string export_directory;

public:
	PLUGININFORMATIONS("Regularization p-TV (Vector version)","Cyrille FAUCHEUX","2013/01/17","Perform a p-TV regularization (Vector version).","1.0","Graph")
	Reg_pTV_Vector(tlp::PluginContext* context);
	~Reg_pTV_Vector() {}
	bool run();
	bool check(std::string &);

private:
	void exportGraph(const int i);
};

PLUGIN(Reg_pTV_Vector)

#endif
