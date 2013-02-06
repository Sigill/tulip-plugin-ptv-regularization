#ifndef _CV_TA_H
#define _CV_TA_H

#include <tulip/TulipPlugin.h>
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
	Reg_pTV_Vector(const tlp::AlgorithmContext& context);
	~Reg_pTV_Vector() {}
	bool run();
	bool check(std::string &);

private:
	void exportGraph(const int i);
};

#endif
