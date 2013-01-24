#include "Reg-pTV.h"
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cfloat>
#include <stdexcept>

#include <QFileInfo>
#include <QDir>

ALGORITHMPLUGIN(Reg_pTV, "Regularisation p-TV", "Cyrille FAUCHEUX", "17-01-2012", "Alpha", "1.0");

using namespace std;
using namespace tlp;

template <class T>
std::string to_string(T t, std::ios_base & (*f)(std::ios_base&))
{
	std::ostringstream oss;
	oss << f << t;
	return oss.str();
}

std::string random_string(const size_t len);

namespace {
	const char * paramHelp[] = {
		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "DoubleProperty" ) \
			HTML_HELP_DEF( "default", "f0" ) \
			HTML_HELP_BODY() \
			"The property to regularize." \
			HTML_HELP_CLOSE(),

		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "DoubleProperty" ) \
			HTML_HELP_DEF( "default", "fn" ) \
			HTML_HELP_BODY() \
			"The regularized property (result of the regularization)." \
			HTML_HELP_CLOSE(),

		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "double" ) \
			HTML_HELP_DEF( "default", "p" ) \
			HTML_HELP_BODY() \
			"Penalization coefficient (penalizes high variations). Choose p = q for anisotropic model." \
			HTML_HELP_CLOSE(),

		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "double" ) \
			HTML_HELP_DEF( "default", "q" ) \
			HTML_HELP_BODY() \
			"Norm power (choose 2 for isotropic model)." \
			HTML_HELP_CLOSE(),

		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "Unsigned int" ) \
			HTML_HELP_DEF( "default", "1000" ) \
			HTML_HELP_BODY() \
			"The number of iterations to perform." \
			HTML_HELP_CLOSE(),

		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "DoubleProperty" ) \
			HTML_HELP_DEF( "default", "weight" ) \
			HTML_HELP_BODY() \
			"The property holding the weight associated to each edge and node." \
			HTML_HELP_CLOSE(),

		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "Unsigned int" ) \
			HTML_HELP_DEF( "default", "0" ) \
			HTML_HELP_BODY() \
			"Export interval. (0 to disable intermediate export)" \
			HTML_HELP_CLOSE(),

		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "Directory pathname" ) \
			HTML_HELP_BODY() \
			"Export directory." \
			HTML_HELP_CLOSE(),

		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "double" ) \
			HTML_HELP_BODY() \
			"Arbitrary small positive value used when computing gradients in order to avoid 0 values." \
			HTML_HELP_CLOSE(),
	};
}

//======================================================
Reg_pTV::Reg_pTV(const tlp::AlgorithmContext &context):Algorithm(context), f0(NULL), wl(NULL), fn_result(NULL) {
	addParameter< DoubleProperty >      ("f0",                    paramHelp[0]);
	addParameter< DoubleProperty >      ("fn",                    paramHelp[1]);
	addParameter< double >              ("p",                     paramHelp[2], "2");
	addParameter< double >              ("q",                     paramHelp[3], "2");
	addParameter< unsigned int >        ("number of iterations",  paramHelp[4], "1000");
	addParameter< DoubleProperty >      ("weight/lambda",         paramHelp[5]);
	addParameter< double >              ("epsilon",               paramHelp[8], "0.01");
	addParameter< unsigned int >        ("export interval",       paramHelp[6], "0");
	addParameter< string >              ("dir::export directory", paramHelp[7], "", false);
}

#define CHECK_PROP_PROVIDED(PROP, STOR) \
	do { \
		if(!dataSet->get(PROP, STOR)) \
			throw std::runtime_error(std::string("No \"") + PROP + "\" provided."); \
	} while(0)

bool Reg_pTV::check(std::string &err)
{
	try {
		if(dataSet == NULL)
			throw std::runtime_error("No dataset provided.");

		CHECK_PROP_PROVIDED("f0", this->f0);

		CHECK_PROP_PROVIDED("fn", this->fn_result);

		CHECK_PROP_PROVIDED("p", this->p);
		
		CHECK_PROP_PROVIDED("q", this->q);

		CHECK_PROP_PROVIDED("number of iterations", this->iter_max);

		CHECK_PROP_PROVIDED("weight/lambda", this->wl);

		CHECK_PROP_PROVIDED("export interval", this->export_interval);

		CHECK_PROP_PROVIDED("epsilon", this->epsilon);

		if(this->epsilon <= 0)
			throw std::runtime_error("epsilon must be grater than 0.");

		// No need for an export directory if export is disabled
		if(this->export_interval > 0) {
			CHECK_PROP_PROVIDED("dir::export directory", this->export_directory);

			// Checking if we can write in the export directory
			QString qstring_export_directory = QString::fromStdString(this->export_directory);
			QFileInfo info_export_directory(qstring_export_directory);
			QDir qdir_export_directory(qstring_export_directory);

			if(info_export_directory.exists()) {
				if(info_export_directory.isDir()) {
					if(qdir_export_directory.entryInfoList(QDir::NoDotAndDotDot | QDir::AllEntries).count() != 0) {
						throw std::runtime_error("Export directory (" + qdir_export_directory.absolutePath().toStdString() + ") is not empty.");
					}
				} else {
					throw std::runtime_error("Export directory (" + qdir_export_directory.absolutePath().toStdString() + ") already exists but is not a directory.");
				}
			} else {
				if(!qdir_export_directory.mkpath(".")) {
					throw std::runtime_error("Export directory (" + qdir_export_directory.absolutePath().toStdString() + ") cannot be created.");
				}
			}
		}

		if(this->iter_max < 0) {
			std::ostringstream m;
			m << "Invalid number of iterations: " << this->iter_max;
			throw std::runtime_error(m.str());
		}

		if( (this->p <= 0) || (this->q <= 0) )
			throw std::runtime_error("p and q must be greater than 0.");

		this->is_isotropic = (this->q == 2.0);
		this->is_anisotropic = (this->p == this->q);

		if(!this->is_isotropic && !this->is_anisotropic)
			throw std::runtime_error("You need to either choose (p > 0 and q = 2) or (p = q > 0).");

		if(this->is_anisotropic && !this->is_isotropic) { // Anisotropic
			std::cout << "Running the anisotropic model." << std::endl;
		} else if(this->is_isotropic && !this->is_anisotropic) { // Isotropic
			std::cout << "Running the isotropic model." << std::endl;
		} else {
			std::cout << "Running both isotropic and anisotropic models." << std::endl;
		}

		std::cout << "Processing graph \"" << graph->getName() <<  "\"" << std::endl;
		std::cout << "p: " << p << std::endl;
		std::cout << "q: " << q << std::endl;
		std::cout << "Number of iterations: " << this->iter_max << std::endl;
		std::cout << "epsilon: " << epsilon << std::endl;
		std::cout << "Export interval: " << this->export_interval << std::endl;
		std::cout << "Export directory: " << this->export_directory << std::endl;
	} catch (std::runtime_error &ex) {
		err.assign(ex.what());
		return false;
	}

	return true;
}

//======================================================
bool Reg_pTV::run() {
	{
		DoubleProperty *tmp, *fn = NULL, *fnp1 = NULL, *grad = NULL;
		Iterator<node> *itNodesU;
		Iterator<edge> *itEdges;
		node u, v;
		edge e;
		double num, denum, b, grad_u, grad_v, lambda_u, fn_v;
		bool continueProcess = true;
		std::string fn_name, fnp1_name, grad_name;

		if(pluginProgress)
			pluginProgress->setComment("Initializing");

		/*
		 * Find an unused property name for fn & fnp1
		 */
		do {
			fn_name = random_string(6);
		} while(graph->existLocalProperty(fn_name));
		do {
			fnp1_name = random_string(6);
		} while(graph->existLocalProperty(fnp1_name));

		fn = graph->getLocalProperty< DoubleProperty >(fn_name);
		fnp1 = graph->getLocalProperty< DoubleProperty >(fnp1_name);

		/*
		 * Find an unused property name for grad
		 * Only if we are using the isotropic model
		 */
		if(this->is_isotropic && !this->is_anisotropic) {
			do {
				grad_name = random_string(6);
			} while(graph->existLocalProperty(grad_name));

			grad = graph->getLocalProperty< DoubleProperty >(grad_name);
		}

		/*
		 * Initializing the function to regularize
		 */
		fn->copy(this->f0);

		if(pluginProgress)
			pluginProgress->setComment("Processing");

		/*
		 * Main loop
		 */
		for(unsigned int i = 0; i < iter_max; ++i) {
			if(this->is_anisotropic && !this->is_isotropic) {
				/*
				 * If using the anisotropic model
				 */

				itNodesU = graph->getNodes();
				while(itNodesU->hasNext()) {
					u = itNodesU->next();

					lambda_u = this->wl->getNodeValue(u);
					num = lambda_u * this->f0->getNodeValue(u);
					denum = lambda_u;

					itEdges = graph->getInOutEdges(u);
					while(itEdges->hasNext()) {
						e = itEdges->next();
						v = graph->opposite(e, u);

						fn_v = fn->getNodeValue(v);

						/*
						 * fabs() is regularized to avoir 0
						 */
						b = pow(this->wl->getEdgeValue(e), this->p / 2) * pow(fabs(fn->getNodeValue(u) - fn_v) + this->epsilon, this->p - 2);
						num += b * fn_v;
						denum += b;
					}
					delete itEdges;

					fnp1->setNodeValue(u, num / denum);
				}
				delete itNodesU;

			} else if(this->is_isotropic && !this->is_anisotropic) {
				/*
				 * If using the anisotropic model
				 *
				 * Precomputing gradient norms
				 * Norm q of a weighted gradient of the function f at node u
				 * [ sum_{v~u} w_{uv}^{q/2} . |f(v) - f(u)|^q ]^{1/q}
				 * We are computing it for q = 2, and at power p-2
				 */
				{
					double epsilon_squared = this->epsilon * this->epsilon,
					       gradient_power = (this->p - 2) / 2;

					itNodesU = graph->getNodes();
					while(itNodesU->hasNext()) {
						u = itNodesU->next();

						grad_u = 0;

						itEdges = graph->getInOutEdges(u);
						while(itEdges->hasNext()) {
							e = itEdges->next();
							v = graph->opposite(e, u);
							grad_u += this->wl->getEdgeValue(e) * pow(fabs(fn->getNodeValue(u) - fn->getNodeValue(v)), 2);
						}
						delete itEdges;

						// grad_u is the gradient^2

						/*
						 * Regularized to avoid divisions by 0
						 * Computing it at power p-2
						 */
						grad->setNodeValue(u, pow(grad_u + epsilon_squared, gradient_power)); 
					}
					delete itNodesU;
				}

				itNodesU = graph->getNodes();
				while(itNodesU->hasNext()) {
					u = itNodesU->next();

					lambda_u = this->wl->getNodeValue(u);
					num = lambda_u * this->f0->getNodeValue(u);
					denum = lambda_u;

					itEdges = graph->getInOutEdges(u);
					while(itEdges->hasNext()) {
						e = itEdges->next();
						v = graph->opposite(e, u);

						grad_u = grad->getNodeValue(u);
						grad_v = grad->getNodeValue(v);

						b = this->wl->getEdgeValue(e) * (grad_u + grad_v);
						num += b * fn->getNodeValue(v);
						denum += b;
					}
					delete itEdges;

					fnp1->setNodeValue(u, num / denum);
				}
				delete itNodesU;

			} else {
				/*
				 * Isotropic & anisotropic model
				 * p = q = 2
				 * Some computations are simpliified
				 */

				itNodesU = graph->getNodes();
				while(itNodesU->hasNext()) {
					u = itNodesU->next();

					lambda_u = this->wl->getNodeValue(u);
					num = lambda_u * this->f0->getNodeValue(u);
					denum = lambda_u;

					itEdges = graph->getInOutEdges(u);
					while(itEdges->hasNext()) {
						e = itEdges->next();
						v = graph->opposite(e, u);

						fn_v = fn->getNodeValue(v);

						b = this->wl->getEdgeValue(e);
						num += b * fn_v;
						denum += b;
					}
					delete itEdges;

					fnp1->setNodeValue(u, num / denum);
				}
				delete itNodesU;
			}

			tmp = fn;
			fn = fnp1;
			fnp1 = tmp;

			if(pluginProgress) {
				if(pluginProgress->state() != TLP_CONTINUE)
					continueProcess = false;

				if((i + 1) % 10 == 0) {
					pluginProgress->progress(i+1, iter_max);
					std::cerr << "Iteration " << (i+1) << "/" << iter_max << std::endl;
				}

				if((this->export_interval > 0) && ((i + 1) % this->export_interval) == 0) {
					try {
						exportGraph(i+1);
					} catch (std::runtime_error &ex) {
						std::cerr << "Export failed et iteration #" << i+1 << ": " << ex.what() << std::endl;
					}
				}
			}

			if(!continueProcess)
				break;
		} // End main loop

		this->fn_result->copy(fn);

		graph->delLocalProperty(fn_name);
		graph->delLocalProperty(fnp1_name);
		if(this->is_isotropic && !this->is_anisotropic)
			graph->delLocalProperty(grad_name);
	}

	if(pluginProgress && iter_max > 0)
		pluginProgress->progress(iter_max, iter_max);

	pluginProgress->setComment("Computing selection");

	if( this->export_interval > 0 )
	{
		try {
			exportGraph(iter_max);
		} catch (std::runtime_error &ex) {
			std::cerr << "Export failed et iteration #" << iter_max << ": " << ex.what() << std::endl;
		}
	}

	return true;
}
//=======================================================================

void Reg_pTV::exportGraph(const int i) {
	std::ostringstream directory_name;
	directory_name << this->export_directory << "/" << std::setfill('0') << std::setw(6) << i << ".tlp";

	if(!tlp::saveGraph(graph, directory_name.str()))
		throw std::runtime_error(directory_name.str() + " cannot be written");
}

std::string random_string(const size_t len) {
	static const char alphanum[] =
		"0123456789"
		"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		"abcdefghijklmnopqrstuvwxyz";

	std::string result(len, '-');

	for (size_t i = 0; i < len; ++i) {
		result[i] = alphanum[rand() % 61];
	}

	return result;
}
