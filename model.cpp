/*
 *  model.cpp
 *
*  This code computes an optimal set of video representations in an Adaptive Streaming scenario using an ILP model.
 *  The model assumes that a CDN entity prepare multiple representations of the videos in its catalogue. Each representation
 *  is characterized by only 3 parameters: the own video, the encoding rate and the spatial resolution.
 *  The set is optimal in terms of maximizing the sum of individual satisfactions of a given population of users, where
 *  users are defined by their video requests, bandwidth of internet connections and display sizes.
 *  Finally, the CDN entity is constrained by a limited CDN budget to buy bandwidth to ISPs and a limited number of representation
 *  that it can prepare. The ILP model is solved using IBM ILOG CPLEX solver. More details are in [1].
 *
 *  This code loads from a folder named 'scenarios/' a file containing all the input parameters characterizing a given scenario:
 *  video catalogue, encoding settings, satisfaction model, users population, .... These parameters are the input of the model.
 *  This code saves 3 files to 3 folders named 'optimalValues/', 'status/' and 'representations/'. The file saved in 'optimalValues/'
 *  contains the optimal valued found by the CPLEX solver for the decision variables of the model. The file saved in 'status/'
 *  contains the status of solution found by CPLEX solver, for example, if the solution is optimal, the status is the word 'optimal'.
 *  Finally, The file saved in 'representations/' contains a list of the optimal representations found by the solver. Each row in the
 *  file corresponds to a representation, where first, second and third column corresponds to video index, encoding rate (in kbps)
 *  and resolution (e.g. 360p or 1080p).
 *
 *  NOTE: This code uses IBM ILOG CPLEX optimizer [2] to solve the ILP formulation in [1]. Thus, you need to install IBM ILOG CPLEX
 *  and to prepare a makefile according to the IBM ILOG CPLEX indications to be able to compile and run this code.
 *
 *  [1] Laura Toni, Ramon Aparicio-Pardo, Gwendal Simon, Alberto Blanc, Pascal Frossard, "Optimal Set of Video Representations in Adaptive Streaming,"
 *  in Proc. ACM Multimedia Systems conference, MMSys 2014, Singapore, March 2014
 *  [2] http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/
 *
 */

/*#include <vector.h>
 #include <string.h>
 #include <stdio.h>
 #include <stdlib.h>*/

#include <ilcplex/ilocplex.h>
using namespace std;
ILOSTLBEGIN

typedef IloArray<IloArray<IloNumArray> > Float3DArray; // 3D matrix of float coefficients
typedef IloArray<IloNumArray>    FloatMatrix;  // 2D matrix of float coefficients
typedef IloArray<IloArray<IloNumVarArray> > NumVar3DArray; // 3D matrix of variables
typedef IloArray<IloNumVarArray> NumVarMatrix; // 2D matrix of variables

static void
programme(IloNum optimalityGap, IloInt maxNbRepres, IloNum budgetCapacity, IloNum minRatioOfSatisfiedUsers, string scenarioString);

int main(){

    IloNum optimalityGap = 1e-4;
    string scenarioString = "scenario";
    IloInt maxNbRepres = 60;
    IloInt budgetCapacity =1500; // kbps
    IloNum minRatioOfSatisfiedUsers = 0.95;
    programme(optimalityGap, maxNbRepres, budgetCapacity, minRatioOfSatisfiedUsers, scenarioString);
	return 0;
} // END main

static void programme(IloNum optimalityGap, IloInt maxNbRepres, IloNum budgetCapacity, IloNum minRatioOfSatisfiedUsers, string scenarioString)
{
	ILOSTLBEGIN
	IloEnv   env;
	try {
        IloInt  u, l; // user index 'u', video representation index 'l'

        //NOTE: the model considers that a video representation 'l' corresponds to one unique triple (v_l, r_l, s_l),
        // where 'v_l' is the index of the coded video stream, 'r_l' is the encoding rate and 's_l' is the spatial resolution.

        // U: nb Users
        // L: nb Representations
        // V: nb Videos
		// INPUT PARAMETERS 
        FloatMatrix satisfaction_2D(env);//  satisfacction level experienced by user 'u' watching representation 'l'; UxL matrix
        FloatMatrix demand(env);// d_uv: 1 if user 'u' demands video 'v'; 0, otherwise; UxV binary matrix
        FloatMatrix allowedRepresentations_2D(env); // a_ul: 1 if video representation 'l' is compatible with user 'u' features, that is, if user display can support resolution and rate of representation 'l'; UxL matrix
        IloNumArray bitRate(env); //b_l in R+ : encoding rate (kbps) of video representation 'l'; 1xL vector
        IloNumArray repRates_ids(env); //r_l in Z : encoding rate index of video representation 'l'; 1xL vector
        IloNumArray repResolutions(env); //p_l : spatial resolution value (e.g. 1080p) of video representation 'l'; 1xL vector
        IloNumArray repResolutions_ids(env); //s_l in Z : spatial resolution index of video representation 'l'; 1xL vector
        IloNumArray repVideos(env); //v_l in Z: video index of video representation 'l'; 1xL vector
        IloNumArray linkCapacity(env); //c_u in R+ : capacity (kbps) of internet connection of user  'u'; 1xU vector
		
        string scenarioFullPathString = "scenarios/" +  scenarioString + ".dat";
		const char* scenarioFilename = scenarioFullPathString.c_str();
		ifstream infile(scenarioFilename);
		
		if (!infile) {
			cerr << "ERROR: could not open file '" << scenarioFilename
			<< "' for reading" << endl;
			throw(-1);
		}
		infile >> satisfaction_2D >> demand >> allowedRepresentations_2D >> bitRate >> linkCapacity >> repRates_ids >> repResolutions >> repResolutions_ids >> repVideos;

		IloInt nbUsers  = demand.getSize();
		IloInt nbRepres = bitRate.getSize();
		
		cout << "nbUsers " << nbUsers << endl;
        cout << "nbRepres " << nbRepres << endl;
				
		// DECISION VARIABLES 
        // 1) alpha_ul variables: 1 if user 'u' is served by a representation 'l'; 0, otherwise
		NumVarMatrix alpha(env, nbUsers);
		for(u = 0; u < nbUsers; u++){
			alpha[u] = IloNumVarArray(env, nbRepres, 0, 1, ILOINT);
		}
        // 2) beta_l variables:  1 if representation 'l' is played by any user; 0, otherwise
		IloNumVarArray beta(env, nbRepres, 0, 1, ILOINT);
		
        // 3) gamma_u variables: 1 if user 'u'; 0, otherwise
		IloNumVarArray gamma(env, nbUsers, 0, 1, ILOINT);
		
		cout << "var done " << endl;
		
		
		// CONSTRAINTS
		IloModel model(env);
		
        // a) Constraint: link between beta and alpha (lower bound on beta)
		for(u = 0; u < nbUsers; u++)
			for(l = 0; l < nbRepres; l++)
                model.add(alpha[u][l] <= beta[l]); // (u in U, l in L)
		
        // b) Constraint: link between beta and alpha (upper bound on beta)
		for(l = 0; l < nbRepres; l++){
			IloExpr betaUb(env);
			for(u = 0; u < nbUsers; u++)
				betaUb += alpha[u][l]; 
            model.add(beta[l]<=betaUb); // (l in L)
			betaUb.end();
		}
		
        // c) Constraint: link between gamma and alpha (lower bound on gamma)
		for(u = 0; u < nbUsers; u++)
			for(l = 0; l < nbRepres; l++)
                model.add(alpha[u][l] <= gamma[u]); // (u in U, l in L)
		
        // d) Constraint:  link between gamma and alpha (upper bound on gamma)
		for(u = 0; u < nbUsers; u++){
			IloExpr gammaUb(env);
			for(l = 0; l < nbRepres; l++){
				gammaUb += alpha[u][l];
			}
            model.add(gamma[u]<=gammaUb); // (u in U)
		}
        // e) Constraint: user 'u' can only play representations of the requested video 'v_u'
		for(u = 0; u < nbUsers; u++){
			for(l = 0; l < nbRepres; l++){
				int v_l = repVideos[l];
                model.add(alpha[u][l] <= demand[u][v_l]); // (u in U, l in L)
			}
		}		
        // g)  Constraint: user 'u' can only play one video representation
		for(u = 0; u < nbUsers; u++){
			model.add(IloSum(alpha[u]) <= 1); // (u in U)
		}
		
        // h) Constraint: user 'u' cannnot play reprsentations exceeding banwidth of user connection
		IloRangeArray singleLinkCapacities(env);
		for(u = 0; u < nbUsers; u++) {
			IloExpr linkLoad = IloScalProd(bitRate, alpha[u]);
			singleLinkCapacities.add(linkLoad <= linkCapacity[u]);
			model.add(singleLinkCapacities); // u in U
			//model.add(linkLoad <= linkCapacity[u]); // u in U
			linkLoad.end();
		}
        // i) Constraint: CDN cannot exceed its CDN budget in terms of kbps
		IloExpr totalLoad(env);
		for(u = 0; u < nbUsers; u++) {
			totalLoad += IloScalProd(bitRate, alpha[u]);
		}
		model.add(totalLoad <= budgetCapacity*nbUsers); 
		totalLoad.end();
		
        // j) Constraint: CDN does not prepare more than a maximal number of video representations
		model.add(IloSum(beta) <= maxNbRepres); 
		
        // k) Constraint: At least, a fraction of users should be served
		model.add(IloSum(gamma) >= nbUsers*minRatioOfSatisfiedUsers); 
		
        // l) Constraint: a user  'u' only can play represenations with spatial resolutions and rates compatible with the user display
		for(u = 0; u < nbUsers; u++){
			for(l = 0; l < nbRepres; l++){
				model.add(alpha[u][l] <= allowedRepresentations_2D[u][l]); // (u in U, v in V, l in L)
			}
		}
		cout << "ctr done " << endl;
		

        // a) OBJECTIVE : Maximizing the overall user satisfaction.
		IloExpr obj(env);
		for(u = 0; u < nbUsers; u++) {
			obj += IloScalProd(satisfaction_2D[u], alpha[u]);
		}
		model.add(IloMaximize(env, obj));
		obj.end();
		cout << "obj done " << endl;
		
		IloCplex cplex(model);
		cplex.setParam(IloCplex::EpGap,optimalityGap);
		cplex.solve();
		
		cplex.out() << "Solution status " << cplex.getStatus() << endl;
		cplex.out() << "Objective value " << cplex.getObjValue() << endl;
		
		
		// Optimal Values of the decision variables
        // 1) alpha_ul variables
		FloatMatrix opt_alpha_values(env, nbUsers);
		for(u = 0; u < nbUsers; u++){
			opt_alpha_values[u] = IloNumArray(env, nbRepres);
			cplex.getValues(opt_alpha_values[u], alpha[u]);
		}
		// 2) beta_vl variables
		IloNumArray opt_beta_values(env, nbRepres);
		cplex.getValues(opt_beta_values, beta);
		
		// 3) Status & Objective Value
		IloAlgorithm::Status status = cplex.getStatus();
		IloNum opt_obj_value = cplex.getObjValue();
			
		// SAVING RESULTS
		std::stringstream maxNbRepres_char;
		maxNbRepres_char << int(maxNbRepres);
		string maxNbRepres_string (maxNbRepres_char.str());
		
		std::stringstream budgetCapacity_char;
		budgetCapacity_char << int(budgetCapacity);
		string budgetCapacity_string (budgetCapacity_char.str());
		
		std::stringstream minRatioOfSatisfiedUsers_char;
		minRatioOfSatisfiedUsers_char << int(minRatioOfSatisfiedUsers*100);
		string minRatioOfSatisfiedUsers_string (minRatioOfSatisfiedUsers_char.str());		
		
		// Saving in /optimalValues the optimal values of the decision variables
        string optValString = "optimalValues/opt_val_" +  scenarioString + "_"+ minRatioOfSatisfiedUsers_string + "_pct_" + maxNbRepres_string + "_reps_"  + budgetCapacity_string +  "_Kbps.dat";
		const char* optValFilename = optValString.c_str();
		ofstream optValFile(optValFilename);
		if (!optValFile) {
			cerr << "ERROR: could not open file '" << optValFile
			<< "' for writing" << endl;
			throw(-1);
		}
		optValFile << opt_obj_value << opt_alpha_values << opt_beta_values << endl;
		
		// Saving in /status the status of the solution
        string statusString = "status/status_" +  scenarioString + "_" + minRatioOfSatisfiedUsers_string + "_pct_" + maxNbRepres_string + "_reps_"  + budgetCapacity_string +  "_Kbps.dat";
		const char* statusFilename = statusString.c_str();
		ofstream statusFile(statusFilename);
		if (!statusFile) {
			cerr << "ERROR: could not open file '" << statusFile
			<< "' for writing" << endl;
			throw(-1);
		}
		statusFile << status << endl;
		
		// Saving in /representations the representations
        string representationsString = "representations/representations_" +  scenarioString + "_" + minRatioOfSatisfiedUsers_string + "_pct_" + maxNbRepres_string + "_reps_"  + budgetCapacity_string +  "_Kbps.txt";
		const char* representationsFilename = representationsString.c_str();
		ofstream representationsFile(representationsFilename);
		if (!representationsFile) {
			cerr << "ERROR: could not open file '" << representationsFile
			<< "' for writing" << endl;
			throw(-1);
		}
		representationsFile << "videoID" << " " << "bitRate" << " " << "resolutions" << "\n"; //Outputs array to txtFile 
		for(l = 0; l < nbRepres; l++){
			if (opt_beta_values[l] > 0) {
				representationsFile << repVideos[l] << " " << bitRate[l] << " " << repResolutions[l] << "\n"; //Outputs array to txtFile 
			}
		}
	}
	catch (IloException& e) {
		cerr << "Concert exception caught: " << e << endl;
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
	}
	
	env.end();
} 
