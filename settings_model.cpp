/*
 *  This code generates a testing scenario file for the ILP model in model.cpp using the default settings in Section 6 of [1].
 *  The scenario file contains an instance of values for the input parameters of the ILP model described in [2].
 *  This file is saved in a folder named 'scenarios/'. The objective of this code is merely to illustrate how to
 *  assign values to the input parameters of model.cpp to build an instance of the optimisation problem.
 *
 *  NOTE: This code assumes that you will use the scenario file to feed model.cpp code. model.cpp solves
 *  the ILP formulation in [2] by using IBM ILOG CPLEX Optimizer [3]. This CPLEX solver makes use of the ILOG
 *  Concert Technology, that offers a C++ library of classes and functions enabling to define models for
 *  optimization problems and to apply algorithms to those models. Thus, you need to install IBM ILOG CPLEX
 *  and to prepare a makefile according to the IBM ILOG CPLEX indications to be able to compile and run this code.
 *
 *
 *  [1] Laura Toni, Ramon Aparicio-Pardo, Gwendal Simon, Alberto Blanc, Pascal Frossard, "Optimal Set of Video
 *  Representations in Adaptive Streaming," in Proc. ACM Multimedia Systems conference, MMSys 2014, Singapore, March 2014
 *  [2] Laura Toni, Ramon Aparicio-Pardo, Karine Pires, Alberto Blanc, Gwendal Simon, Pascal Frossard, ìOptimal Selection
 *  of Adaptive Streaming Representations,î ACM Transactions on Multimedia Computing, Communications and Applications
 *  (ACM TOMM), vol. 11, no. 2s, article 43, February 2015
 *  [3] http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/
 *
 */


#include <stdlib.h>
#include <vector.h>
#include <iostream.h>
#include <fstream.h>
#include <string.h>
#include <math.h> 
#include <time.h> 
//#include <sstream.h>
//using  std::vector;

#include <ilcplex/ilocplex.h>
using namespace std;

typedef struct resolution {
	int resolutionValue;
	double minBitRate;
	double maxBitRate;
};

typedef struct representation {
	string videoName;
	int videoIndex;
	int resolution;
	int resolutionIndex;
	double bitRate;
	int bitRateIndex;
};

typedef struct user {
	int videoIndex;
	int connectionIndex;
	int deviceIndex;
	float capacity;
	int dispResolIndex;
};

typedef struct userDevice {
//	string deviceName;
	int distribution;
};

typedef IloArray<IloNumArray>    FloatMatrix;  // 2D matrix of float coefficients
typedef IloArray<IloArray<IloNumArray> > Float3DArray; // 3D matrix of float coefficients

static void 
write(string scenarioString,  int HDTVrate, int sportRate, int resolSwitching);


int main(){
	int resolSwitching = 1;
    string scenarioString = "scenario";
    int HDTVrate = -1;
    int sportRate = -1;
    write( scenarioString,  HDTVrate,  sportRate, resolSwitching);
	return 0;
} // END main



static void write(string scenarioString, int HDTVrate, int sportRate, int resolSwitching)
{
	IloEnv   env;
	int nbResol = 4;
	int nbRates = 17;
	int nbDevices = 4;
	int nbConnections =5;
	int nbUsers = 500;
	int nbVideos = 4;
	
	int r,s,l,c,d,u, v;
	
	// Definition of the Rate-distortion functions parameters 
	int resolvalues[]  = { 224,  360,  720, 1080};
	//	RDfunctionParamStruct RDfunctionParamSet [nbVideos][nbResol][nbResol]; //v x display s x s	
	string videoNames[] = {"cartoon", "documentary", "sport" , "movie"}; 
	//			       big_buck_bunny  SnowMnt  RushFieldCuts  old_town_cross

    // parameters Fitting (video v x display resolution d x represenation resolution s)
    float RDparam_minRate[4][4][4]; //(v,d,s)
    float RDparam_maxRate[4][4][4]; //(v,d,s)
    float RDparam_a[4][4][4]  = { //(v,d,s)
		{
			{-0.024564, 0.106675,	-1,		-1},   // disp 224p
			{0.046629, -0.023172, 0.0444,	-1},   // disp 360p
			{-1,	    0.088225, -0.031672, -0.003387}, // disp 720p
			{-1,		 -1,	 -0.067402, -0.010134}  // disp 1080p
		},	
		{
			{-0.013738, 0.000497,	-1,		-1},   // disp 224p 
			{0.095128, -0.019052, 0.006332,	-1},   // disp 360p
			{-1,	   0.038282, -0.017795, 0.00246}, // disp 720p
			{-1,		-1,	     -0.045269, -0.034516}  // disp 1080p
		},
		{
			{-0.096091, -0.040061,	-1,			-1},   // disp 224p 
			{0.040243, -0.12335, -0.062405,	-1},   // disp 360p
			{		-1,	 0.06258, -0.103019, -0.029194}, // disp 720p
			{		-1,	    -1,	   -0.031548, -0.074582}  // disp 1080p			
		},
		{
			{-0.036885, 0.017509,	-1,		-1},   // disp 224p 
			{0.070397, -0.042149, -0.04338,	-1},   // disp 360p
			{-1,		 0.092607, -0.013527, 0.03745}, // disp 720p
			{-1,			-1,		-0.038311, 0.018633}  // disp 1080p
		}
	};
    float RDparam_b[4][4][4]  = { //(v,d,s)
		{
			{35.600774, 2.044954,	 -1,		-1},   // disp 224p 
			{14.456387, 49.199101, 23.969837,	-1},   // disp 360p
			{-1,		23.369155, 166.446371, 80.939813}, // disp 720p
			{-1,			-1,		511.036587, 127.785575}  // disp 1080p
		},
		{
			{19.500928, 21.320402,	-1,		-1},   // disp 224p 
			{25.487623, 52.518504, 74.370257,	-1},   // disp 360p
			{-1,	106.184238, 187.436979, 204.11794}, // disp 720p
			{-1,	-1,			414.674132, 372.063754}  // disp 1080p
		},	
		{
			{188.633385, 167.484789,	-1,		    -1},   // disp 224p 
			{219.794369, 445.593377, 339.127191,	-1},   // disp 360p
			{-1,	     447.377634, 1348.641838, 852.278741}, // disp 720p
			{-1,			-1,		 1137.040541, 1548.175792}  // disp 1080p
		},
		{
			{77.868807, 65.493685,	-1,		-1},   // disp 224p 
			{112.799801, 136.259005, 462.165884,	-1},   // disp 360p
			{-1,		226.491157, 119.494934, 148.76041}, // disp 720p
			{-1,			-1,		270.344906, 148.378612}  // disp 1080p
		}
	};
    float RDparam_c[4][4][4]  = { //(v,d,s)
		{
			{31.63, -87.7,	-1,		-1},   // disp 224p 
			{-60.65, 116.24, -800.08,	-1},   // disp 360p
			{-1,	22.26, -65.56, -1156.78}, // disp 720p
			{-1,	-1,		 1834.94, -523.06}  // disp 1080p	
		},
		{
			{-68.49, -120.68,	-1,		-1},   // disp 224p 
			{-55.62, -105.32, -371.8,	-1},   // disp 360p
			{-1,	  89.47, -74.22, -636.24}, // disp 720p
			{-1,	    -1,		704.83, -165.76}  // disp 1080p				
		},	
		{
			{196.92, 62.29,	-1,		-1},   // disp 224p 
			{235.89, 422.25, -164.01,	-1},   // disp 360p
			{-1,	 426.25, 1574.48, 262.06}, // disp 720p	
			{-1,	-1,			1025.2, 1286.62}  // disp 1080p
		},
		{
			{150.03, 86,	-1,		-1},   // disp 224p 
			{243.34, 259.1, 4214.38,	-1},   // disp 360p
			{-1,	477.13, -543.77, -288.9}, // disp 720p
			{-1,	-1,		-61.45, -1498.73}  // disp 1080p
		}
	};
	
	// Definition of the Resolution Set \mathcal{S}
	resolution resolutionsSet [nbResol]; 
	for (s=0; s<nbResol; s++){
		resolutionsSet[s].resolutionValue = resolvalues[s];
	}	
	// Definition of the QoE Set \mathcal{L}	
	double lowerQoE = 0.6;
	double higherQoE = 1;
	double QoEStep = (higherQoE-lowerQoE)/(nbRates-1);
	
	double QoEset[nbRates];
	double val = lowerQoE; 
	for (r=0; r<nbRates; r++){
		QoEset[r] = val; 
		val += QoEStep;
		cout << QoEset[r] << "\n";
	}
	// Definition of representation_set(s_disp, v) and allowedRepresentations_aux(s_disp, v, l)
	// Definition of the Representation Set \mathcal{L}
	vector <representation> representationsSet;	
	for(v = 0; v < nbVideos; v++){
		for (s=0; s<nbResol; s++){
			float maxRate_this_vs = 0;
			float minRate_this_vs = INFINITY;
			for (r=0; r<nbRates; r++){
				
				representation repres;
				
				repres.videoName = videoNames[v];
				repres.videoIndex = v;
				
				float a = RDparam_a[v][s][s];
				float b = RDparam_b[v][s][s];
				float c = RDparam_c[v][s][s];
				
				repres.bitRate = -c+(b/(1-QoEset[r]-a));
				repres.bitRateIndex = r;
				
				repres.resolution = resolutionsSet[s].resolutionValue;
				repres.resolutionIndex = s;
				
				if (repres.bitRate > 50){
					representationsSet.push_back (repres);
					cout << "v: " << v << " s: " <<  repres.resolution  << "p "<< repres.bitRate << "\n";
					
					if (repres.bitRate >= maxRate_this_vs) {
						maxRate_this_vs = repres.bitRate;
					}
					if (repres.bitRate <= minRate_this_vs ) {
						minRate_this_vs = repres.bitRate;
					}
				}
			}
			IloInt s_max = IloMin(s+1,nbResol);
			IloInt s_min = IloMax(s-1,0);
			for (IloInt s_disp = s_min; s_disp <= s_max; s_disp++) {
				RDparam_minRate[v][s_disp][s] = minRate_this_vs;
				RDparam_maxRate[v][s_disp][s] = maxRate_this_vs;
			}
		}		
	}		
	int nbRepres = representationsSet.size();
	cout << "nbRepres " << nbRepres << endl;
		
	// Definition of the User Devices Set \mathcal{D}
	userDevice userDevicesSet [nbDevices];
	
	string devices[]  = {"smartphone","tablet", "laptop", "HDTV"};
	int devicesDistributions [nbDevices]; 
	
	int otherDevsPopulationRate;
	int smartpPopulationRate;
	int HDTVpopulationRate;
    if (HDTVrate < 0){
		otherDevsPopulationRate = 25;
		smartpPopulationRate = 25;
		HDTVpopulationRate = 25;
	}
	else if (HDTVrate > 80) {
		cerr << "ERROR: HDTVrate is larger than 70 %. \n";
		throw(-1);
	}
	else {
		otherDevsPopulationRate = 10;
		smartpPopulationRate = 80 - HDTVrate;
		HDTVpopulationRate = HDTVrate;
	}
	
	devicesDistributions[0] = smartpPopulationRate;
	cout << devices[0].c_str() << " " << devicesDistributions[0] << "\n";
	for (d=1; d<nbDevices; d++){
		int populationRate;
		if (d==3){
			populationRate = HDTVpopulationRate;
		}
		else {
			populationRate = otherDevsPopulationRate;
		}
		devicesDistributions[d] = devicesDistributions[d-1] + populationRate;
		
		cout << devices[d].c_str() << " " << devicesDistributions[d] << "\n";
	}
	
	cout  << "\n";
	
	for (d=0; d<nbDevices; d++){
//		userDevicesSet[d].deviceName = devices[d];
		userDevicesSet[d].distribution = devicesDistributions[d];
	}
	
	// Definition of the User Connections Set \mathcal{C}
	struct userConnections {
	//	string connectionName;
		int minBandwidth;
		int maxBandwidth;
		int distribution;
	} userConnectionsSet [nbConnections];
	
	string connections[] = {"ADSL-slow", "ADSL-fast", "FTTH",  "Wifi", "3G"};
	int connecMinBandwidths[] = { 300,       700,      1500,      400,      150};
	int connecMaxBandwidths[] = {3000,     10000,     25000,     4000,      800};
	int connecDistributions[] = {  10,        30,        10,       30,       20};
	for (c=0; c<nbConnections; c++){
//		userConnectionsSet[c].connectionName = connections[c];
		userConnectionsSet[c].minBandwidth   = connecMinBandwidths[c];
		userConnectionsSet[c].maxBandwidth   = connecMaxBandwidths[c];
		userConnectionsSet[c].distribution   = connecDistributions[c];			
	}			
	
	// Definition of the User Set \mathcal{U}	
	user usersSet [nbUsers];
	
	int videosDistributions [nbVideos]; 
	
	int otherVidsPopularityRate;
	int cartoonPopularityRate;
	int sportPopularityRate;
	if (sportRate < 0){
		otherVidsPopularityRate = 25;
		cartoonPopularityRate = 25;
		sportPopularityRate = 25;
	}
	else if (sportRate > 80) {
		cerr << "ERROR: sportRate is larger than 80 %. \n";
		throw(-1);
	}
	else {
		otherVidsPopularityRate = 10;
		cartoonPopularityRate = 80 - sportRate;
		sportPopularityRate = sportRate;
	}
	
	videosDistributions[0] = cartoonPopularityRate;
	cout << videoNames[0].c_str() << " " << videosDistributions[0] << "\n";
	for (v=1; v<nbVideos; v++){
		int videoRate;
		if (v==2){
			videoRate = sportPopularityRate;
		}
		else {
			videoRate = otherVidsPopularityRate;
		}
		videosDistributions[v] = videosDistributions[v-1] + videoRate;
		
		cout << videoNames[v].c_str() << " " << videosDistributions[v] << "\n";
	}
	cout  << "\n";
	
	IloNumArray videosIndeces(env,nbUsers);
	IloNumArray connectionsIndeces(env,nbUsers);
	IloNumArray devicesIndeces(env,nbUsers);
	IloNumArray userMaxRateDueUserResol(env,nbUsers);
	IloNumArray userMinRateDueUserResol(env,nbUsers);
	IloNumArray userMaxRate(env,nbUsers); 	
	
	srand (time(NULL));
	
	for (u=0; u<nbUsers; u++){
		
		int randVideo = rand() % 100 + 1;
		for (v=0; v<nbVideos; v++){
			if (randVideo <= videosDistributions[v]) {
				videosIndeces[u] = int(v);
				usersSet[u].videoIndex = v;
				break;
			}
		}
		
		// device assignment --> allowed resolutions assignment
		int randDevice = rand() % 100 + 1;
		for (d=0; d<nbDevices; d++){
			if (randDevice <= userDevicesSet[d].distribution) {
				devicesIndeces[u] = d;
				usersSet[u].deviceIndex = d;
				usersSet[u].dispResolIndex =d;  
				break;
			}
		}
		
		int s_disp = usersSet[u].dispResolIndex;
		int s_max = s_disp;
		int s_min = s_disp;
		if (resolSwitching) {
			s_max = IloMin(s_disp+1,nbResol);
			s_min = IloMax(s_disp-1,0);
		}
		int v_aux = videosIndeces[u];
		userMinRateDueUserResol[u] = RDparam_minRate[v_aux][s_disp][s_min];
		userMaxRateDueUserResol[u] = RDparam_maxRate[v_aux][s_disp][s_max];
				
		// connection assignment --> allowed rates assignment
		int allowedUserConnections[] ={0, 0, 0, 0, 0};
		float allowedConnectionProbability = 0;
		float notAllowedConnectionProbability = 0;
		for (c=0; c<nbConnections; c++){
			if (((userConnectionsSet[c].maxBandwidth-IloMax(userMinRateDueUserResol[u],userConnectionsSet[c].minBandwidth))/(userConnectionsSet[c].maxBandwidth-userConnectionsSet[c].minBandwidth)) >= 0.1) {
				allowedUserConnections[c] = 1;
				allowedConnectionProbability += userConnectionsSet[c].distribution;  
			}
			else {
				notAllowedConnectionProbability += userConnectionsSet[c].distribution;  
			}
		}
		float cdf = 0;		
		double thisUserConnecDistribution[] ={-1, -1, -1, -1, -1}; 
		for (c=0; c<nbConnections; c++){
			if (allowedUserConnections[c] == 1) {
				thisUserConnecDistribution[c] = (userConnectionsSet[c].distribution/allowedConnectionProbability)*notAllowedConnectionProbability + userConnectionsSet[c].distribution + cdf;
				cdf = thisUserConnecDistribution[c]; 
				//	cout << "user " << u << ", c " << connections[c] << ", thisUserConnecDistribution[c] " << thisUserConnecDistribution[c] << "\n";
			}
		}
		
		int randConnection = rand() % 100 + 1; 
		for (c=0; c<nbConnections; c++) {
			if (randConnection <= thisUserConnecDistribution[c]) {
				connectionsIndeces[u] = c;
				usersSet[u].connectionIndex = c;
				break;
			}
		}
		int upperBoundCap = userConnectionsSet[c].maxBandwidth;
		int lowerBoundCap = IloMax(userMinRateDueUserResol[u],userConnectionsSet[c].minBandwidth);			
		usersSet[u].capacity = (rand() % (upperBoundCap - lowerBoundCap+1)) + lowerBoundCap + 1;
		
		IloNum realUserCapacity = IloMin(userMaxRateDueUserResol[u],usersSet[u].capacity);
		int counter = 0;
		for (l=0; l<nbRepres; l++){
			if ((representationsSet[l].bitRate <= realUserCapacity) && (representationsSet[l].bitRate >= realUserCapacity)) {
				userMaxRate[u] = representationsSet[l].bitRate; 
				counter++;
			}
		}
		//cout << "user " << u << ", c " << connections[c] << ", d " << devices[d] << ", v " << usersSet[u].videoIndex << ", user capacity " << usersSet[u].capacity << ", real user capacity " << userMaxRate[u] << ", counter " << counter <<"\n\n";
		
		if (upperBoundCap < lowerBoundCap){
			cerr << "ERROR: upperBoundCap < lowerBoundCap at user:  " << u << "\n";
			cout << "user " << u << " " << lowerBoundCap << " " << upperBoundCap << " " << usersSet[u].capacity <<"\n";
			cout << userConnectionsSet[c].maxBandwidth << " " << userConnectionsSet[c].minBandwidth << " " << userMinRateDueUserResol[u]  <<"\n";
			//	cout << userConnectionsSet[c].maxBandwidth-userMinRateDueUserResol[u] << " " << userConnectionsSet[c].maxBandwidth-userConnectionsSet[c].minBandwidth << " " << (userConnectionsSet[c].maxBandwidth-userMinRateDueUserResol[u])/(userConnectionsSet[c].maxBandwidth-userConnectionsSet[c].minBandwidth) <<"\n";
			throw(-1);
		}
	}
	
	// Definition of auxiliar QoE function f_aux(s_disp, v, l) and allowedRepresentations_aux(s_disp, v, l)
	FloatMatrix satisfaction_aux_3D(env, nbResol); // s_disp x v_l
	FloatMatrix allowedRepresentations_aux_3D(env, nbResol); // s_disp x v_l
	int negativeSatisCounter = 0;
	for (int s_disp=0; s_disp<nbResol; s_disp++){
		satisfaction_aux_3D[s_disp] = IloNumArray(env, nbRepres);
		allowedRepresentations_aux_3D[s_disp] = IloNumArray(env, nbRepres);
		
		for (l=0; l<nbRepres; l++){

			int v = representationsSet[l].videoIndex;
			int s_aux = representationsSet[l].resolutionIndex;
			int r = representationsSet[l].bitRateIndex;
			float rate = representationsSet[l].bitRate;
			
			float a = RDparam_a[v][s_disp][s_aux];
			float b = RDparam_b[v][s_disp][s_aux];
			float c = RDparam_c[v][s_disp][s_aux];
			
			float QoE_fitting= 1-(a+(b/(rate+c)));
			if (QoE_fitting > 1) {
				QoE_fitting = 1;
			}
			
			if (a == -1) {
				continue;
			}
			if (QoE_fitting < 0) {
				continue;
			}
			
			satisfaction_aux_3D[s_disp][l] = QoE_fitting;
			allowedRepresentations_aux_3D[s_disp][l] = 1;
			
			if (satisfaction_aux_3D[s_disp][l] < 0){ 
                cerr << "negative sat[s_disp=" << s_disp << "][v=" << v<< "][s=" << s_aux << "]: " << satisfaction_aux_3D[s_disp][l] << "\n";
				negativeSatisCounter ++;
			}
			if (satisfaction_aux_3D[s_disp][l] > 1) 
				cout << "over 1 sat[s_disp=" << s_disp << "][v=" << v<< "][s=" << s_aux << "]: " << satisfaction_aux_3D[s_disp][l] << "\n";			
		}
	}
	cout << " negativeSatisCounter " <<  negativeSatisCounter  << "\n";
	cout << " satisfaction_aux_3D computed"  << "\n";
	
	// Input Parameters to the ILP model
	
    // satisfaction_2D:  satisfacction level experienced by user 'u' watching representation 'l'; UxL matrix
    // allowedRepresentations_2D: 1 if video representation 'l' is compatible with user 'u' features, that is, if user display can support resolution and rate of representation 'l'; UxL matrix
	FloatMatrix satisfaction_2D(env, nbUsers); // u x vl
	FloatMatrix allowedRepresentations_2D(env, nbUsers); // u x vl
	for (u=0; u<nbUsers; u++){
		satisfaction_2D[u] = IloNumArray(env, nbRepres);
		allowedRepresentations_2D[u] = IloNumArray(env, nbRepres);
		for (l=0; l<nbRepres; l++){
			int this_s_disp = usersSet[u].dispResolIndex;
			if (satisfaction_aux_3D[this_s_disp][l] > 0) {
				
				if (resolSwitching){
					satisfaction_2D[u][l] = satisfaction_aux_3D[this_s_disp][l];
					allowedRepresentations_2D[u][l] = allowedRepresentations_aux_3D[this_s_disp][l];
				}
				else {
					if (this_s_disp == representationsSet[l].resolutionIndex) {
						satisfaction_2D[u][l] = satisfaction_aux_3D[this_s_disp][l];
						allowedRepresentations_2D[u][l] = allowedRepresentations_aux_3D[this_s_disp][l];
					}
				}
			}
		}
		if (IloSum(allowedRepresentations_2D[u]) == 0){
			cerr << "user: " << u << " has not any allowed representation \n";	
			throw(-1);
		}
	}	
    //bitRate (b_l in R+) : encoding rate (kbps) of video representation 'l'; 1xL vector
    //repRates_ids (r_l in Z) : encoding rate index of video representation 'l'; 1xL vector	IloNumArray bitRate(env,nbRepres);
	IloNumArray repRates_ids(env,nbRepres);
	for (l=0; l<nbRepres; l++){	
		bitRate[l] =representationsSet[l].bitRate;
		repRates_ids[l] =representationsSet[l].bitRateIndex;
	}
    //linkCapacity (c_u parameters in R+) : capacity (kbps) of internet connection of user  'u'; 1xU vecto
	IloNumArray linkCapacity(env,nbUsers); 
	for(u = 0; u < nbUsers; u++)
		linkCapacity[u] = usersSet[u].capacity;
	
    // demand(d_uv parameters): 1 if user 'u' demands video 'v'; 0, otherwise; UxV binary matrix
	FloatMatrix demand(env, nbUsers);
	IloNumArray videoAssignment(env, nbUsers);
	for (u = 0; u < nbUsers; u++){
		demand[u] = IloNumArray(env, nbVideos);
		demand[u][usersSet[u].videoIndex] = 1;
		videoAssignment[u] = usersSet[u].videoIndex;
		if (userMinRateDueUserResol[u] > usersSet[u].capacity) {
			cout << "user " << u << "\n";
			cout << "c " << connections[usersSet[u].connectionIndex] << "\n";
			cout << "d " << devices[usersSet[u].deviceIndex] << "\n";
			cout << "usersSet[u].capacity " << usersSet[u].capacity << "\n";
			cout << "userMinRateDueUserResol[u] " << userMinRateDueUserResol[u] << "\n";

			cerr << "user: " << u << " is not feasible \n";	
			throw(-1);
		}
	}	
	cout << "nbRepres: " << nbRepres << "\n";

    // repResolutions_ids (s_l in Z): spatial resolution index of video representation 'l'; 1xL vector
    // repResolutions (p_l): spatial resolution value (e.g. 1080p) of video representation 'l'; 1xL vector
    // repVideos(v_l in Z): video index of video representation 'l'; 1xL vector
    IloNumArray repResolutions(env,nbRepres);
    IloNumArray repResolutions_ids(env,nbRepres);
    IloNumArray repVideos(env,nbRepres);
    for (l=0; l<nbRepres; l++){
        repResolutions_ids[l] = representationsSet[l].resolutionIndex;
        repResolutions[l] = representationsSet[l].resolution;
        repVideos[l] = representationsSet[l].videoIndex;
    }

	// Saving in scenario.dat the input parameters of the model.	
	string scenarioFullPathString = "scenarios_TOMCAP/" +  scenarioString + "_OPT.dat";
	const char* filename = scenarioFullPathString.c_str();
	ofstream file(filename);
	if (!file) {
		cerr << "ERROR: could not open file '" << filename
		<< "' for writing" << endl;
		throw(-1);
	}
    file << satisfaction_2D << demand << allowedRepresentations_2D << bitRate << linkCapacity << repRates_ids << repResolutions << repResolutions_ids <<repVideos << videosIndeces <<endl;

	env.end();
}
