//#pragma comment(linker, "/SUBSYSTEM:windows /ENTRY:mainCRTStartup")

#pragma once
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <random>
#include <fstream>
#include "Cell.h"
#include "DevelopingCell.h"
#include "Utility.h" 
#include "Button.h"
#include "Slider.h"
#include <numeric>
#include <string>
#include <functional>
#include <algorithm>
#include <iterator>

using namespace std;
using namespace sf;

#define E(i,j)
#define EMPTY_CELL 0,0,0,0,0
#define EMPTYPARAMS 0,vector<bool>(),0.0f,false,0.0f,0,0,0.0f,0,0,0.0, vector<vector<bool>>(), vector<int>(),0,0,0,true,0

// settings. Will be overwritten when the settings file is read. 
bool includePredator = false;
int example = 2;
float biomassProduction = 0.0008; // 0.0002, 0.001, 0.002, 0.004, 0.008
float mutRate = 0.05;
float mutRatePred = mutRate;
string filename = "thus";
int startWithOptimal = 0; // 0 random, 1 midle, 2 optimal
bool createDeathGradient = false;
bool createMRGradient = false;
float preyMassLoss = 0.01; // how much mass the prey loses every time step
float predMassLoss = 0.01;
float repRate = 0.2; //0.2; // 0.02
bool evolveCenter = true;
bool displayVegetation = true;
float preyEatEfficiency = 0.1;
float initialBiomass = 0.4;
// mortality
bool initialMortalityState = false; // immortal = 0
bool evolveMortality = false;
float mortalMutRate = 0.001;
// max age
bool useMaxAge = false; 
bool evolveMaxAge = false;
float maxAgeMutRate = 0.02;
int initialMaxAge = 1000;
// Death Rate
bool evolveDeathRate = false;
bool useDeathRate = false; 
float deathRate = 0.001;
float deathMR = 0.01;
// Development
bool evolveDevelopmentTiming = false; 
int devSwaps = 0;
float costOfReproduction = 0.4; 
float maxPreyMass = 1.0;
float minimumPreyMass = 0.01; 
bool cannibalism = true;

float maxRewardFactor = 2.0;
bool headless = false; 
int totalInds = 0;
bool maxFound = false;
float movementSpeed = 0.2;

// settings that were in main:
int width = 200;
int height = 200;
int genesize = 32;
int maxFit = 80;
float ratio = width / (float)height;
int rowSize = 200;
int columnSize = (int)ratio * rowSize;
float dx = height / rowSize;

int iter1 = 0;
int iter2 = 0;
int stopAt = 1000;
ofstream myfile;

float colorFactor = 0.008;

int totalCycles = 0;
vector<vector<DevelopingCell>> board;

#ifdef __linux__ 
// linux starts in headless mode by default
headless = false;
#elif _WIN32
#include <windows.h>
#include <SFML/Graphics.hpp>
#else
#endif

float randFloat(float low, float high) {
	return (low + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (high - low))));
}

inline bool fileExists(const std::string& name) {
	ifstream f(name.c_str());
	return f.good();
}

void newCycle() {
	vector<vector<int>> skippers;
	for (int i = 0; i < columnSize; i++)
	{
		for (int j = 0; j < rowSize; j++)
		{
			// make sure that you don't move the same cell twice. skippers keeps track of the cells that moved. 
			bool skip = false;
			for (int n = 0; n < skippers.size(); n++) {
				if (i == skippers[n][0]) {
					if (j == skippers[n][1]) {
						skip = true;
					}
				}
			}
			if (skip == false) {
				//	cout << "type " << i << "," << j << "="<<board[i][j].type << endl;
				if (board[i][j].type == 0) // primary energy producers
				{
					board[i][j].updateMass(biomassProduction);
				}
				if (board[i][j].type == 1) // prey
				{
					// update development of gene
					if (devSwaps > 0) {
						if (board[i][j].cStage < devSwaps && board[i][j].dAge > board[i][j].developTime[board[i][j].cStage]) {
							board[i][j].cStage = board[i][j].cStage + 1;
							board[i][j].dAge = 0;
							board[i][j].developGene(board[i][j].cStage);
							board[i][j].fitness = board[i][j].evaluateGene();
							if (board[i][j].fitness >= maxFit) {
								maxFound = true;
							}
						}
					}
					board[i][j].age += 1;
					board[i][j].dAge += 1;
					vector<int> mov;
					if (randFloat(0, 1.0) < movementSpeed) {
						mov = board[i][j].getRandomMovement();
					}
					else {
						mov.push_back(0);
						mov.push_back(0);
					}
					board[i][j].updatePreyMass(-preyMassLoss, maxPreyMass);
					if (board[i][j].immortal == false) {
						if (board[i][j].deathRate > randFloat(0, 1)) {
							board[i][j].setParams(EMPTYPARAMS);
							board[i][j].updateMass(board[i][j].preyMass * 0.5);
							board[i][j].preyMass = 0.0;
							board[i][j].type = 0;
							board[i][j].fitness = 0;
							//board[i][j].mass = 0;
						}
						else if (useMaxAge == true && board[i][j].maxAge < board[i][j].age) {
							board[i][j].setParams(EMPTYPARAMS);
							board[i][j].updateMass(board[i][j].preyMass * 0.5);
							board[i][j].preyMass = 0.0;
							board[i][j].type = 0;
							board[i][j].fitness = 0;
						}
					}
					if (board[i][j].preyMass < minimumPreyMass) { // dies when not enough mass
						board[i][j].setParams(EMPTYPARAMS);
						board[i][j].updateMass(board[i][j].preyMass * 0.5);
						board[i][j].preyMass = 0.0;
						board[i][j].type = 0;
						board[i][j].fitness = 0;
					}
					// dies when more than maxAge

					else if ((i + mov[0]) < rowSize &&
						(i + mov[0]) >= 0 &&
						(j + mov[1]) < columnSize &&
						(j + mov[1]) >= 0)
					{
						vector<int> XY;
						XY.push_back(i + mov[0]);
						XY.push_back(j + mov[1]);
						skippers.push_back(XY);
						XY.clear();

						if (board[i + mov[0]][j + mov[1]].type == 0) {
							board[i + mov[0]][j + mov[1]].setParams(board[i][j].senseRange, board[i][j].cellGene,
								board[i][j].mutationRate, board[i][j].mutateDeathRate, board[i][j].deathRate,
								board[i][j].movementRange, board[i][j].type, board[i][j].reproductionRate,
								board[i][j].foodType, board[i][j].predatorType, board[i][j].fitness, board[i][j].genes,
								board[i][j].developTime, board[i][j].age, board[i][j].cStage, board[i][j].dAge, board[i][j].immortal, board[i][j].maxAge);
							board[i + mov[0]][j + mov[1]].preyMass = board[i][j].preyMass;
							float plantMass = 0;
							if (board[i + mov[0]][j + mov[1]].fitness >= maxFit) {
								plantMass = board[i + mov[0]][j + mov[1]].updatePreyMass(board[i + mov[0]][j + mov[1]].mass * (board[i + mov[0]][j + mov[1]].fitness) / maxFit * preyEatEfficiency * maxRewardFactor, maxPreyMass);
							}
							else {
								plantMass = board[i + mov[0]][j + mov[1]].updatePreyMass(board[i + mov[0]][j + mov[1]].mass * (board[i + mov[0]][j + mov[1]].fitness) / maxFit * preyEatEfficiency, maxPreyMass);
							}
							board[i + mov[0]][j + mov[1]].mass = plantMass;
							// reproduce 
							if (board[i][j].reproductionRate > randFloat(0, 1) && board[i][j].preyMass > (2*costOfReproduction)) {
								board[i + mov[0]][j + mov[1]].updatePreyMass(-costOfReproduction, maxPreyMass);
								board[i][j].preyMass = costOfReproduction * 0.8;

								board[i][j].mutate(mutRate);
								if (evolveMortality == true) {
									board[i][j].mutMortality(mortalMutRate);
								}
								if (evolveMaxAge == true) {
									board[i][j].mutMaxAge(maxAgeMutRate);
								}
								if (evolveDeathRate == true) {
									if (board[i][j].immortal == false) {
										board[i][j].mutDeathRate(deathMR);
									}
								}
								if (evolveDevelopmentTiming == true) {
									board[i][j].mutDevelopmentTime(mutRate);
								}
								board[i][j].cStage = 0;
								board[i][j].age = 0;
								board[i][j].dAge = 0;
								board[i][j].fitness = board[i][j].evaluateGene();
								if (board[i][j].fitness >= maxFit) {
									maxFound = true;
								}
								totalInds++; 
							}
							else {
								board[i][j].setParams(EMPTYPARAMS); // no reproduction
								board[i][j].preyMass = 0.0;
							}
						}
						else if (cannibalism == true && board[i + mov[0]][j + mov[1]].type == 1 && board[i + mov[0]][j + mov[1]].fitness <= board[i][j].fitness) {
							if (mov[0] != 0 || mov[1] != 0) {
								board[i][j].updatePreyMass(board[i + mov[0]][j + mov[1]].preyMass * (board[i + mov[0]][j + mov[1]].fitness) / maxFit * preyEatEfficiency, maxPreyMass);
								board[i + mov[0]][j + mov[1]].setParams(board[i][j].senseRange, board[i][j].cellGene,
									board[i][j].mutationRate, board[i][j].mutateDeathRate, board[i][j].deathRate,
									board[i][j].movementRange, board[i][j].type, board[i][j].reproductionRate,
									board[i][j].foodType, board[i][j].predatorType, board[i][j].fitness, board[i][j].genes,
									board[i][j].developTime, board[i][j].age, board[i][j].cStage, board[i][j].dAge, board[i][j].immortal, board[i][j].maxAge);
								board[i + mov[0]][j + mov[1]].preyMass = board[i][j].preyMass;
								//cout << "CANNIBALISM EVENT: " << cannibalism << endl;
								// reproduce 
								if (board[i][j].reproductionRate > randFloat(0, 1) && board[i][j].preyMass > (2 * costOfReproduction)) {
									board[i + mov[0]][j + mov[1]].updatePreyMass(-costOfReproduction, maxPreyMass);
									board[i][j].preyMass = costOfReproduction * 0.8;

									board[i][j].mutate(mutRate);
									if (evolveMortality == true) {
										board[i][j].mutMortality(mortalMutRate);
									}
									if (evolveMaxAge == true) {
										board[i][j].mutMaxAge(maxAgeMutRate);
									}
									if (evolveDeathRate == true) {
										if (board[i][j].immortal == false) {
											board[i][j].mutDeathRate(deathMR);
										}
									}
									if (evolveDevelopmentTiming == true) {
										board[i][j].mutDevelopmentTime(mutRate);
									}
									board[i][j].cStage = 0;
									board[i][j].age = 0;
									board[i][j].dAge = 0;
									board[i][j].fitness = board[i][j].evaluateGene();
									if (board[i][j].fitness >= maxFit) {
										maxFound = true;
									}
									totalInds++;
								}
								else {
									board[i][j].setParams(EMPTYPARAMS); // no reproduction
									board[i][j].preyMass = 0.0;
								}
							}
						}
					}
				}
				else if (board[i][j].type == 2)
				{
					vector<int> mov = board[i][j].getRandomMovement();
					board[i][j].updatePredMass(-0.001);
					if (board[i][j].deathRate > randFloat(0, 1)) {
						board[i][j].setParams(EMPTYPARAMS);
						board[i][j].updateMass(board[i][j].predMass * 0.2);
						board[i][j].predMass = 0.0;
						board[i][j].type = 0;
						board[i][j].fitness = 0;
					}
					else if (board[i][j].predMass < 0.1) {
						board[i][j].setParams(EMPTYPARAMS);
						board[i][j].updateMass(board[i][j].predMass * 0.2),
							board[i][j].predMass = 0.0;
						board[i][j].type = 0;
						board[i][j].fitness = 0;

					}
					//					r1 = distrib(gen)-1 ;
					//					r2 = distrib(gen)-1 ;

					else if ((i + mov[0]) < rowSize &&
						(i + mov[0]) >= 0 &&
						(j + mov[1]) < columnSize &&
						(j + mov[1]) >= 0)
					{
						if (board[i + mov[0]][j + mov[1]].type == 0) {
							board[i + mov[0]][j + mov[1]].setParams(board[i][j].senseRange, board[i][j].cellGene,
								board[i][j].mutationRate, board[i][j].mutateDeathRate, board[i][j].deathRate,
								board[i][j].movementRange, board[i][j].type, board[i][j].reproductionRate,
								board[i][j].foodType, board[i][j].predatorType, board[i][j].fitness,
								board[i][j].genes, board[i][j].developTime, board[i][j].age, board[i][j].cStage,
								board[i][j].dAge, board[i][j].immortal, board[i][j].maxAge);
							board[i + mov[0]][j + mov[1]].predMass = board[i][j].predMass;
							board[i][j].setParams(EMPTYPARAMS);
							board[i][j].predMass = 0.0;
						}
						else if (board[i + mov[0]][j + mov[1]].type == 1) {
							board[i + mov[0]][j + mov[1]].setParams(board[i][j].senseRange, board[i][j].cellGene,
								board[i][j].mutationRate, board[i][j].mutateDeathRate, board[i][j].deathRate,
								board[i][j].movementRange, board[i][j].type, board[i][j].reproductionRate,
								board[i][j].foodType, board[i][j].predatorType, board[i][j].fitness,
								board[i][j].genes, board[i][j].developTime, board[i][j].age, board[i][j].cStage,
								board[i][j].dAge, board[i][j].immortal, board[i][j].maxAge);
							board[i + mov[0]][j + mov[1]].predMass = board[i][j].predMass;

							board[i + mov[0]][j + mov[1]].updatePredMass(((board[i + mov[0]][j + mov[1]].preyMass) * 0.1));
							board[i + mov[0]][j + mov[1]].preyMass = 0.0;
							board[i + mov[0]][j + mov[1]].mass = board[i + mov[0]][j + mov[1]].preyMass * 0.3;

							if (board[i][j].reproductionRate > randFloat(0, 1) && board[i][j].predMass > 0.6) {
								board[i + mov[0]][j + mov[1]].updatePredMass(-0.2);
								board[i][j].predMass = 0.2;
							}
							else {
								board[i][j].setParams(EMPTYPARAMS); // no reproduction
								board[i][j].predMass = 0;
							}
						}
					}
				}
			}
		}
	}
	skippers.clear();
}

int getGenomeDif(vector<bool> gen) {
	int divValue = 0;
	for (int j = 0; j < gen.size(); j++) {
		if (gen[j] == 1) {
			divValue++;
		}
	}
	return divValue;
}

bool log() {
	totalCycles += 1;
	vector<int> totalDevTime;
	for (int m = 0; m < devSwaps;m++) {
		totalDevTime.push_back(0);
	}
	vector<int> avFit;
	for (int m = 0; m < devSwaps;m++) {
		avFit.push_back(0);
	}
	vector<bool> maxFitAchieved;
	for (int m = 0; m < devSwaps + 1;m++) {
		maxFitAchieved.push_back(false);
	}
	int totalAge = 0;
	int maxAge = 0;
	int amountPlants = 0;
	int amountImmortals = 0;
	int amountMortals = 0;
	int amountPrey = 0;
	int amountPredator = 0;
	float amountPlantBiomass = 0;
	float amountPreyBiomass = 0;
	float amountPredBiomass = 0;
	float fittestind = 0; 
	// DIVERSITY MEASURE (NEW)
	vector<int> divValues;
	// -----------------------

	bool displayMaxGenome = false; 
	vector<float> allFitnesses;
	vector<int> allAges;
	vector<float> allDeathRates;
	vector<float>allmaxAge;
	//	vector<float> allDeathRates;
	vector<vector<int>> allDevelopmentalTimes;
	for (int i = 0; i < rowSize; i++) {
		for (int j = 0; j < rowSize; j++) {
			if (board[i][j].type == 0) {
				amountPlants += 1;
			}
			if (board[i][j].type == 1) {
				divValues.push_back(getGenomeDif(board[i][j].cellGene));

				if (board[i][j].fitness > fittestind) {
					fittestind = board[i][j].fitness;
				}
				allFitnesses.push_back(board[i][j].fitness);
				allAges.push_back(board[i][j].age);
				allmaxAge.push_back(board[i][j].maxAge);
				if (board[i][j].age > maxAge) {
					maxAge = board[i][j].age;
				}
				amountPrey += 1;
				if (board[i][j].immortal == true) {
					amountImmortals += 1;
				}
				else {
					amountMortals += 1;
					allDeathRates.push_back(board[i][j].deathRate);
				}
				if (board[i][j].cStage == 0) {
					if (board[i][j].fitness >= maxFit) {
						 if (displayMaxGenome == false) {
							for (int x = 0; x < board[i][j].cellGene.size(); x++) {
								cout << "[" << board[i][j].cellGene[x] << "]";
							}
							cout << endl; 
							displayMaxGenome = true;
						}
						maxFitAchieved[0] = true;
					}
				}
				for (int m = 0; m < board[i][j].developTime.size();m++) {
					totalDevTime[m] += board[i][j].developTime[m];
					if (board[i][j].cStage == m + 1) {
						if (board[i][j].fitness >= maxFit) {
							if (displayMaxGenome == false) {
								cout << "D:" << board[i][j].cStage << ","; 
								for (int x = 0; x < board[i][j].cellGene.size(); x++) {
									cout << "[" << board[i][j].cellGene[x] << "]";
								}
								cout << endl;
								displayMaxGenome = true;
							}
							maxFitAchieved[m + 1] = true;
						}
					}
				}
			}
			if (board[i][j].type == 2) {
				amountPredator += 1;
			}
			amountPlantBiomass += board[i][j].mass;
			amountPreyBiomass += board[i][j].preyMass;
			amountPredBiomass += board[i][j].predMass;
		}
	}
	stringstream ss;
	ss << "cycle:," << totalCycles << ",";
	ss << "Amount Mortals," << amountMortals << "," << "Amount Immortals," << amountImmortals << ",";
	ss << "PopSize," << amountPrey << ", ";

	float totBiomass = amountPlantBiomass + amountPreyBiomass + amountPredBiomass;
	if (example == 1) {
		ss << "Amount Creatures," << amountPrey << ", ";
		//	cout << "Biomass: Plant = " << amountPlantBiomass << " , creature = " << amountPreyBiomass << endl;
	}
	else if (example == 2) {
		//	cout << "Amount Prey " << amountPrey << " and amount Predators = " << amountPredator << ", ";
		//	cout << "Biomass: Plant = " << amountPlantBiomass << " , Prey = " << amountPreyBiomass << " , Pred = " << amountPredBiomass << " , tot = " << totBiomass << endl;
	}
	if (amountPrey > 0) {
		vector<int> averageDevelopTime;
		for (int m = 0; m < devSwaps;m++) {
			averageDevelopTime.push_back(totalDevTime[m] / amountPrey);
		}
		if (maxFitAchieved[0] == true) {
			ss << "mf,1,";
		}
		else {
			ss << "mf,0,";
		}
		ss << "AvDevT,";
		for (int m = 0; m < devSwaps;m++) {
			ss << "dt," << averageDevelopTime[m] << ",";

			if (maxFitAchieved[m + 1] == true) {
				ss << "mf,1,";
			}
			else {
				ss << "mf,0,";
			}
		}
		maxAge = 0; 
		for (size_t i = 0; i < allAges.size(); i++)
		{
			if (allAges[i] > maxAge) {
				maxAge = allAges[i];
			}
		}
		int agesum = std::accumulate(allAges.begin(), allAges.end(), 0.0);
		double agemean = agesum / allAges.size();
		double agesq_sum = std::inner_product(allAges.begin(), allAges.end(), allAges.begin(), 0.0);
		double agestdev = std::sqrt(agesq_sum / allAges.size() - agemean * agemean);
		double fitsum = std::accumulate(allFitnesses.begin(), allFitnesses.end(), 0.0);
		double fitmean = fitsum / allFitnesses.size();
		double fitsq_sum = std::inner_product(allFitnesses.begin(), allFitnesses.end(), allFitnesses.begin(), 0.0);
		double fitstdev = std::sqrt(fitsq_sum / allFitnesses.size() - fitmean * fitmean);

		vector<int> sortedAges = allAges;
		sort(sortedAges.begin(), sortedAges.end());
		int sortASize = sortedAges.size();
		double fp0A = sortedAges[0];
		double fp25A = sortedAges[int(sortASize / 4)];
		double fp50A = sortedAges[int(sortASize / 2)];
		double fp75A = sortedAges[int(sortASize / 4 * 3)];;
		double  fp100A = sortedAges[sortASize-1];

		vector<float> sortedFitnesses = allFitnesses;
		sort(sortedFitnesses.begin(), sortedFitnesses.end());
		int sortSize = sortedFitnesses.size();
		double fp0 = sortedFitnesses[0];
		double fp25 = sortedFitnesses[int(sortSize / 4)];
		double fp50 = sortedFitnesses[int(sortSize / 2)];
		double fp75 = sortedFitnesses[int(sortSize / 4 * 3)];;
		double fp100 = sortedFitnesses[sortSize-1];

		//cout << "death log: " << allDeathRates.size() << endl;
		if (allDeathRates.size() == 0) {
			allDeathRates.push_back(0);
		}
		double deathsum = std::accumulate(allDeathRates.begin(), allDeathRates.end(), 0.0);
		double deathmean = deathsum / allDeathRates.size();
		double deathsq_sum = std::inner_product(allDeathRates.begin(), allDeathRates.end(), allDeathRates.begin(), 0.0);
		double deathstdev = std::sqrt(deathsq_sum / allDeathRates.size() - deathmean * deathmean);
	
		vector<float> sortedDeath = allDeathRates;
		sort(sortedDeath.begin(), sortedDeath.end());
		int sortDSize = sortedDeath.size();
		double fp0D = sortedDeath[0];
		double fp25D = sortedDeath[int(sortDSize / 4)];
		double fp50D = sortedDeath[int(sortDSize / 2)];
		double fp75D = sortedDeath[int(sortDSize / 4 * 3)];;
		double  fp100D = sortedDeath[sortDSize - 1];

		double maxAgesum = std::accumulate(allmaxAge.begin(), allmaxAge.end(), 0.0);
		double maxAgemean = maxAgesum / allmaxAge.size();
		double maxAgesq_sum = std::inner_product(allmaxAge.begin(), allmaxAge.end(), allmaxAge.begin(), 0.0);
		double maxAgestdev = std::sqrt(maxAgesq_sum / allmaxAge.size() - maxAgemean * maxAgemean);

		vector<float> sortedmaxAge = allmaxAge;
		sort(sortedmaxAge.begin(), sortedmaxAge.end());
		int sortmaxAgeSize = sortedmaxAge.size();
		double fp0maxAge = sortedmaxAge[0];
		double fp25maxAge = sortedmaxAge[int(sortmaxAgeSize / 4)];
		double fp50maxAge = sortedmaxAge[int(sortmaxAgeSize / 2)];
		double fp75maxAge = sortedmaxAge[int(sortmaxAgeSize / 4 * 3)];;
		double  fp100maxAge = sortedmaxAge[sortmaxAgeSize - 1];


		// DIVERSITY MEASURE (NEW)
		double divsum = std::accumulate(divValues.begin(), divValues.end(), 0.0);
		double divmean = divsum / divValues.size();
		double divsq_sum = std::inner_product(divValues.begin(), divValues.end(), divValues.begin(), 0.0);
		double divstdev = std::sqrt(divsq_sum / divValues.size() - divmean * divmean);

		vector<int> sortedDiv = divValues;
		sort(sortedDiv.begin(), sortedDiv.end());
		int sortSizeDiv = sortedDiv.size();
		double fp0Div = sortedDiv[0];
		double fp25Div = sortedDiv[int(sortSize / 4)];
		double fp50Div = sortedDiv[int(sortSize / 2)];
		double fp75Div = sortedDiv[int(sortSize / 4 * 3)];;
		double  fp100Div = sortedDiv[sortSize - 1];
		// END DIVERSITY MEASURE

		//ss << "AvAge," << averageAge << "," << "mAge," << maxAge << ",";
		ss << "AvAge," << agemean << ",std," << agestdev << ",pct," << fp0A << "," << fp25A << "," << fp50A << "," << fp75A << "," << fp100A << ",";
		ss << "AvFit," << fitmean << ",std," << fitstdev << ",pct," << fp0 << "," << fp25 << "," << fp50 << "," << fp75 << "," << fp100 << ",";
		ss << "AvDR," << deathmean << ",std," << deathstdev << ",pct," << fp0D << "," << fp25D << "," << fp50D << "," << fp75D << "," << fp100D << ",";
		ss << "AvMaxAge," << maxAgemean << ",std," << maxAgestdev << ",pct," << fp0maxAge << "," << fp25maxAge << "," << fp50maxAge << "," << fp75maxAge << "," << fp100maxAge << ",";
		ss << "t_n," << totalInds << ",maxFound," << maxFound;
		ss << ",fittestInd," << fittestind << ",";
		ss << ",diversity," << divmean << ",std," << divstdev << ",pct," << fp0Div << "," << fp25Div << "," << fp50Div << "," << fp75Div << "," << fp100Div << ",";
		ss << endl; 

		cout << ss.str();

		myfile.open(filename + ".csv", std::ios::out | std::ios::app);
		myfile << ss.str();
		myfile.close();
		vector<int>().swap(allAges);
		vector<float>().swap(allFitnesses);
	}
	else {
		cout << "Population Extinct" << endl;
		stringstream ss;
		ss << totalCycles << endl; 
		cout << ss.str();
		//myfile.open(filename + ".csv", std::ios::out | std::ios::app);
		//myfile << ss.str();
		//myfile.close();
		return true; // causes program to stop 
	}
	return false; // continue program 
}



void initializeBoard() {

	myfile.open(filename + ".csv", std::ofstream::out | std::ofstream::trunc);
	myfile.close();

	srand(time(NULL));
	vector<vector<int>> cellPositions;
	vector<int> team;
	team.push_back((int)(rowSize / 2));
	team.push_back((int)(rowSize / 2));
	cellPositions.push_back(team);
	team[0] = int(((rowSize / 2)*0.5));
	team[1] = int(((rowSize / 2)*0.5));
	cellPositions.push_back(team);
	team[0] = int(((rowSize / 2)*1.5));
	team[1] = int(((rowSize / 2)*1.5));
	cellPositions.push_back(team);
	team[0] = int(((rowSize / 2)*1.5));
	team[1] = int(((rowSize / 2)*0.5));
	cellPositions.push_back(team);
	team[0] = int(((rowSize / 2)*0.5));
	team[1] = int(((rowSize / 2)*1.5));
	cellPositions.push_back(team);
	team[0] = int(((rowSize / 2)));
	team[1] = int(((rowSize / 2)*0.5));
	cellPositions.push_back(team);
	team[0] = int(((rowSize / 2)));
	team[1] = int(((rowSize / 2)*1.5));
	cellPositions.push_back(team);
	team[0] = int(((rowSize / 2)*0.5));
	team[1] = int(((rowSize / 2)));
	cellPositions.push_back(team);
	team[0] = int(((rowSize / 2)*1.5));
	team[1] = int(((rowSize / 2)));
	cellPositions.push_back(team);
	team[0] = int(((rowSize / 2)*0.25));
	team[1] = int(((rowSize / 2)*0.25));
	cellPositions.push_back(team);

	if (includePredator == true) {
		vector<int> team2;
		team2.push_back(20);
		team2.push_back(20);
		cellPositions.push_back(team2);
	}

	for (int i = 0; i < columnSize; i++)
	{
		vector<DevelopingCell> cellRow;
		for (int j = 0; j < rowSize; j++)
		{
			bool none = true;
			for (int k = 0; k < cellPositions.size(); k++) {
				if (i == cellPositions[k][0] && j == cellPositions[k][1]) {
					//if (k != 1) {
						cellRow.push_back(DevelopingCell(mutRate, deathRate, repRate, 1, devSwaps));
						none = false;
						cellRow[cellRow.size() - 1].preyMass = 1.0;
						cellRow[cellRow.size() - 1].movementRange = 1;
						cellRow[cellRow.size() - 1].mutationRate = mutRate;
						cellRow[cellRow.size() - 1].maxAge = initialMaxAge; 
						if (initialMortalityState == false) {
							cellRow[cellRow.size() - 1].immortal = true;
						}
						else {
							cellRow[cellRow.size() - 1].immortal = false;
						}
						if (evolveCenter == true) {
							cellRow[cellRow.size() - 1].createCellGene(genesize, devSwaps); // predefined
						}
						else {
							cellRow[cellRow.size() - 1].initializeGene(genesize, devSwaps);
						}
						//					cout << cellRow[cellRow.size() - 1].cellGene.size() << endl;
						cellRow[cellRow.size() - 1].fitness = cellRow[cellRow.size() - 1].evaluateGene();
						//					cout << "fit: " << cellRow[cellRow.size() - 1].fitness;
					/*}
					else {
						cellRow.push_back(DevelopingCell(mutRatePred, 0.0, 0.2, k + 1, devSwaps));
						cellRow[cellRow.size() - 1].predMass = 2.0;
						cellRow[cellRow.size() - 1].movementRange = 2;
						none = false;
					}*/
				}
			}
			if (none == true) {
				cellRow.push_back(DevelopingCell(EMPTY_CELL));
				cellRow[cellRow.size() - 1].mass = initialBiomass;
			}
		}
		if (createDeathGradient = true) {
			for (int n = 0; n < cellRow.size(); n++) {
				cellRow[n].deathGradient = (float)i / (float)columnSize * 0.25;
				//cellRow[n].mutationRate = cellRow[n].deathGradient;
				//		cout << (float) i / (float) columnSize << ',' << columnSize << ',' << i <<  endl;
			}
		}
		board.push_back(cellRow);
		cellRow.clear();
	}
	totalInds++;
}

void run() {
	bool stop = false;
	initializeBoard();
	while (stop == false) {
		newCycle();
		if (iter1 == 100)
		{
			stop = log();
			iter1 = 0;
			iter2++;
		}
		iter1++;
		if (iter2 == stopAt) {
			stop = true;
		}
	}
}

Image visualize(Image image) {
	image.create(width, height);
	for (int i(0); i < columnSize; i++)
	{
		for (int j(0); j < rowSize; j++)
		{
			if (board[i][j].type == 0) {
				int cf = (int)(board[i][j].mass * 254);
				//int cf = (int)(board[i][j].deathGradient * 100);
				for (int k(0); k < dx; k++)
				{
					for (int n(0); n < dx; n++)
					{
						if (displayVegetation == true) {
							image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(0, cf, (int)(cf / 2)));
						}
						else {
							image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(0, 0, 0));
						}
					}
				}
			}
			if (board[i][j].type == 1) {
				int cf = (int)(board[i][j].preyMass * 254);
				if ((int)(board[i][j]).fitness >= maxFit) {
					vector<bool> gene = board[i][j].cellGene;
					if (board[i][j].cellGene[0] == 1)
					{
						for (int k(0); k < dx; k++)
						{
							for (int n(0); n < dx; n++)
							{
								if (board[i][j].cStage == 0) {
									image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(255, 255, 255));
								}
								else if (board[i][j].cStage == 1) {
									image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(255, 0, 0));
								}
								else if (board[i][j].cStage == 2) {
									image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(0, 255, 0));
								}
								else if (board[i][j].cStage == 3) {
									image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(0, 0, 255));
								}
								else if (board[i][j].cStage == 4) {
									image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(255, 0, 255));
								}
								else if (board[i][j].cStage == 5) {
									image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(255, 255, 0));
								}
							}
						}
					}
					else {
						for (int k(0); k < dx; k++)
						{
							for (int n(0); n < dx; n++)
							{
								if (board[i][j].cStage == 0) {
									image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(255, 255, 255));
								}
								else if (board[i][j].cStage == 1) {
									image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(255, 0, 0));
								}
								else if (board[i][j].cStage == 2) {
									image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(0, 255, 0));
								}
								else if (board[i][j].cStage == 3) {
									image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(0, 0, 255));
								}
								else if (board[i][j].cStage == 4) {
									image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(255, 0, 255));
								}
								else if (board[i][j].cStage == 5) {
									image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(255, 255, 0));
								}
							}
						}
					}
				}
				else if ((int)(board[i][j]).immortal == false) {
					for (int k(0); k < dx; k++)
					{
						for (int n(0); n < dx; n++)
						{
							image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(0, 0, 255));
						}
					}
				}

				else {
					int fit = (int)(board[i][j].fitness / maxFit * 170);
					for (int k(0); k < dx; k++)
					{
						for (int n(0); n < dx; n++)
						{
							if (board[i][j].immortal == false) {
								image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(fit, fit, fit));
							}
							else {
								image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(fit, fit, fit));
							}
							//								if (board[i][j].cStage < devSwaps) {
							//									image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(100, 100, 100));
							//								}
							//								else if (board[i][j].cStage == devSwaps) {
							//									image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(100, 100, 100));
							//								}
						}
					}
				}
			}
			else if (board[i][j].type == 2 /*&& board[i][j].deathRate == 0.0*/) {
				int cf = (int)(board[i][j].predMass * 127);
				for (int k(0); k < dx; k++)
				{
					for (int n(0); n < dx; n++)
					{
						if (example != 1) {
							image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(255, 0, 0));
						}
						else {
							int cf = (int)(board[i][j].mass * 254);
							image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(cf, cf, cf));
						}
					}
				}
			}
			//		else if (board[i][j].type == 2 && board[i][j].deathRate > 0.005) {
			//			int colorFactor = (int)(board[i][j].deathRate * 255 * 1000);
			//			if (colorFactor > 255) {
			//				colorFactor = 255;
			//			}
			//			for (int k(0); k < dx; k++)
			//			{
			//				for (int n(0); n < dx; n++)
			//				{
			//					image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(255, colorFactor, colorFactor));
			//				}
			//			}
			//		}
			else if (board[i][j].type == 2) {
				for (int k(0); k < dx; k++)
				{
					for (int n(0); n < dx; n++)
					{
						image.setPixel((i * dx) + k, (j * dx) + n, Color::Color(255, 0, 0));
					}
				}
			}

			//for(int k(0) ; k < dx ; k++)
			//{
			//	for(int n(0) ; n < dx ; n++)
			//	{
			//		if(board[i][j].type==1)
			//		{
			//			image.setPixel(i*dx+k,j*dx+n,Color::Color(255,255,0)) ;
			//			//image.setPixel(i*dx+k,j*dx+n,Color::Color((board[E(i,j)].getForce()/colorFactor/(float)maxForce),0,0)) ;
			//		}
			//		else if(board[i][j].type==2)
			//		{
			//		}
			//	}
			//}
		}
	}
	return image;
}

void runWindow() {
	RenderWindow window(VideoMode(width, height, maxFit), "Cellular automaton");
	Image image; image.create(width, height);
	Texture texture; texture.loadFromImage(image);
	Sprite sprite; sprite.setTexture(texture, true);
	Vector2f mousePos;
	Event event;
	initializeBoard();

	while (window.isOpen())
	{
		while (window.pollEvent(event))
		{
			if (event.type == Event::Closed || event.type == Event::KeyPressed)
				window.close();
		}

		newCycle();
		if (iter1 == 100)
		{
			if (log() == true) {
				window.close();
			}
			iter1 = 0;
			iter2++;

		}
		//if (iter2 == stopAt || totalInds >= 65536) {
		//	window.close();
		//}
		iter1++;
		image = visualize(image);
		texture.update(image);
		sprite.setTexture(texture, true);
		window.clear();
		window.draw(sprite);
		window.display();
	}
}

int main(int argc, char* argv[])
{
	//	vector<int> organismType(0);
	int process = 0;
	if (process == 1) { // deprecated
		stringstream ss;
		ss << "Biomass production: " << argv[1] << ", Mutation Rate: " << argv[2] << ",Dev Swaps: " << argv[3] << ",Evolve Center:" << argv[4] << ",maxRewardFactor:" << argv[5] << "saveFileName:" << argv[6] << endl;
		string s = ss.str();
		//	SetConsoleTitle(TEXT(s));
		//	SetConsoleTitle(TEXT("Biomass: ") + argv[1] + ", Mutation Rate: " + argv[2]));
		biomassProduction = atof(argv[1]);
		mutRate = atof(argv[2]);
		devSwaps = atoi(argv[3]);
		evolveCenter = atoi(argv[4]);
		maxRewardFactor = atof(argv[5]);
		filename = argv[6];
	}
	string settingsname = "settings.csv";
	if (argc > 1) {
		settingsname = argv[1];
	}
	if (fileExists(settingsname)) {
		ifstream file(settingsname); // declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/
		vector<string> values;
		string value;
		while (file.good())
		{
			getline(file, value, ','); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
			values.push_back(string(value)); // display value removing the first and the last character from it
			cout << value;
		}
		cout << endl;
		for (int i = 0; i < values.size(); i++) {
			
			if (values[i] == "Biomass Production") {
				i++;
				biomassProduction = atof(values[i].c_str());
			}
			else if (values[i] == "Initial Biomass") {
				i++;
				initialBiomass = atof(values[i].c_str());
			}
			else if (values[i] == "Display Vegetation") {
				i++;
				displayVegetation = atoi(values[i].c_str());
			}
			else if (values[i] == "Mutation Rate") {
				i++;
				mutRate = atof(values[i].c_str());
			}
			else if (values[i] == "Developmental Swaps") {
				i++;
				devSwaps = atoi(values[i].c_str());
			}
			else if (values[i] == "Evolve Center") {
				i++;
				evolveCenter = atoi(values[i].c_str());
			}
			else if (values[i] == "Mass Loss") {
				i++;
				preyMassLoss = atof(values[i].c_str());
			}
			else if (values[i] == "Prey Eat Efficiency") {
				i++;
				preyEatEfficiency = atof(values[i].c_str());
			}
			else if (values[i] == "Max Reward Factor") {
				i++;
				maxRewardFactor = atof(values[i].c_str());
			}
			else if (values[i] == "Reproduction Rate") {
				i++;
				repRate = atof(values[i].c_str());
			}
			else if (values[i] == "Cost Of Reproduction") {
				i++;
				costOfReproduction = atof(values[i].c_str());
			}
			else if (values[i] == "Max Prey Mass") {
				i++;
				maxPreyMass = atof(values[i].c_str());
			}
			else if (values[i] == "Min Prey Mass") {
				i++;
				minimumPreyMass = atof(values[i].c_str());
			}
			else if (values[i] == "File Name") {
				i++;
				if (argc > 1) {
					filename = values[i] + argv[2];
				}
				else {
					filename = values[i];
				}
				cout << "filename = " << filename << endl; 
			}
			else if (values[i] == "Stop After") {
				i++;
				stopAt = atoi(values[i].c_str());
				cout << "stopAt = " << stopAt << endl;
			}
			else if (values[i] == "Headless") {
				i++;
				headless = atoi(values[i].c_str());
			}
			else if (values[i] == "Gene Size") {
				i++;
				genesize = atoi(values[i].c_str());
			}
			else if (values[i] == "Initial Mortality State") {
				i++;
				initialMortalityState = atoi(values[i].c_str());
			}
			else if (values[i] == "Evolve Mortality") {
				i++;
				evolveMortality = atoi(values[i].c_str());
			}
			else if (values[i] == "Mortality Mutation Rate") {
				i++;
				mortalMutRate = atof(values[i].c_str());
			}
			else if (values[i] == "Use Max Age") {
				cout << "A";
				i++;
				useMaxAge = atoi(values[i].c_str());
			}
			else if (values[i] == "Evolve Max Age") {
				cout << "B";
				i++;
				evolveMaxAge = atoi(values[i].c_str());
			}
			else if (values[i] == "Initial Max Age") {
				cout << "C";
				i++;
				initialMaxAge = atoi(values[i].c_str());
			}
			else if (values[i] == "Max Age Mutation Rate") {
				cout << "D";
				i++;
				maxAgeMutRate = atof(values[i].c_str());
			}
			else if (values[i] == "Use Death Rate") {
				i++;
				useDeathRate = atoi(values[i].c_str());
				//cout << "evolveDeathRate " << evolveDeathRate << endl;
			}
			else if (values[i] == "Initial Death Rate") {
				i++;
				deathRate = atof(values[i].c_str());
			}
			else if (values[i] == "Death Mutation Rate") {
				i++;
				deathMR = atof(values[i].c_str());
				//cout << "deathMR " << deathMR << endl;
			}
			else if (values[i] == "Evolve Death Rate") {
				i++;
				evolveDeathRate = atoi(values[i].c_str());
				//cout << "evolveDeathRate " << evolveDeathRate << endl;
			}
			else if (values[i] == "Movement Speed") {
				i++;
				movementSpeed = atof(values[i].c_str());
			}
			else if (values[i] == "Window Size") {
				i++;
				width = atoi(values[i].c_str());
				height = width;
				ratio = (float)width / (float)height;
				rowSize = width;
				columnSize = (int)ratio * rowSize;
				dx = height / rowSize;
				cout << "width,height: " << width << "," << height << "," << rowSize << "," << columnSize << endl;

			}
			else if (values[i] == "Evolve Developmental Timing") {
				i++;
				evolveDevelopmentTiming = atoi(values[i].c_str());
			}
			else if (values[i] == "Cannibalism") {
				cout << "herpsk" << endl; 
				i++;
				cannibalism = atoi(values[i].c_str());
			}
		}
	}
	if (genesize == 16) {
		maxFit = 32;
	}
	else if (genesize == 8) {
		maxFit = 12;
	}
	if (headless == true) {
		run();
	}
	else {
		runWindow();
	}

}
