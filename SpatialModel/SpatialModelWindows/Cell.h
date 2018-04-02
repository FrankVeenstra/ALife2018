#pragma once

#include <string>
#include "Utility.h" 

using namespace std ;

class Cell
{
public:
	Cell(float mr, float dr, float rr,  int const& tp, int swaps);
	~Cell();
//	void setHealth(int const& health) ;
//	void setMaxHealth(int const& maxHealth) ;
//	void setForce(int const& force) ;
//	void setMaxForce(int const& maxForce) ;
//	void setTeam(int const& team) ;
//	int getHealth() const ;
//	int getMaxHealth() const ;
//	int getForce() const ;
//	int getMaxForce() const ;
//	int getTeam() const ;
	void mutate(float mr);
	void mutMortality(float mr);
	void mutMaxAge(float mr);
	void mutDeathRate(float mr);
	void initalizeGene(int geneLength);
	void initalizeGene(int geneLength, int swaps) {};
	void setParams(int senseRange, vector<bool> gene, float mutationRate,
		bool mutationDeathRate, float deathRate, int movementRange, int type, 
		float reproductionRate, int foodType, int predatorType, float fit, vector<vector<bool>> genes,
		vector<int> dTime, int age_v, int cs, int dAge, bool im, int maxAge_v);

	vector<int> getRandomMovement(); 
	void updatePredMass(float amount);
	float updatePreyMass(float amount, float maxMass);
	void updateMass(float amount);
	
	float deathGradient = 0.0; // probability of a cell to die in a specific region.
	//vector<int> gene;
	float mutationRate = 0.001;
	bool mutateDeathRate = false;
	float deathRate = 0.0;
	bool immortal = true;
	int movementRange = 1; // movementSteps per turn
	int type = 0;
	float reproductionRate;
	int foodType = 0; 
	int predatorType = 0; 
	float randFloat(float low, float high);

	vector<bool> cellGene; // 16 bit // currently expressed
	void createCellGene(int length);

	float evaluateGene();
	int senseRange = 0; 
	float fitness = 0.0;
	float predMass = 0.0; 
	float preyMass = 0.0; 
	float mass; 

	// temporary
	int cStage = 0;
	// bool developed = false;
	vector<int> developTime;
	int age = 0;
	int dAge = 0;
	int maxAge = 1000;
	vector<vector<bool>> genes;

private:
	//int m_health ;
	//int m_maxHealth ;
	//int m_force ;
	//int m_maxForce ;
	//int m_team ;
};

