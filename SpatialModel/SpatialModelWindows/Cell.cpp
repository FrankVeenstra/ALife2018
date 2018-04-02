#include "Cell.h"
#include <math.h> 
#include <iostream>

using namespace std;

//Cell::Cell(int const& maxHealth, int const& maxForce, int const& team)
//{
//	m_health = maxHealth ;
//	m_force = maxForce ;
//	m_maxHealth = maxHealth ;
//	m_maxForce = maxForce ;
//	m_team = team ;
//}

Cell::Cell(float mr, float dr, float rr, int const& tp, int devSwaps) {
//	mutationRate = mr;
	deathRate = dr;
	reproductionRate = rr; 
	type = tp;
}

Cell::~Cell()
{

}


float Cell::randFloat(float low, float high) {
	return (low + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (high - low))));
}

void Cell::mutate(float mr)
{
	for (int i = 0; i < cellGene.size(); i++) {
		if (randFloat(0.0, 1.0) < mr) {
			cellGene[i] = rand() % 2;
		}
	}

}

void Cell::mutMortality(float mr)
{
	if (randFloat(0, 1) < mr) {
		immortal = rand() % 2;
	}
}

void Cell::mutDeathRate(float mr)
{
	// mutating deathRate..
	if (randFloat(0, 1) < mr) {
		deathRate += randFloat(-deathRate,deathRate);
		if (deathRate > 1.0) {
			deathRate = 1.0;
		}
		else if (deathRate < 0.000000001) {
			deathRate = 0.000000001;
		}
	}
	if (randFloat(0, 2) < mr) {
		deathRate += randFloat(0,0.01);
	}
}

void Cell::mutMaxAge(float mr) {
	if (randFloat(0, 1) < mr) {
		maxAge += rand() % (maxAge + 1) - ((int)0.5*maxAge);
	}
	else if (randFloat(0, 2) < mr) {
		maxAge += rand() % ((maxAge * 2) + 1) - maxAge;
	}
}

void Cell::initalizeGene(int geneLength)
{
	cellGene.clear();
	for (int i = 0; i < geneLength; i++) {
		cellGene.push_back(rand() % 1);
	}
}

template <typename T>
void delete_pointed_to(T* const ptr)
{
	delete ptr;
}

void Cell::setParams(int sr, vector<bool> gn, float mr, bool mDr, float dr, 
	int moveR, int tp, float rr, int ft, int pt, float fit, vector< vector<bool> > gs, 
	vector<int> dTime, int age_v, int currentStage, int da, bool im, int maxAge_v)
{
	vector<vector<bool>>().swap(genes);
	vector<bool>().swap(cellGene);
	vector<int>().swap(developTime);
	senseRange = sr;
	cellGene = gn;
	mutationRate = mr; 
	mutateDeathRate = mDr;
	deathRate = dr; 
	movementRange = moveR; 
	type = tp; 
	reproductionRate = rr;
	foodType = ft;
	predatorType = pt; 
	fitness = fit;
	cStage = currentStage;
	developTime = dTime;
	age = age_v;
	maxAge = maxAge_v;
	genes = gs;
	dAge = da; 
	immortal = im;
	vector<vector<bool>>().swap(gs);
	vector<bool>().swap(gn);
	vector<int>().swap(dTime);
	gs.clear();
	gn.clear();
	dTime.clear();
}



vector<int> Cell::getRandomMovement()
{
	int x = 0;
	int y = 0; 
	for (int i = 0; i < movementRange; i++) {
		int movementType = rand() % 5;
		switch (movementType) {
		case 0:
			x -= 1;
			break;
		case 1:
			x += 1;
			break;
		case 2:
			y -= 1;
			break;
		case 3:
			y += 1;
			break;
		case 4:
			break;
		}
	}
	vector<int> mov; 
	mov.push_back(x);
	mov.push_back(y);
	return mov;
}

void Cell::updatePredMass(float amount)
{
	if (predMass > 20) {
		predMass = 20;
	}
	else {
		predMass += (amount);
		if (predMass > 20) {
			predMass = 20;
		}
	}
}

float Cell::updatePreyMass(float amount, float maxMass)
{
	float newPreyMass = preyMass + amount;
	if (newPreyMass > maxMass) {
		preyMass = maxMass;
		return (newPreyMass -maxMass);
	}
	else {
		preyMass = newPreyMass;
		return 0.0; 
	}
}

void Cell::updateMass(float amount)
{
	if (mass > 1) {
		mass = 1;
	}
	else {
		mass += (amount);
		if (mass > 1) {
			mass = 1;
		}
	}
}

void Cell::createCellGene(int length) {
	cellGene.clear();
	for (int i = 0; i < length; i++) {
		if (i < length / 2) {
			cellGene.push_back(1);
		}
		else {
			cellGene.push_back(0);
		}
	}
}

float Cell::evaluateGene() {
	float f = 0;
	int layers = 0;
	int gs = cellGene.size();
	if (gs == 2) {
		layers = 1;
	}
	else if (gs == 4) {
		layers = 2;
	}
	else if (gs == 8) {
		layers = 3;
	}
	else if (gs == 16) {
		layers = 4;
	}
	else if (gs == 32) {
		layers = 5;
	}
	else if (gs == 64) {
		layers = 6;
	}
	else if (gs == 128) {
		layers = 7;
	}

	int currentLayer = 0;
	while (currentLayer < layers) {
		currentLayer += 1;
		vector<vector<int>> compareArray;
		int bitlength = pow(2, (currentLayer - 1));
		int i = 0;
		while (i < (cellGene.size() - bitlength + 1)) {
			vector<int> cBitString;
			for (int j = 0; j < bitlength; j++) {
				cBitString.push_back(cellGene[j + i]);
			}
			i += bitlength;
			compareArray.push_back(cBitString);
		}
		int amountChecks = (int)compareArray.size() / 2;
		i = 0;
		while (i < (compareArray.size() - 1)) {
			i += 1;
			if (compareArray[i] == compareArray[i - 1]) {
				f += pow(2, (currentLayer - 1));
			}
			i += 1;
		}
	}
	return f;
}