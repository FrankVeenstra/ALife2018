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
	if (randFloat(1, 0) < mr) {
		immortal = rand() % 2;
	}
}

void Cell::mutDeathRate(float mr)
{
	// mutating deathRate..
	if (randFloat(1, 0) < mr) {
		deathRate += randFloat(-deathRate,deathRate);
		if (deathRate > 1.0) {
			deathRate == 1.0;
		}
		else if (deathRate < 0.000000001) {
			deathRate == 0.000000001;
		}
	}
	if (randFloat(1, 0) < mr) {
		deathRate += randFloat(0,0.1);
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
	vector<int> dTime, int age_v, int currentStage, int da, bool im)
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
	genes = gs;
	dAge = da; 
	immortal = im;
	//cout << &gs << endl;
	//cout << &genes << endl;
	//if (gs.size() > 0) {
	//	cout << &gs[0] << endl;
	//	cout << &gs[0][0] << endl;
	//	cout << &genes[0] << endl;
	//	cout << &genes[0][0] << endl;
	//}
	//
	//cout << &da << endl;
	//cout << &dAge << endl;

	//if (gs.size() > 0) {
	//	cout << "mje" << endl;
	//}
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

float Cell::updatePreyMass(float amount)
{
	//	if (preyMass > 10) {
	//		preyMass = 10;
	//	}
	//	else {
	float newPreyMass = preyMass + amount;
	if (newPreyMass > 10) {
		preyMass = 10;
		return (newPreyMass -10);
	}
	else {
		preyMass = newPreyMass;
	}
	//	}
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
//		cellGene.push_back(rand() % 2);
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


	//for (int i = 0; i < 10; i++) { // cannot compute layers larger than 10 (1024 bit)
	//	i += 1;
	//	float l = pow(2, i);
	//	if (pow(2, i) == cellGene.size()) {
	//		layers = i;
	//		break;
	//	}
	//	if (layers == 0) {
	//		break;
	//	}
	//}
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
//	if (f == 32) {
//		f = 64;
//	}
	return f;
}

//void Cell::setHealth(int const& health)
//{
//	m_health = clamp(health,0,m_maxHealth) ;
//}
//void Cell::setMaxHealth(int const& maxHealth)
//{
//	m_maxHealth = maxHealth ;
//}
//void Cell::setForce(int const& force)
//{
//	m_force = clamp(force,0,m_maxForce) ;
//}
//void Cell::setMaxForce(int const& maxForce)
//{
//	m_maxForce = maxForce ;
//}
//void Cell::setTeam(int const& team)
//{
//	m_team = team ;
//}
//int Cell::getHealth() const
//{
//	return m_health ;
//}
//int Cell::getMaxHealth() const
//{
//	return m_maxHealth ;
//}
//int Cell::getForce() const
//{
//	return m_force ;
//}
//int Cell::getMaxForce() const
//{
//	return m_maxForce ;
//}
//int Cell::getTeam() const
//{
//	return m_team ;
//}
