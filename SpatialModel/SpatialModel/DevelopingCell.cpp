#include "DevelopingCell.h"
#include <iostream>

DevelopingCell::DevelopingCell(float mr, float dr, float rr, int const& tp, int devSwaps) : Cell(mr, dr, rr, tp, devSwaps)
{
	mutationRate = mr;
	deathRate = dr;
	reproductionRate = rr;
	type = tp;
	for (int m = 0; m < devSwaps; m++) {
		developTime.push_back(100);
	}
}

DevelopingCell::~DevelopingCell()
{
}

void DevelopingCell::mutate(float mr)
{


	Cell::mutate(mr);
	for (int n = 0; n < genes.size(); n++) {
		for (int i = 0; i < genes[0].size(); i++) {
			if (randFloat(0.0, 1.0) < mr) {
				genes[n][i] = rand() % 2;
			}
		}
	}
	for (int n = 0; n < genes[0].size(); n++) {
		cellGene[n] = genes[0][n];
	}	
}

void DevelopingCell::mutDevelopmentTime(float mr) {
	int minDevTime = 10;
	int maxDevTime = 1000;
	for (int m = 0; m < developTime.size();m++) {
		if (randFloat(0.0, 1.0) < mr) {
			developTime[m] += (rand() % (2 * devTimeMR)) - devTimeMR;
			if (developTime[m] < minDevTime) {
				developTime[m] = minDevTime;
			}
			if (developTime[m] > maxDevTime) {
				developTime[m] = maxDevTime;
			}
		}
	}
}

void DevelopingCell::developGene(int geneNr) {
	for (int i = 0; i < cellGene.size(); i++) {
		int value = cellGene[i] + genes[geneNr][i];
		if (value == 2) {
			cellGene[i] = false;
		}
		else if (value == 1) {
			cellGene[i] = true;
		}
		else {
			cellGene[i] = false;
		}
	}
	//fitness = evaluateGene();
	/*cout << "deveGene:";
	for (int i = 0; i < developedGene.size(); i++) {
		 cout << developedGene[i];
	} 
	cout << endl;
	cout << "cellGene: ";
	for (int i = 0; i < developedGene.size(); i++) {
		cout << cellGene[i];
	}
	cout << endl;
	cout << "cdevGene: ";
	for (int i = 0; i < developedGene.size(); i++) {
		cout << cellDevelopGene[i];
	}
	cout << endl;
	cellGene = developedGene;*/
}

void DevelopingCell::initializeGene(int geneLength, int swaps)
{
	cellGene.clear();
	genes.clear(); // can stay empty
	for (int n = 0; n < swaps + 1; n++) {
		vector<bool> gene;
		for (int i = 0; i < geneLength; i++) {
			gene.push_back(rand() & 1);
		}
		genes.push_back(gene);
	}
	for (int i = 0; i < geneLength; i++) {
		cellGene.push_back(genes[0][i]);
	}
	cStage = 0;
	age = 0;
}

void DevelopingCell::createCellGene(int length, int swaps) {
	cellGene.clear();
	genes.clear();
	vector <bool> gene;
	for (int i = 0; i < length; i++) {
		if (i < length / 2) {
			gene.push_back(1);
		}
		else {
			gene.push_back(0);
		}
	}
	genes.push_back(gene);
	for (int n = 0;n < swaps; n++) {
		vector <bool> gene;
		for (int i = 0; i < length; i++) {
			gene.push_back(0);
		}
		genes.push_back(gene);
	}

	for (int i = 0; i < length; i++) {
		cellGene.push_back(genes[0][i]);
	}
	cStage = 0;
	age = 0;
}
