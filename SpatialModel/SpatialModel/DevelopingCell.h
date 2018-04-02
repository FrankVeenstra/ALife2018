#pragma once
#include "Cell.h"
class DevelopingCell :
	public Cell
{
public:
	DevelopingCell(float mr, float dr, float rr, int const& tp, int devSwaps);
	~DevelopingCell();
	/*int developTime = 20;
	int age = 0;
	vector<int> cellDevelopGene;*/
	void mutate(float mr);
	void mutDevelopmentTime(float mr);
	void developGene(int geneNR);
	void initializeGene(int geneLength, int swaps);
	void createCellGene(int length, int swaps);
	int devTimeMR = 100;
};

