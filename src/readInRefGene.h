/*
 * readInRefGene.h
 *
 *  Created on: Jan 25, 2014
 *      Author: joeyazo
 */

#ifndef READINREFGENE_H_
#define READINREFGENE_H_
#include <string>
#include <iterator>
#include "split.h"
using namespace std;
class BT{
public:
	int start, stop;
	BT * left;
	string GENE_ID;
	vector<int> stop_overlaps;
	BT * right;
	void insert(int, int, string);
	BT(int val1, int val2, string val3){
		start = val1;
		stop = val2;
		GENE_ID = val3;
		left = NULL;
		right = NULL;
	}
	string search(int);
	void clear(BT *);

};

string BT::search(int query_pos){
	if(query_pos >= start && query_pos <= stop){
		return GENE_ID;
	}
	else if (query_pos < start){
		if (left){
			return left->search(query_pos);
		}
		else{
			return "";
		}
	}

	else if(query_pos > start){
		if (right){
			return right->search(query_pos);
		}
		else{
			return "";
		}
	}
	else{
		return "";
	}
}

void BT::clear(BT * node){
	if (node->left != NULL){
		clear(node->left);
	}
	if (node->right != NULL){
		clear(node->right);
	}
	delete(node);

}



void BT::insert(int ins_start, int ins_stop,string gene){
	if (ins_start < start){
		if(left){
			left->insert(ins_start, ins_stop, gene);
		}
		else{
			left = new BT(ins_start, ins_stop, gene);
		}
	}
	else if (ins_start > start){
		if (right){
			right->insert(ins_start, ins_stop, gene);
		}
		else{
			right = new BT(ins_start, ins_stop, gene);
		}

	} else{
	}



}




class chromosome{
	//Goal of this is class is be a root of a Binary Tree for the interval
public:
	string ID;
	BT * interval;
	chromosome*next;
	chromosome(string value, int start, int stop,string gene){
		ID = value;
		next = NULL;
		interval = new BT(start, stop,gene);
	}
	void insert(string , int , int,string);
	void clear_chromosome(chromosome*);

};

void chromosome::insert(string ins_chrom, int start, int stop, string gene_ID){
	if (ins_chrom == ID){
		interval->insert(start, stop, gene_ID);
	}
	else if (next != NULL){
		if ( next -> ID == ins_chrom){
			next->interval->insert(start, stop, gene_ID);
		}
		else{
			chromosome * temp_ptr = next;
			while (temp_ptr->next){
				if (temp_ptr->ID == ins_chrom ){
					break;
				}
				temp_ptr=temp_ptr->next;
			}
			if (temp_ptr->ID == ins_chrom){
				temp_ptr->interval->insert(start, stop, gene_ID);
			}
			else if (not temp_ptr->next)
			{
				temp_ptr->next = new chromosome(ins_chrom, start, stop, gene_ID);
				temp_ptr->next->interval->insert(start, stop, gene_ID);
			}
		}
	}

	else{
		next = new chromosome(ins_chrom, start, stop, gene_ID);
		next->interval->insert(start, stop, gene_ID);
	}
}


void chromosome::clear_chromosome(chromosome*node){
	if (node->next != NULL){
		clear_chromosome(node->next);
	}
	if (node->interval != NULL){
		node->interval->clear(node->interval);
	}
	delete(node);
}
class strand{
	//Goal of this class is store a vector of chromosome classes
public:
	string ID;
	chromosome*chromosome_ptr;
	strand *next;
	strand(){
		ID = "";
		chromosome_ptr = NULL;
		next = NULL;
	}

	strand(string value){
		ID = value;
		chromosome_ptr = NULL;
		next = NULL;

	}
	void clear(strand*);
	void insert(string , string, int, int,string);
};

void strand::insert(string ins_strand, string ins_chromosome, int start, int stop, string ins_gene_ID){
	if (ins_strand == ID){
		if (chromosome_ptr == NULL){
			chromosome_ptr = new chromosome(ins_chromosome, start, stop, ins_gene_ID);
		}

		chromosome_ptr->insert(ins_chromosome, start, stop, ins_gene_ID);
	}
	else{
		if (next){
			if (next->ID == ins_strand){
				next->chromosome_ptr->insert(ins_chromosome, start, stop, ins_gene_ID);
			}
		}
		else{
			next  = new strand(ins_strand);
			next->chromosome_ptr = new chromosome(ins_chromosome, start, stop, ins_gene_ID);
			next->chromosome_ptr->insert(ins_chromosome, start, stop, ins_gene_ID);
		}


	}

}

void strand::clear(strand*node){
	if (node->next != NULL){
		clear(node->next);
	}
	if (node->chromosome_ptr != NULL){
		node->chromosome_ptr->clear_chromosome(node->chromosome_ptr);
	}
	delete(node);
}

class tree{
public:
	string ID;
	strand * strand_ptr;
	tree(){
		ID="ROOT";
		strand_ptr = NULL;
	}
	void insert(string ins_strand, string ins_chromosome, int start, int stop, string ins_gene_ID){
		if (strand_ptr != NULL){
			strand_ptr->insert(ins_strand, ins_chromosome, start, stop, ins_gene_ID);
		}
		else{
			strand_ptr = new strand(ins_strand);
			strand_ptr->insert(ins_strand, ins_chromosome, start, stop, ins_gene_ID);
		}
	}
	string search(string, string,  int);
	void clear(tree *);
};


void tree::clear(tree * node){
	if (strand_ptr!= NULL){
		strand_ptr->clear(strand_ptr);
	}
	delete(node);
}




string search(tree * root, string query_strand , string query_chrom, int query_pos){
	strand * ptr = root->strand_ptr;
	if (ptr == NULL)
	{
		cout<<"Reference Gene Data Structure was not initialized..."<<endl;
		return "";
	}

	chromosome * chrom_ptr = NULL;
	BT * node_ptr = NULL;

	while(ptr != NULL){
		if (ptr->ID == query_strand){
			chrom_ptr = ptr->chromosome_ptr;
			break;
		}
		ptr=ptr->next;
	}
	if (chrom_ptr == NULL){
		return "";
	}

	while(chrom_ptr != NULL){
		if (chrom_ptr->ID == query_chrom){
			node_ptr=chrom_ptr->interval;
			break;
		}
		chrom_ptr=chrom_ptr->next;
	}
	if (node_ptr == NULL){
		return "";
	}
	string result = node_ptr->search(query_pos);
	return result;
}












tree * readRefGene(string refGeneFileName){
	ifstream myfile(refGeneFileName);

	tree *root = new tree();
	int pos;
	string line;
	string token;
	string gene;
	string chrom, strand;
	int start, stop;
	string item;
	vector<string> lineInfo;
	int i = 0;
	if (myfile.is_open()){
		while ( getline(myfile,line) ){

			lineInfo=split(line, "\t");

			//===============================================//
			// Initialize Information from RefGene.Bed file	 //
			// like chromosome, start, stop, strand etc.	 //
			//===============================================//
			chrom = lineInfo[2];
			strand = lineInfo[3];
			start = atoi(lineInfo[4].c_str());
			stop = atoi(lineInfo[5].c_str());
			gene = lineInfo[1];
			root->insert(strand, chrom , start, stop,gene);
			lineInfo.clear();
		}
	}
	else{
		cout<<"Unable to open reference Gene File: "<<refGeneFileName<<"...try specifying full path/spelling error correction"<<endl;
		return root;
	}
	//=============================================================
	// Some Unit Test Cases
		return root;
}


bool refGeneTest(string refGeneFileName){
	tree *root = readRefGene(refGeneFileName);
	cout<<"\nRunning Unit Tests on Reference Gene Data Structure..."<<endl;
	int correct = 0;
	string FOUND = search(root, "-", "chr16", 20370491);
	if (FOUND.empty()){
		cout<<"Test1: **FAILED**"<< endl;
	}
	else{
		cout<<"Test1: **PASSED**"<< endl;
		correct++;
	}
	FOUND = search(root, "-", "chr16", 6);
	if (FOUND.empty()){
		cout<<"Test2: **PASSED**"<< endl;
		correct++;
	}
	else{
		cout<<"Test2: **FAILED**" << endl;
	}
	FOUND = search(root, "+", "chr10", 64893011);
	if (not FOUND.empty()){
		cout<<"Test3: **PASSED**"<<endl;
		correct++;
	}
	else{
		cout<<"Test3: **FAILED**"<<endl;
	}
	if (correct!= 3){
		return 0;
	}
	return 1;
};
#endif /* READINREFGENE_H_ */
