#ifndef interval_tree_H
#define interval_tree_H
#include "read.h"
#include <vector>
#include <string>
using namespace std;
class node{
public:
	int start;
	int stop; 
	vector<interval> INFO;
	node * next;
	node(int, int, std::vector<interval>);
};
class T{
public:
	T();
	T(node *);
	T(interval *);
	node * nodes;
	T * left = NULL;
	T * right= NULL;
	void assemble(interval *);	
	void build();
	void traverse();
	vector<interval> search_interval(int, int);
	vector<interval> search_point(int);
};


#endif