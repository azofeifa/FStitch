#include "interval_tree.h"
#include "read.h"
#include <iostream>

using namespace std;
void sort_intervals(interval * intervals){
	bool SWAPED 	= 1;
	interval * C 	= intervals;
	interval * root = intervals;

	int save_start, save_stop;
	string save_info;
	while (SWAPED){
		SWAPED 	= 0;
		C 		= root;
		while (C->next != NULL){
			if (C->next->start < C->start){
				save_start 		= C->start;
				save_stop 		= C->stop;
				save_info 		= C->info;
				C->start 		= C->next->start;
				C->stop 		= C->next->stop;
				C->info 		= C->next->info;
				C->next->start 	= save_start;
				C->next->stop 	= save_stop;
				C->next->info 	= save_info;
				SWAPED 			= 1;
			}
			C 	= C->next;
		}
	}
}
node::node(int st, int sp, vector<interval> INFOS){
	start 	= st;
	stop 	= sp;
	INFO 	= INFOS;
	next 	= NULL;
}
T::T(){};
void T::assemble(interval * intervals){
	interval * C 	= intervals;
	node * root 	= NULL;
	node * N 		= NULL;
	int o_st, o_sp ;
	vector<interval> INFOS;
	while (C){
		o_st 	= C->start, o_sp 	= C->stop;
		INFOS.push_back(*C);
		while(C->next != NULL and C->next->start < o_sp and C->next->stop >o_st){
			if (C->next->start < o_st){
				o_st 	= C->next->start;
			}
			if (C->next->stop > o_sp){
				o_sp 	= C->stop;
			}
			C 	= C->next;
			INFOS.push_back(*C);
		}
		if (N== NULL){
			N 			= new node(o_st, o_sp, INFOS);	
			root 		= N;
		}else{
			N->next 	= new node(o_st, o_sp, INFOS);
			N 			= N->next;
		}
		INFOS.clear();
		C 		= C->next;
		
	}
	nodes 	= root;
}
T::T(node * NODES){
	nodes 	= NODES;
}
void T::build(){
	node * C 		= nodes;
	node * root 	= NULL;
	node * LEFT 	= NULL;
	if (nodes != NULL){
		//find middle node
		node * one 		= nodes;
		node * two 		= nodes; 
		while (one != NULL){
			one 		= one->next;
			if (one != NULL){
				one 	= one->next;
			}
			if (LEFT == NULL and one !=NULL){
				LEFT 		= new node(two->start, two->stop, two->INFO);
				root 		= LEFT;
			}else if (one != NULL and LEFT!=NULL){
				LEFT->next 	= new node(two->start, two->stop, two->INFO);
				LEFT 		= LEFT->next;
			}
			two 		= two->next;
		}
		//so two is the middle now
		if (two!=NULL){
			if (two->next){
				right 			= new T(two->next);
				right->build();
			}
			if (LEFT){
				left 			= new T(root);
				left->build();
			}
			two->next 		= NULL;
			nodes 			= two;
		}
	}
}
void T::traverse(){
	if (nodes->next !=NULL){
		cout<<"WHAT?"<<endl;
	}
	if (left){
		left->traverse();
	}
	if (right){
		right->traverse();
	}
}


vector<interval> T::search_interval(int st, int sp){
	if (sp < nodes->start and left != NULL){
		return left->search_interval(st, sp);
	}else if(st > nodes->stop and right != NULL){
		return right->search_interval(st, sp);
	}
	vector<interval> FOUND;
	
	for (int i = 0; i < nodes->INFO.size(); i++){
		if ((nodes->INFO[i].start <st && st <nodes->INFO[i].stop) || (nodes->INFO[i].start <sp && sp <nodes->INFO[i].stop)
			||  (st < nodes->INFO[i].start and sp > nodes->INFO[i].stop) ){
			FOUND.push_back(nodes->INFO[i]);
		}
	}
	return FOUND;
}
vector<interval> T::search_point(int pt){
	if (pt < nodes->start and left != NULL){
		return left->search_point(pt);
	}else if(pt > nodes->stop and right != NULL){
		return right->search_point(pt);
	}
	vector<interval> FOUND;
	
	for (int i = 0; i < nodes->INFO.size(); i++){
		if (nodes->INFO[i].start <pt && pt <nodes->INFO[i].stop ){
			FOUND.push_back(nodes->INFO[i]);
		}
	}
	return FOUND;
}
T::T(interval * intervals){
	sort_intervals(intervals);

	assemble(intervals);
	build();	
}
