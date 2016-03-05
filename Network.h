/* 
 * File:   Network.h
 * Author: qfeuille
 *
 * Created on 26 May 2012, 00:16
 */

#ifndef NETWORK_H
#define	NETWORK_H


#include<vector>
#include "mtrand.h"
#include <functional>
namespace NumMethod{
    
template<typename NodeT>   
class Node;
class BooleanNetwork;


class ConnectionFunction {
protected:
	int N; int n; std::vector <int> Coordinates;
	ConnectionFunction (int iN);
public:
	int getN (){return N;}
        int mod (int x, int modulo);
	//virtual void getConnections (int in, NodeIt Start, std::vector<NodeIt>* oConnections)=0;
};

inline int ConnectionFunction::mod (int x, int modulo){
    return (x%modulo + modulo)%modulo;
}

template<typename NodeT>
class Node {
protected: const unsigned int id;
          typedef typename std::vector<NodeT>::iterator NodeIt;
          typedef typename std::vector<NodeIt>::iterator NodeItIt;

public:	std::vector <NodeIt> Connections;
	Node(int);
	Node(int, std::vector<NodeIt>* iConnections);
	void setConnections (std::vector<NodeIt>* iConnections);
	NodeIt getConnection (int numConnection);
        const unsigned int getid (){return id;}
};



class BooleanNode : public Node<BooleanNode>{
	bool State;
public:	
	BooleanNode(int iid):Node<BooleanNode>(iid){};
	BooleanNode(int iid, std::vector<NodeIt>* iConnections): Node<BooleanNode>(iid, iConnections){};
	void setState(bool istate);
	bool getState();
	void flipState();
	int sumXOR ();
	
};


inline void BooleanNode::setState (bool iState){
	State=iState;}

inline void BooleanNode::flipState (){
	State=!State;}

inline bool BooleanNode::getState (){
	return State;}


inline int BooleanNode::sumXOR(){
	int Sum=0;
	for ( NodeItIt itConnection=Connections.begin() ; itConnection < Connections.end(); ++itConnection ){
			if ( ((*itConnection)->getState()) ^ State)
				Sum++;
	}
	return Sum;
}



template<typename NodeT>
class Network {

protected:
   
    std::vector <NodeT> Nodes;
           
public:
     typedef typename std::vector<NodeT>::iterator NodeIt;
    typedef typename std::vector<NodeIt>::iterator NodeItIt;
	Network (){};
        template<typename FunctPtr>
        Network (int in, FunctPtr&);
	~Network(){};
        int addNodes (int in);
        NodeIt node (int in) {return Nodes.begin()+ in;}
         template<typename FunctPtr>
	int setNodeConnections (FunctPtr&);
	int getNumConnections (int n);
        int getSize () {return Nodes.size();};
	

};

template<typename NodeT>
 template<typename FunctPtr>
Network<NodeT>::Network(int in, FunctPtr& setter){
    addNodes (in);
    setNodeConnections (setter); 
}

template<typename NodeT>
 template<typename FunctPtr>
int Network<NodeT>::setNodeConnections (FunctPtr& setter){
	int NumConnections=0;
	

	for (unsigned int n=0; n<Nodes.size(); n++ ){
		std::vector <NodeIt> Connections;
                std::vector <int> out;
                setter (n,out);
                for (unsigned int i=0;i<out.size();i++){
                    Connections.push_back(Nodes.begin() + out[i]);
                }
		
		Nodes [n].setConnections (&Connections);
		NumConnections+=Connections.size ();
	}
	return (NumConnections);
}


 
class BooleanNetwork : public Network<BooleanNode> {

	MTRand * drand; 
public:
        BooleanNetwork ();
        template<typename FunctPtr>
	 BooleanNetwork (int in, FunctPtr& setter) : Network<BooleanNode> (in,setter){drand =new MTRand ((unsigned) time (0));};
	~BooleanNetwork();
	int sumXOR ();
	int sumNodeXOR (int n);
        template<typename FunctPtr>
        void setNodes (FunctPtr&& setter);
	void randomiseStates ();
	void flipNode (int n);
	

};


 template <typename FunctPtr >
void BooleanNetwork::setNodes (FunctPtr&& setter){
    for (unsigned int n=0; n<Nodes.size(); n++ ){
        Nodes [n].setState (setter(n));
	}
}



inline void BooleanNetwork::flipNode (int n){
	Nodes [n].flipState();}


	


inline int BooleanNetwork::sumNodeXOR(int n){
	return Nodes[n].sumXOR ();
}



ConnectionFunction::ConnectionFunction (int iN){
	N=iN;

}



template<typename NodeT>
Node<NodeT>::Node (int iid, std::vector<NodeIt>* iConnections): id(iid){
		setConnections (iConnections);
	}
	template<typename NodeT>
Node<NodeT>::Node (int iid): id(iid){}


int BooleanNetwork::sumXOR (){
	NodeIt itNode;
	NodeItIt itConnection;
	int Sum=0;
	for ( itNode=Nodes.begin() ; itNode < Nodes.end(); itNode++ )
		Sum+= itNode->sumXOR();
	return (int) (Sum/2);
}


template<typename NodeT>
void Node<NodeT>::setConnections(std::vector<NodeIt>* iConnections){
		Connections=*iConnections;
}


template<typename NodeT>
typename std::vector<NodeT>::iterator Node<NodeT>::getConnection (int numConnection){
	return Connections [numConnection];
}

template<typename NodeT>
int Network<NodeT>::addNodes (int in){
    for (int n=0; n<in; n++ ){
		Nodes.push_back(NodeT (n));}
    return Nodes.size();    
}



void BooleanNetwork::randomiseStates (){
	for ( NodeIt itNode=Nodes.begin() ; itNode < Nodes.end(); itNode++ ){
		double random=(*drand)();
		if (random>0.5) itNode->setState (true);
		else itNode->setState(false);}
	
	
}

BooleanNetwork::~BooleanNetwork(){
		delete drand;}

BooleanNetwork::BooleanNetwork (){
		drand =new MTRand ((unsigned) time (0));
	}

}

#endif	/* NETWORK_H */

