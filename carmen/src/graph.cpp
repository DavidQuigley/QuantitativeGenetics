#include <fstream>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <string>
#include <boost/unordered_map.hpp>
#include <boost/random.hpp>
#include "DataStructures.h"
#include "graph.h"


bool Graph::has_node_label(int node_a){
	return this->idx2str.find(node_a) != this->idx2str.end();
}

std::string Graph::get_node_label_by_index(int node_a){
	if( this->idx2str.find(node_a) == this->idx2str.end() )
		throw std::string("Node not found in get_node_label_by_index");
	else
		return this->idx2str[node_a];
}

int Graph::get_node_index_by_label(std::string label){
	if( this->str2idx.find(label) == this->str2idx.end() )
		return -1;
	else
		return this->str2idx[label];
}

int Graph::node_type(int node_a){
	if( this->node_types.find(node_a) == this->node_types.end() )
		return -1;
	else
		return this->node_types[node_a];
}

void Graph::set_node_label(int node_a, std::string label){
	if( this->adj.find(node_a) == this->adj.end() )
		throw std::string("node not found in call to set_node_label");
	else{
		this->idx2str[node_a] = label;
		this->str2idx[label] = node_a;
	}
}

void Graph::set_node_type(int node_a, int node_type ){
	if( this->adj.find(node_a) == this->adj.end() )
		throw std::string("node not found in call to set_node_type");
	else{
		this->node_types[node_a] = node_type;
	}
}

bool Graph::has_locus_neighbor( int node_a ){
	if( adj.find(node_a) != adj.end() ){
		for(boost::unordered_map<int, double>::iterator node_itr=adj[node_a]->begin(); node_itr != adj[node_a]->end(); node_itr++){
			if( this->node_types[node_itr->first] == NODE_LOCUS ){
				return true;
			}
		}
	}
	return false;
}


void Graph::add_edge( int node_a, int node_b ){
    add_edge(node_a, node_b, 0, NODE_GENE, NODE_GENE);
}

void Graph::add_edge(int node_a, int node_b, double weight){
	add_edge(node_a, node_b, weight, NODE_GENE, NODE_GENE);
}

void Graph::add_edge(int node_a, int node_b, double weight, int node_a_type, int node_b_type){
	if( adj.find(node_a) == adj.end() )
		adj[node_a] = new boost::unordered_map<int, double>();
	(*adj[node_a])[node_b] = weight;
	if( adj.find(node_b) == adj.end() )
		adj[node_b] = new boost::unordered_map<int, double>();
	(*adj[node_b])[node_a] = weight;
	node_types[node_a] = node_a_type;
	node_types[node_b] = node_b_type;
}

void Graph::remove_edge(int node_a, int node_b){
	adj[node_a]->erase(node_b);
	adj[node_b]->erase(node_a);
	if( adj[node_a]->empty() ){
		adj.erase(node_a);
		node_types.erase(node_a);
		if( this->has_node_label(node_a) ){
			str2idx.erase( this->get_node_label_by_index(node_a) );
			idx2str.erase(node_a);
		}
	}
	if( adj[node_b]->empty() ){
		adj.erase(node_b);
		node_types.erase(node_b);
		if( this->has_node_label(node_b) ){
			str2idx.erase( this->get_node_label_by_index(node_b) );
			idx2str.erase(node_b);
		}
	}
}


void Graph::print(){
	std::vector<std::vector<int>*>* e = new std::vector<std::vector<int>*>();
	this->edges(e);
	int g1, g2;
	std::string lbl1, lbl2;
	for(int i=0; i<(int)e->size(); i++){
		g1 = e->at(i)->at(0);
		g2 = e->at(i)->at(1);
		lbl1 = this->idx2str[g1];
		lbl2 = this->idx2str[g2];
		std::cout << g1 << "," << g2 << " (" << lbl1 << ", " << lbl2 << ")\n";
		delete e->at(i);
	}
}

void Graph::print_nodes(){
    for( graph_itr itr = adj.begin(); itr != adj.end(); itr++){
        std::cout << itr->first;
        if( this->has_node_label(itr->first))
        	std::cout << " " << this->get_node_label_by_index(itr->first);
        std::cout << "\n";
    }
}

void Graph::copy( Graph* g ){
	this->clear();
	for (boost::unordered_map <int, std::string>::iterator i=g->idx2str.begin();
		                                         i != g->idx2str.end(); ++i){
		this->idx2str[ i->first ] = i->second;
		this->str2idx[ i->second ] = i->first;
	}
	std::vector<std::vector<int>*>* e = new std::vector<std::vector<int>*>();
	g->edges(e);
	int g1, g2;
	// copy all idx2str entries
	for(int i=0; i<(int)e->size(); i++){
		g1 = e->at(i)->at(0);
		g2 = e->at(i)->at(1);
		this->add_edge(g1, g2, g->edge_weight(g1, g2));
		delete e->at(i);
	}
	std::vector<int> source_nodes;
	g->nodes(source_nodes);
	for(int i=0; i<int(source_nodes.size()); i++)
		g->set_node_label(source_nodes.at(i), g->get_node_label_by_index( source_nodes.at(i) ) );

}


Graph::Graph(){
    this->clique_graph = NULL;
    this->verbose = false;
}

void Graph::add_qtls(std::vector<QTL*>& qtls){
	QTL* Q;
	int node_a, node_b;
	int cur_idx = 0;
	for(int i=0; i<int(qtls.size()); i++){
		Q = qtls.at(i);
		std::string LID( Q->locus_id );
		if( this->str2idx.find(LID) != this->str2idx.end() )
			node_a = this->str2idx[LID];
		else{
			node_a = cur_idx;
			cur_idx += 1;
			this->str2idx[LID] = node_a;
			this->idx2str[node_a] = LID;
		}
		std::string PID( Q->probe_name );
		if( this->str2idx.find(PID) != this->str2idx.end() )
			node_b = this->str2idx[PID];
		else{
			node_b = cur_idx;
			cur_idx += 1;
			this->str2idx[PID] = node_b;
			this->idx2str[node_b] = PID;
		}
		this->add_edge( node_a, node_b, Q->pval, NODE_LOCUS, NODE_GENE);
	}
}


Graph::~Graph(){
	this->clear();
}

void Graph::set_verbose(bool verbose){
	this->verbose = verbose;
}

void Graph::clear(){
    delete this->clique_graph;
    this->idx2str.clear();
    this->str2idx.clear();
    this->node_types.clear();
    for( graph_itr itr = adj.begin(); itr != adj.end(); itr++){
        delete itr->second;
    }
    adj.clear();
}

bool Graph::has_node(int node_a){
	return adj.find(node_a) != adj.end();
}


bool Graph::has_edge(int node_a, int node_b, double& weight){
	if( adj.find(node_a) != adj.end() ){
		if( adj[node_a]->find( node_b ) != adj[node_a]->end() ){
			weight = (*adj[node_a])[node_b];
			return true;
		}
	}
	return false;
}

double Graph::edge_weight(int node_a, int node_b){
	if( adj.find(node_a) == adj.end() || adj[node_a]->find( node_b ) == adj[node_a]->end() )
		throw std::string("edge not found");
	else
		return (*adj[node_a])[node_b];
}


int Graph::order(){
    return (int)adj.size();
}

bool comp_degrees(const std::pair<int, int>*a, const std::pair<int, int>*b){
	return a->second > b->second;
}

void Graph::hubs(std::vector<int>& h, int percentile){

    h.clear();
    std::vector<std::pair<int, int>*> degrees;
    std::pair<int, int>* p;
    int a;
    for( graph_itr itr = adj.begin(); itr != adj.end(); itr++){
        a = itr->first;
        p=new std::pair<int, int>();
        p->first =a;
        p->second = degree(a);
        degrees.push_back( p );
    }
    if( degrees.size()==0 )
        return;
    std::sort(degrees.begin(), degrees.end(), comp_degrees );
    int cut_off = (int) ( (float)degrees.size() *  ( (float)percentile / 100.0 ) );
    if( cut_off==0 )
        cut_off=1;
    for(int i=0; i<cut_off; i++){
        h.push_back( degrees.at(i)->first );
    }
    for(int i=0; i<(int)degrees.size(); i++){
        delete degrees.at(i);
    }
}

int find_dupe(std::vector<int>* V){
    for(int i=1; i<(int)V->size(); i++){
        if(V->at(i-1)==V->at(i))
            return V->at(i);
    }
    return -1;
}

void Graph::diff( std::vector<int>* A, std::vector<int>* B, std::vector<int>* C){
    // assumes sorted.  Place into C anything found in B but not A
    C->clear();
    int a=0, b=0;
    int n_a = (int)A->size();
    int n_b = (int)B->size();
    int d = find_dupe(A);
    if(d!=-1)
        std::cout << "dupe in A incoming\n";

    if( n_b == 0 )
        return;
    if(n_a==0 && n_b != 0){ // A is empty
        for(int i=0; i<n_b; i++){ C->push_back( B->at(i) ); }
        return;
    }
    while(a<n_a && b<n_b){
        if(A->at(a)==B->at(b)){
            ++a; ++b;
        }
        else if( A->at(a) > B->at(b) ){
            C->push_back(B->at(b) );
            ++b;
        }
        else{
            ++a;
        }
    }
    while(b<n_b){
        C->push_back(B->at(b++));
    }
}

int Graph::n_diff( std::vector<int>* A, std::vector<int>* B){
    // assumes sorted
    int a=0, b=0, n_diff=0;
    int n_a = (int)A->size();
    int n_b = (int)B->size();

    if(n_a==0 && n_b != 0)
        return n_b;
    else if(n_b==0 && n_a != 0 )
        return n_a;

    while(a<n_a && b<n_b){
        if(A->at(a)==B->at(b)){
            ++a;
            ++b;
        }
        else if( A->at(a)< B->at(b) ){
            ++n_diff;
            ++a;
        }
        else{
            ++n_diff;
            ++b;
        }
    }
    if( a != b ){
        if(b>a)
            return n_diff + b-a;
        else
            return n_diff + a-b;
    }
    else{
        return n_diff;
    }
}

void Graph::subgraph( Graph& g, std::vector<int>& ids ){
    // copy all edges in ids and node properties of those ids into g
	double w;
    g.clear();
    int s = (int)ids.size();

    for(int i=0; i<s; i++){
        if( !has_node(ids.at(i)))
        	throw std::string( "Requested id for subgraph does not exist in graph" );
        for(int j=i+1; j<s; j++){
            if( has_edge(ids.at(i), ids.at(j), w) ){
                g.add_edge(ids.at(i), ids.at(j), w);
            }
        }
        if( g.has_node(ids.at(i)) ){ // checks if i is in ids but not connected to any other id in ids
        	g.set_node_type(ids.at(i), this->node_type(ids.at(i)));
            if( this->has_node_label( ids.at(i) ) )
           		g.set_node_label( ids.at(i), this->get_node_label_by_index( ids.at(i) ) );
        }
    }
}


void Graph::subgraph( std::vector<int>& ids ){
    // mutates this graph; reduces to edges in ids
	boost::unordered_map<int, int> id2id;
	for(int i=0; i<int(ids.size()); i++){
		id2id[ids.at(i)] = 1;
	}
	std::vector<std::vector<int>*>* e = new std::vector<std::vector<int>*>();
	this->edges(e);
	int idx1=0, idx2=0;
	for(int i=0; i<(int)e->size(); i++){
		idx1 = e->at(i)->at(0);
		idx2 = e->at(i)->at(1);
		if(id2id.find(idx1)==id2id.end() || id2id.find(idx2)==id2id.end() )
			this->remove_edge(idx1, idx2);
	}
}

int Graph::degree(int node_a){
    if( adj.find(node_a) == adj.end() )
        return -1;
    else
        return (int)adj[node_a]->size();
}

void Graph::edges(std::vector<std::vector<int>*>* edges){
    int a,b;
    boost::unordered_map<int, double>::iterator edge_itr;
    std::vector<int>* v;
    for( int i=0; i<(int)edges->size(); i++)
    	if( edges->at(i) != NULL )
    		delete edges->at(i);
    edges->clear();

    for( graph_itr itr = adj.begin(); itr != adj.end(); itr++){
        a = itr->first;
        for(edge_itr = itr->second->begin(); edge_itr != itr->second->end(); edge_itr++){
            b = edge_itr->first;
            if( a<b ){
                v = new std::vector<int>();
                v->push_back(a);
                v->push_back(b);
                edges->push_back(v);
            }
        }
    }
}


void Graph::nodes( std::vector<int>& N ){
    N.clear();
    for( graph_itr itr = adj.begin(); itr != adj.end(); itr++)
        N.push_back( itr->first );
}


bool Graph::is_neighbor( int node_a, int node_b ){
    if( adj.find(node_a) != adj.end() && adj[node_a]->find(node_b) != adj[node_a]->end() )
        return true;
    return false;
}

void Graph::shared_neighbors(std::vector<int>* C, std::vector<int>& N){
    // given vector of nodes C, find all nodes not in C that are 
    // a neighbor of each node in C
    N.clear();
    boost::unordered_map<int, int> neighbor_cnt;
    boost::unordered_map<int, int> in_C;
    for(int i=0; i<(int)C->size(); i++)
        in_C[C->at(i)] = 1;
    int new_neighbor;
    std::vector<int> node_neighbors;
    for(int i=0; i<(int)C->size(); i++){
        neighbors(C->at(i), node_neighbors);
        for(int j=0; j<(int)node_neighbors.size(); j++){
            new_neighbor = node_neighbors.at(j);
            if( in_C.find(new_neighbor)==in_C.end() ){
                // new neighbor not in C already
                if( neighbor_cnt.find( new_neighbor ) == neighbor_cnt.end() )
                    neighbor_cnt[new_neighbor] = 1;
                else
                    neighbor_cnt[new_neighbor] = neighbor_cnt[new_neighbor] + 1;
            }
        }
    }
    int n_required = (int)C->size();
    for(boost::unordered_map<int, int>::iterator n_itr = neighbor_cnt.begin(); n_itr != neighbor_cnt.end(); n_itr++){
        if( (*n_itr).second == n_required )
            N.push_back( (*n_itr).first );
    }
}


void Graph::reduce_to_triangles( Graph& g ){
	// for each node c
	//   create subgraph N of c's neighbors.
	//   Any edge in this graph must be a triangle, since both nodes
	//   in the edge connect to c and they also connect to each other.
	std::vector<int> candidates;
	this->nodes(candidates);
	boost::unordered_map<int,int> status;
	for(int i=0; i<(int)candidates.size(); i++)
		status[candidates.at(i)]=0;

	Graph neighbors;
	std::vector<int> N;
	std::vector<std::vector<int>*>* edges = new std::vector<std::vector<int>* >();
	if( this->verbose ){
		std::cout << "MESSAGE: Reducing to 3-cliques with " << candidates.size() << " candidates\n";
		std::cout.flush();
	}
	for(int i=0; i<(int)candidates.size(); i++){
		int c = candidates.at(i);
		if( status[c]==1 )
			continue; // already know this is in a triangle
		this->neighbors(c, N);
		this->subgraph(neighbors, N); // deliberately did not include c
		neighbors.edges( edges );
		for( int e=0; e<(int)edges->size(); e++){
			status[ edges->at(e)->at(0) ] = 1;
			status[ edges->at(e)->at(1) ] = 1;
			status[c] = 1;
		}
		if( this->verbose && i % 1000==0){
			std::cout << "MESSAGE: completed " << i + 1 << " of " << candidates.size() <<  " candidates\n";
			std::cout.flush();
		}
	}
	std::vector<int> triangles;
	for(int i=0; i<(int)candidates.size(); i++){
		if( status[ candidates.at(i) ] == 1 )
			triangles.push_back(candidates.at(i) );
	}
	this->subgraph(g, triangles);
}


void Graph::evaluate_quad( boost::unordered_map<int, int>& status, int c, std::vector<int>& neighbors ){
	// find triangles in neighbors of c
	// for each neighbor subnode
	//    iterate through pairs of subnode's neighbors
	//    if any of these pairs are connected, we've found a triangle.
	Graph sub;
	this->subgraph(sub, neighbors);
	std::vector<int> subnodes;
	std::vector<int> subneighbors;
	sub.nodes(subnodes);
	double w;
	for( int x=0; x<(int)subnodes.size(); x++){
		int subnode = subnodes.at(x);
		sub.neighbors(subnode, subneighbors);
		for( int i=0; i<(int)subneighbors.size(); i++){
			for( int j=i+1; j<(int)subneighbors.size(); j++){
				if( this->has_edge( subneighbors.at(i), subneighbors.at(j), w ) ){
					status[c] = 1;
					status[subnode] = 1;
					status[subneighbors.at(i)] = 1;
					status[subneighbors.at(j)] = 1;
					return;
				}
			}
		}
	}
	return;
}

void Graph::reduce_to_quads( Graph& g ){
	std::vector<int> candidates;
	this->nodes(candidates);
	boost::unordered_map<int,int> status;
	for(int i=0; i<(int)candidates.size(); i++)
		status[candidates.at(i)]=0;

	if( this->verbose ){
		std::cout << "MESSAGE: Reducing to 4-cliques with " << candidates.size() << " candidates\n";
		std::cout.flush();
	}
	std::vector<int> neighbors;
	for(int i=0; i<(int)candidates.size(); i++){
		int c = candidates.at(i);
		this->neighbors(c, neighbors);
		if( ( (int)neighbors.size() ) < 3 ){
			continue;
		}
		else{
			if( status[c]==0 )
				evaluate_quad(status, c, neighbors);
		}
		if( this->verbose && i % 1000==0){
			std::cout << "MESSAGE: completed " << i + 1 << " of " << candidates.size() <<  " candidates\n";
			std::cout.flush();
		}
	}
	std::vector<int> quads;
	for(int i=0; i<(int)candidates.size(); i++){
		if( status[ candidates.at(i) ] == 1 )
			quads.push_back(candidates.at(i) );
	}
	this->subgraph(g, quads);
}

void Graph::clear_redundant( std::vector<std::vector<int>*>* results, std::vector<int>* target ){
    // clear any vector in results whose elements are a subset of the elements in vector target
    int i, j;
    bool clear;
    boost::unordered_map<int, int> in_target;
    for(i=0; i<(int)target->size(); i++)
        in_target[ target->at(i) ] = 1;
    for(i=0; i<(int)results->size(); i++){
        if( results->at(i)->size() < target->size() ){
            clear = true;
            for(j=0; j<(int)results->at(i)->size(); j++){
                if( in_target.find( results->at(i)->at(j) ) == in_target.end() ){
                    clear = false;
                    break;
                }
            }
            if(clear){
                results->at(i)->clear();
            }
        }
    }
}


void Graph::remove_vectors_by_size( std::vector<std::vector<int>*>* V, int min_size, bool verbose ){
    // Remove any vectors of size 0 from vector of vectors V.  
    std::vector<std::vector<int>*> O;
    for(int i=0; i<(int)V->size(); i++){
        if( (int)V->at(i)->size() < min_size )
            delete V->at(i);
        else
            O.push_back(V->at(i));
    }
    if( verbose ){
        std::cout << "MESSAGE: Removed " << V->size() - O.size() << " candidates below size threshold.\n";
        std::cout.flush();
    }
    V->clear();
    for(int i=0; i<(int)O.size(); i++)
        V->push_back( O.at(i) );
}


void Graph::semi_cliques(double completeness, std::vector<std::vector<int>*>* results, bool verbose){
    cliques(3,3,results, verbose);
}


void Graph::cliques( int min_k, int max_k, std::vector<std::vector<int>*>* results, bool verbose){
    // works but doesn't scale.  Do not bother to call with non-trivial graphs
    edges( results ); // initialize with all edges
    int n_this_round;
    std::vector<int> N;
    std::vector< std::vector<int>*>::iterator c_itr;
    Node root = Node(Node::ROOT);
    for(int i=0; i<(int)results->size(); i++){
        root.add_path(results->at(i));
    }

    for(int k=3;k<=max_k;k++){
        
        n_this_round = (int)results->size();
        if( verbose ){
            std::cout << "MESSAGE: Finding cliques of size " << k  << "\n";
            std::cout << "MESSAGE: Up to " << n_this_round  << " candidates\n";
            std::cout.flush();
        }
        for(int c=0; c<n_this_round; c++){
            if( (int)results->at(c)->size()==k-1 ){
                shared_neighbors(results->at(c), N);            
                for(int i=0; i<(int)N.size(); i++){
                    std::vector<int>* v = new std::vector<int>();
                    for(int j=0; j<(int)results->at(c)->size(); j++){
                        v->push_back( results->at(c)->at(j) );
                    }
                    v->push_back( N.at(i) );
                    if( root.has_path(v) ){
                        delete v;
                    }
                    else{
                        root.add_path(v);
                        results->push_back( v ); 
                    }
                }
            }
            if(verbose && c>1 && c%1000==0){
                std::cout << "MESSAGE: completed " << c << " of up to " << n_this_round  << " candidates\n";
                std::cout.flush();
            }
        }
    }
    remove_vectors_by_size(results, min_k, verbose);
}


void Graph::neighbors(int node_a, std::vector<int>& ids){
	ids.clear();
	if( adj.find(node_a) != adj.end() ){
		for(boost::unordered_map<int, double>::iterator node_itr=adj[node_a]->begin(); node_itr != adj[node_a]->end(); node_itr++){
			ids.push_back( node_itr->first );
		}
	}
}


void Graph::neighbors( int node_a, std::vector<int>& ids, std::vector<double>& weights){
	ids.clear();
	weights.clear();
	if( adj.find(node_a) != adj.end() ){
		for(boost::unordered_map<int, double>::iterator node_itr=adj[node_a]->begin(); node_itr != adj[node_a]->end(); node_itr++){
			ids.push_back( node_itr->first );
			weights.push_back( node_itr->second );
		}
	}
}


void Graph::read( std::string file_name ){
	FILE * pFile;
	pFile=fopen(file_name.c_str(), "r");

	//int line_no = 1;
	int a, b, out_f;
	float weight;
	while( pFile ){
		out_f = fscanf(pFile, "%d %d %f", &a, &b, &weight);
		if(out_f<0)
			break;
		this->add_edge(a, b, (double)weight);
	}
	fclose(pFile);
}

void Graph::read( std::string file_name, double min_abs_weight ){
	FILE * pFile;
	pFile = fopen(file_name.c_str(), "r");

	//int line_no = 1;
	int a, b, out_f;
	float weight;
	//float min_weight = (float) min_abs_weight;
	while( pFile ){
		out_f = fscanf(pFile, "%d %d %f", &a, &b, &weight);
		if(out_f<0)
			break;
		if( abs(weight) >= min_abs_weight )
			this->add_edge(a, b, (double)weight);
	}
	fclose(pFile);
}

void Graph::write( std::string fn_out, boost::unordered_map<int, std::string>& translate ){
    std::ofstream f_out(fn_out.c_str());
	if( !f_out.is_open() )
        throw new std::string( "Unable to open file for writing:"  + fn_out);
    int a,b;
    double weight;
    boost::unordered_map<int, double>::iterator edge_itr;
    for( graph_itr itr = adj.begin(); itr != adj.end(); itr++){
        a = itr->first;
        for(edge_itr = itr->second->begin(); edge_itr != itr->second->end(); edge_itr++){
            b = edge_itr->first;
            if( a<b ){
                if( this->has_edge(a,b,weight) )
                    f_out << translate[a] << "\t" << translate[b] << "\t" << weight << "\n";
                else{
                    std::cout << "ERROR: Missing edge with indices " << a << " and " << b << "\n";
                }
            }
        }
    }
    f_out.close();
}

Node::Node(int id){
	this->id=id;
}

Node::~Node(){
	boost::unordered_map<int, Node*>::iterator iter;
	for(iter=this->children.begin(); iter != this->children.end(); iter++)
		delete (*iter).second;
}

bool Node::has_path(std::vector<int>* path){
	// A path is a set of numbers representing the indices of genes in a rule.
	// The paths 1,2,3 and 2,3,1 are equivalent, so we must sort the test path 
	// to see if it is present in the tree.  This function is called to see 
	// if a newly enumerated rule has already been found.
	Node* c = this;
	std::vector<int> cp = std::vector<int>( path->size() );
	std::copy(path->begin(), path->end(), cp.begin());
	std::sort(cp.begin(), cp.end());

	for(int i=0; i<(int)cp.size(); i++){
		c = c->has_child(cp.at(i));
		if(c==NULL)
			return false;
	}
	return true;
}

void Node::get_kids(std::vector<int>* c){
	c->clear();
	for( boost::unordered_map<int, Node*>::iterator iter = this->children.begin(); iter != this->children.end(); iter++ ){
		c->push_back( (*iter).first );
	}
}

void Node::add_path(std::vector<int>* path){
	// A path is a set of numbers representing the indices of genes in a rule.
	// The paths 1,2,3 and 2,3,1 are equivalent, so we must sort the test path 
	// before inserting it so that we don't insert redundant paths.  This code 
	// does not check to see if the path being inserted already exists.
	std::vector<int> cp = std::vector<int>( path->size() );
	std::copy(path->begin(), path->end(), cp.begin());
	std::sort(cp.begin(), cp.end());
	Node* t = this;
	Node* c;
	for(int i=0; i<(int)cp.size(); i++){
		c = t->has_child(cp.at(i));
		if(c==NULL)
			c = t->add_child(cp.at(i));
		t = c;
	}
}


Node* Node::add_child(int id){
	Node* n = new Node(id);
	this->children[id] = n;
	return n;
}


void Node::print_kids(){
	for( boost::unordered_map<int, Node*>::iterator iter = this->children.begin(); iter != this->children.end(); iter++ ){
		std::cout << (*iter).first << " ";		
	}
	std::cout << "\n";
}


Node* Node::has_child(int id){
	boost::unordered_map<int, Node*>::iterator iter = this->children.find(id);
	if( iter==this->children.end() )
		return NULL;
	else{
		return (*iter).second;
	}
}

